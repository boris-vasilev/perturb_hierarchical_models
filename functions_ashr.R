library(ashr)
library(parallel)
library(tidyverse)
library(org.Hs.eg.db)


map_ensg_to_symbol <- function(perturbation_effect_df) {
  effects.ensg <- perturbation_effect_df$effect %>% unique
  
  effects.symbols <- mapIds(org.Hs.eg.db, keys = effects.ensg,
                            column = "SYMBOL", keytype = "ENSEMBL",
                            multiVals = "first")
  
  effects.symbols <- effects.symbols %>% na.omit
  
  perturbation_effect_df <- perturbation_effect_df %>%
    dplyr::rename(effect_ensg = "effect") %>%
    mutate(effect = effects.symbols[effect_ensg])
  
  return(perturbation_effect_df)
}

filter_efficient_perturbations <- function(perturbation_effect_df) {
  # Perturbations where self-effect is significantly downregulated (< 30% expression)
  strong_self_effects <- perturbation_effect_df %>%
    filter(effect == perturbation) %>%
    mutate(percent_change = (2^avg_log2FC - 1)) %>%
    filter(percent_change < -0.7 & p_val_adj < 0.05) %>% 
    pull(perturbation)
  valid_perturbations <- strong_self_effects

#  # Perturbations where self-effect is missing (not detected at all)
#  all_perturbed <- perturbation_effect_df %>% pull(perturbation) %>% unique
#  detected_self_effects <- perturbation_effect_df %>%
#    filter(effect == perturbation) %>%
#    pull(perturbation)
#  
#  undetected_self_effects <- setdiff(all_perturbed, detected_self_effects)
#  
#  # Combine both sets
#  valid_perturbations <- union(strong_self_effects, undetected_self_effects)
#  
  perturbation_effect_df <- perturbation_effect_df %>%
    filter(perturbation %in% valid_perturbations)
  
  return(perturbation_effect_df)
}

fit_ash_normal <- function(l2fc, l2fc_se,min_l2fc,max_l2fc, n_sample) {
  ash_model <- ash(
    betahat = l2fc,
    sebetahat = l2fc_se,
    mixcompdist = "normal"
    #pointmass = TRUE
  )
  return(list(fit = ash_model))
}


fit_ash <- function(l2fc,l2fc_se,min_l2fc,max_l2fc,n_sample) {
  #Suppress warning about biased lfsr estimates due to not using null biased prior
#  withCallingHandlers({
#    ash_model = ash(betahat = l2fc,
#                    sebetahat = l2fc_se,
#                    mixcompdist = "halfuniform",
#                    outputlevel = 3,
#                    grange = c(min_l2fc,max_l2fc),
#                    prior = "uniform")  }, warning = function(w) {
#    if (grepl("nullbiased", conditionMessage(w))) {
#      invokeRestart("muffleWarning")  # Suppresses this warning only
#    }
#  })

  ash_model = ash(betahat=l2fc, sebetahat = l2fc_se, prior = "uniform")
  
  weights =  ash_model$fitted_g$pi
  uniform_a = ash_model$fitted_g$a
  uniform_b = ash_model$fitted_g$b


  component = sample(c(1:length(weights)),size = n_sample, prob = weights, replace = TRUE)
  samples = runif(n_sample, min = uniform_a[component], max = uniform_b[component])


  return(list(fit = ash_model,
              samples = samples))
}

apply_ash_single <- function(data) {
  min_lfc <- min(data$log2FoldChange)
  max_lfc <- max(data$log2FoldChange)
  n_sample <- nrow(data)
  
  return(fit_ash_normal(data$log2FoldChange, data$lfcSE, min_lfc, max_lfc, n_sample))
}


#apply_ash <- function(DEG_files, cores=4, significant=TRUE, bulk=F, efficient=TRUE) {
#  print("Applying ash")
#  ash.results <- mclapply(DEG_files, function(file) {
#    # Read the file
#    DEGs <- read.table(file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#    # Extract perturbed gene from filename
#    perturbed_gene <- basename(file) %>% sub("\\.tsv$", "", .)
#    
#    if(!is.numeric(DEGs$lfcSE)) return(NULL)
#    model <- apply_ash_single(DEGs)
#    results <- model$fit$result %>%
#      mutate(perturbation = perturbed_gene,
#	     effect = rownames(DEGs)) %>%
#      rename(avg_log2FC = PosteriorMean)
#    if(significant) {
#      results <- results %>% filter(svalue < 0.05)
#    }
#    return(results)
#  }, mc.cores = cores)
#  ash.results <- bind_rows(ash.results) 
#  return(ash.results)
#}

#extract_DEGs <- function(DEG_files, cores=4, significant=TRUE, bulk=F, efficient=TRUE) {
#  print("Extracting DEGs")
#  
#  summary_list <- list(total_perturbs = length(DEG_files))
#  
#  significant_perturb_effects <- mclapply(DEG_files, function(file) {
#    DEGs <- read.table(file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#    perturbed_gene <- basename(file) %>% sub("\\.tsv$", "", .)
#    if(bulk) {
#      DEGs <- DEGs %>% dplyr::rename(avg_log2FC="log2FoldChange", p_val_adj = "padj")
#    }
#    results <- DEGs %>%
#      mutate(perturbation = perturbed_gene) %>%
#      rownames_to_column(var = "effect")
#    rm(DEGs)
#    gc()
#    return(results)
#  }, mc.cores = cores, mc.preschedule = FALSE)
#  significant_perturb_effects <- bind_rows(significant_perturb_effects)
#  
#  significant_perturb_effects <- map_ensg_to_symbol(significant_perturb_effects)
#  
#  summary_list$total_pairs <- nrow(significant_perturb_effects)
#  
#  if(efficient) {
#    print("Filtering efficient perturbations")
#    significant_perturb_effects <- filter_efficient_perturbations(significant_perturb_effects)
#    summary_list$efficient_perturbs <- significant_perturb_effects %>%
#      pull(perturbation) %>%
#      unique %>%
#      length
#    summary_list$efficient_pairs <- nrow(significant_perturb_effects)
#  }
#  
#  if(significant) {
#    print("Filtering significant DEGs")
#    significant_perturb_effects <- significant_perturb_effects %>%
#      filter(p_val_adj < 0.05)
#    summary_list$significant_pairs <- nrow(significant_perturb_effects)
#  }
#  
#  return (
#    list(
#      perturbation_effect_df = significant_perturb_effects,
#      filter_summary = summary_list
#    )
#  )
#}


extract_DEGs <- function(DEG_files, cores=4, significant=TRUE, bulk=F, efficient=TRUE, ash=F) {
  print("Extracting DEGs")
  
  summary_list <- list(total_perturbs = length(DEG_files))
  
  significant_perturb_effects <- mclapply(DEG_files, function(file) {
    tryCatch({
      message("Reading file: ", file)
      DEGs <- read.table(file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
      
      required_cols <- c("gene", "log2FoldChange", "padj")
      
      if(!all(required_cols %in% colnames(DEGs))) {
        warning(paste("Skipping", file, "- missing required columns"))
        return(NULL)
      }

      perturbed_gene <- basename(file) %>% sub("\\.tsv$", "", .)

      if(ash) {
        if(!is.numeric(DEGs$lfcSE)) return(NULL)
        model <- apply_ash_single(DEGs)

        DEGs <- DEGs %>%
          mutate(log2FoldChange= model$fit$result$PosteriorMean,
                padj= model$fit$result$svalue)
      }

      if(bulk) {
        DEGs <- DEGs %>% dplyr::rename(avg_log2FC="log2FoldChange", p_val_adj = "padj")
      }
      results <- DEGs %>%
        mutate(perturbation = perturbed_gene) %>%
        rownames_to_column(var = "effect")
      return(results)
    }, error = function(e) {
      warning("Error in file ", file, ": ", conditionMessage(e))
      return(NULL)
    })
  }
  , mc.cores = cores, mc.preschedule = FALSE)
  significant_perturb_effects <- bind_rows(significant_perturb_effects)
  significant_perturb_effects <- map_ensg_to_symbol(significant_perturb_effects)
  
  summary_list$total_pairs <- nrow(significant_perturb_effects)
  
  if(efficient) {
    print("Filtering efficient perturbations")
    significant_perturb_effects <- filter_efficient_perturbations(significant_perturb_effects)
    summary_list$efficient_perturbs <- significant_perturb_effects %>%
      pull(perturbation) %>%
      unique %>%
      length
    summary_list$efficient_pairs <- nrow(significant_perturb_effects)
  }
  
  if(significant) {
    print("Filtering significant DEGs")

    significant_perturb_effects <- significant_perturb_effects %>%
      filter(p_val_adj < 0.05)
    summary_list$significant_pairs <- nrow(significant_perturb_effects)
  }
  
  return (
    list(
      perturbation_effect_df = significant_perturb_effects,
      filter_summary = summary_list
    )
  )
}


