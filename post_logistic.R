library(tidyverse)
library(glue)
library(here)
library(brms)
library(argparse)
library(broom)

supported_models <- c("varying_intercept_varying_slope",
                      "no_pooling_intercept_varying_slope",
                      "varying_intercept_fixed_slope",
		                  "vivs_student")

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)
parser$add_argument("--model", type = "character", help = "Model type", required = TRUE, choices = supported_models)

parser$add_argument("--ash", action = "store_true", help="Run adaptive shrinkage (ashr)")
# parser$add_argument("--bulk", action = "store_true", help="Pseudobulk")
parser$add_argument("--efficient", action = "store_true", help="Efficient perturbations >70% only")

args <- parser$parse_args()
cells <- args$cells
model <- args$model
efficient <- args$efficient

ash <- if(args$ash) "ash_" else ""
# bulk <- if(args$bulk) "bulk_" else ""
eff <- if(args$efficient) "_eff" else ""

print(glue("Cells: {cells}            Model: {model}           Efficient: {efficient} "))
dat <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat.csv") %>%
  filter(screen == cells, perturb_eff >= 0.7) %>%
  select(perturb, effect, x, y) %>%
  group_by(perturb) %>%
  filter(any(x == 1 & y == 1)) %>%
  ungroup
fit <- readRDS(here(glue("/rds/project/rds-csoP2nj6Y6Y/biv22/models/{cells}/{ash}logit_{model}{eff}.rds")))

png(here(glue("plots/{cells}/traceplot_logit_{model}.png")), width= 1000, height = 1000)
plot(fit)
dev.off()

partially_pooled_params <-
  # with this line we select each of the perturbation posterior mean (i.e., Estimate)
  # for both `Intercept` and `x - perturb log2FC`
  coef(fit)$perturb[ , 1, 1:2] %>%
  data.frame() %>%              # convert the two vectors to a data frame
  rename(Slope = x) %>%
  rownames_to_column("perturb")

un_pooled_params <- dat %>%
  group_by(perturb) %>%
  do(tidy(glm(y ~ x, data = ., family="binomial"))) %>%  # fit glm per group and tidy the output
  filter(term %in% c("(Intercept)", "x")) %>%  # keep only intercept and slope
  select(perturb, term, estimate) %>%
  spread(key = term, value = estimate) %>%
  rename(Intercept = `(Intercept)`, Slope = x)

# here we combine the partially-pooled and unpooled means into a single data object, 
# which will make plotting easier.
params <-
  # `bind_rows()` will stack the second tibble below the first
  bind_rows(partially_pooled_params, un_pooled_params) %>%
  # index whether the estimates are pooled
  mutate(pooled = rep(c("partially", "not"), each = nrow(.)/2)) 

p1 <-
  ggplot(data = params, aes(x = Intercept, y = Slope)) +
  stat_ellipse(geom = "polygon", type = "norm", level = 1/10, linewidth = 0, alpha = 1/10) +
  stat_ellipse(geom = "polygon", type = "norm", level = 2/10, linewidth = 0, alpha = 1/10) +
  stat_ellipse(geom = "polygon", type = "norm", level = 3/10, linewidth = 0, alpha = 1/10) +
  stat_ellipse(geom = "polygon", type = "norm", level = 4/10, linewidth = 0, alpha = 1/10) +
  stat_ellipse(geom = "polygon", type = "norm", level = 5/10, linewidth = 0, alpha = 1/10) +
  stat_ellipse(geom = "polygon", type = "norm", level = 6/10, linewidth = 0, alpha = 1/10) +
  stat_ellipse(geom = "polygon", type = "norm", level = 7/10, linewidth = 0, alpha = 1/10) +
  stat_ellipse(geom = "polygon", type = "norm", level = 8/10, linewidth = 0, alpha = 1/10) +
  stat_ellipse(geom = "polygon", type = "norm", level = 9/10, linewidth = 0, alpha = 1/10) +
  stat_ellipse(geom = "polygon", type = "norm", level = .99,  linewidth = 0, alpha = 1/10) +
  geom_point(aes(group = perturb, color = pooled), alpha = 0.2, size = 3) +
  geom_line(aes(group = perturb), linewidth = 1/4, alpha = 0.2) +
  scale_color_manual("Pooled?", values = c("blue", "red")) +
  coord_cartesian(xlim = range(params$Intercept),
                  ylim = range(params$Slope))

png(here(glue("plots/{cells}/logit_parameter_shrinkage_{model}.png")), width = 1080, height = 1080)
p1 + theme_bw() + theme(text = element_text(size = 20))
dev.off()


# Extract intercepts and slopes with their 95% credible intervals
partially_pooled_intercepts.confidence <- coef(fit)$perturb[,,1] %>%
	data.frame %>%
	select(Estimate, Q2.5, Q97.5) %>%
	rownames_to_column("perturb")

partially_pooled_slopes.confidence <- coef(fit)$perturb[,,2] %>%
	data.frame %>%
	select(Estimate, Q2.5, Q97.5) %>%
	rownames_to_column("perturb")

p <- ggplot(partially_pooled_intercepts.confidence, aes(x = perturb, y = Estimate)) +
  geom_point(size = 3) +  # plot the estimates as points
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2) +  # plot credible intervals
  labs(x = "Perturbation", y = "Estimate", title = "Intercepts with Credible Intervals (partial pooling)") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")
  #coord_flip()

png(here(glue("plots/{cells}/logit_{model}_intercepts.png")), width = 3080, height = 1080)
p + theme_void() + theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

alpha_pop <- fixef(fit)["Intercept", "Estimate"]
partially_pooled_intercepts.significant <- partially_pooled_intercepts.confidence %>%
         filter((Q2.5 > 0 | Q97.5 < 0) & !(Q2.5 < alpha_pop & Q97.5 > alpha_pop)) 

p <- ggplot(partially_pooled_intercepts.significant, aes(x = perturb, y = Estimate)) +
  geom_point(size = 2) +  # plot the estimates as points
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2) +  # plot credible intervals
  labs(x = "Perturbation", y = "Estimate", title = "Significant Intercepts with Credible Intervals (partial pooling)") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")
  #coord_flip()

png(here(glue("plots/{cells}/logit_{model}_intercepts_significant.png")), width = 580, height = 580)
p + theme_bw() + theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

p <- ggplot(partially_pooled_slopes.confidence, aes(x = perturb, y = Estimate)) +
  geom_point(size = 3) +  # plot the estimates as points
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2) +  # plot credible intervals
  labs(x = "Perturbation", y = "Estimate", title = "Slopes with Credible Intervals (partial pooling)") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")
  #coord_flip()

png(here(glue("plots/{cells}/logit_{model}_slopes.png")), width = 3080, height = 1080)
p + theme_void() + theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

partially_pooled_slopes.significant <- partially_pooled_slopes.confidence %>%
         filter(Q2.5 > 0 | Q97.5 < 0) 

p <- ggplot(partially_pooled_slopes.significant, aes(x = perturb, y = Estimate)) +
  geom_point(size = 2) +  # plot the estimates as points
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2) +  # plot credible intervals
  labs(x = "Perturbation", y = "Estimate", title = "Significant Slopes with Credible Intervals (partial pooling)") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")
  #coord_flip()

png(here(glue("plots/{cells}/logit_{model}_slopes_significant.png")), width = 580, height = 580)
p + theme_bw() + theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

