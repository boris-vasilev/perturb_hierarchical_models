nextflow.enable.dsl=2

params.perturbList = '/Users/borisvasilev/PhD/perturb_hierarchical_models/cross_screen_models/perturblist'

process FIT_CROSS_SCREEN_MODELS {
  time '2h'
  cpus 4

  input:
    val perturb
  
  script:
  """
    perturb_model_comparison.R --perturb ${perturb}
  """
}

workflow {
  perturbList = Channel.fromPath(params.perturbList).splitText()
  FIT_CROSS_SCREEN_MODELS(perturbList)
}