nextflow.enable.dsl=2

params.perturbList = '/Users/borisvasilev/PhD/perturb_hierarchical_models/cross_screen_models/perturblist'
params.runSummary  = '/Users/borisvasilev/PhD/perturb_hierarchical_models/cross_screen_models/run_summary.txt'


/*
 * PROCESS 1 — Precompile models
 */
process PRECOMPILE_MODELS {
  time '20m'
  cpus 4

  output:
    path "precompile_done.flag"

  script:
  """
  precompile_models.R
  touch precompile_done.flag
  """
}

process FIT_NO_ME {
  time '10m'
  cpus 4
  errorStrategy 'ignore'

  input:
    val perturb
    path precompile_flag

  output:
    tuple val(perturb), path("status.txt")

  script:
  """
  fit_no_me.R --perturb ${perturb}
  """
}

process FIT_ME_J {
  time '20m'
  cpus 4
  errorStrategy 'ignore'

  input:
    val perturb
    path precompile_flag

  output:
    tuple val(perturb), path("status.txt")

  script:
  """
  fit_me_j.R --perturb ${perturb}
  """
}

process FIT_FULL_ME {
  time '2h'
  cpus 4
  errorStrategy 'ignore'

  input:
    val perturb
    path precompile_flag

  output:
    tuple val(perturb), path("status.txt")

  script:
  """
  fit_full_me.R --perturb ${perturb}
  """
}

workflow {
  perturbList = Channel.fromPath(params.perturbList).splitText()
  pre = PRECOMPILE_MODELS()

  results_no_me = FIT_NO_ME(perturbList, pre)
  results_me_j = FIT_ME_J(perturbList, pre)
  results_me_both = FIT_FULL_ME(perturbList, pre)

}

