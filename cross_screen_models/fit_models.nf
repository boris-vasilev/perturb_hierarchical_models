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


/*
 * PROCESS 2 — Fit models
 */
process FIT_CROSS_SCREEN_MODELS {
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
  perturb_model_comparison.R --perturb ${perturb}
  exit_code=\$?

  if [[ \$exit_code -eq 0 ]]; then
      echo SUCCESS > status.txt
  else
      echo FAIL > status.txt
  fi

  exit \$exit_code
  """
}


/*
 * PROCESS 3 — Summarise results (runs once)
 */
process RUN_SUMMARY {
  publishDir ".", mode: "copy"

  input:
    tuple val(perturb), path(status)

  output:
    path "run_summary.txt"

  when:
    status.text.contains('FAIL')

  script:
  """
  echo "${perturb}" > run_summary.txt
  """
}

workflow {
  perturbList = Channel.fromPath(params.perturbList).splitText()
  pre = PRECOMPILE_MODELS()

  results = FIT_CROSS_SCREEN_MODELS(perturbList, pre)

  RUN_SUMMARY(results)
    .collectFile(name: "run_summary.txt", mode: 'append')
}

