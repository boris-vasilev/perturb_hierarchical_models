nextflow.enable.dsl=2

params.perturbList = '/Users/borisvasilev/PhD/perturb_hierarchical_models/cross_screen_models/perturblist'
params.runSummary = '/Users/borisvasilev/PhD/perturb_hierarchical_models/cross_screen_models/run_summary.txt'

process FIT_CROSS_SCREEN_MODELS {
  time '2h'
  cpus 4
  errorStrategy 'ignore'

  input:
    val perturb

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

process RUN_SUMMARY {
    publishDir ".", mode: "copy"

    input:
      tuple val(perturb), path(status)

    output:
      path "run_summary.txt"

    script:
    """
    if grep -q FAIL ${status}; then
        echo "${perturb}" >> run_summary.txt
    fi
    """
}


workflow {
  perturbList = Channel.fromPath(params.perturbList).splitText()

  results = FIT_CROSS_SCREEN_MODELS(perturbList)

  RUN_SUMMARY(results)
}
