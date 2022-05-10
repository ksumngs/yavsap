#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { cowsay } from './lib/cowsay.nf'
include { yavsap_logo } from './lib/logo.nf'

cowsay(yavsap_logo())

if (params.help) {
    log.info(
        """\
        YAVSAP (Yet Another Viral Subspecies Analysis Pipeline) - Intra-sample viral
        population analysis

        Usage:

            nextflow run ksumngs/yavsap

        Options:

            --input             Relative or absolute path to directory containing
                                gzipped fastq files
                                    type: path, default: .

            --platform          Type of reads to process. Options are 'illumina' and
                                'nanopore'
                                    type: string, default: none

            --genome            NCBI accession number of the reference genome to
                                align reads against
                                    type: string, default: 'NC_001437.1'

            --kraken2_db        Kraken2-compatible database for classifying reads
                                    type: path, default: none

            --keep_taxid        Space-separated list of NCBI taxids to keep and
                                process after classifying
                                    type: string, default: '0 10239'

            --outdir            Directory in which to place results
                                    type: path, default: ./results

            --help              Print this message and exit

        For more information on usage and parameters, visit the website at
            https://ksumngs.github.io/yavsap
    """.stripIndent()
    )

    exit 0
}

log.info(
    """\
    Input folder:           ${params.input}
    Sequencing platform:    ${params.platform}
    Reference genome:       ${params.genome}
    Kraken2 Database:       ${params.kraken2_db}
    Taxonomic Ids:          '${params.keep_taxid}'
    Output folder           ${params.outdir}
    Diagnostics folder:     ${params.tracedir}
    """.stripIndent()
)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/yavsap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/yavsap
    Website: https://nf-co.re/yavsap
    Slack  : https://nfcore.slack.com/channels/yavsap
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { YAVSAP } from './workflows/yavsap'

//
// WORKFLOW: Run main nf-core/yavsap analysis pipeline
//
workflow NFCORE_YAVSAP {
    YAVSAP ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_YAVSAP ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
