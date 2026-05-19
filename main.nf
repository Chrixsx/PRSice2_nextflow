#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------------------------------------------------------------------
// Parameters — override on the command line with --param value
// ---------------------------------------------------------------------------
params.traits_file = "/vast/projects/Epilepsy_Metabolites/scripts/PRSice2/prsice2_nextflow/traits.txt"
params.suffix      = "_Karjalainen2024_reformatted_prsice2.tsv.gz"
params.cohort      = "ERC_Dec2025"

// Step toggle flags
params.run_reformat = true
params.run_prsice   = false

// Paths
params.gwas_dir  = "/vast/projects/Epilepsy_Metabolites/data/Karjalainen_2024/original"
params.sumstats  = "/vast/scratch/users/le.c/Kajarlainen2024/reformatted_prsice2"

params.app       = "/vast/projects/Epilepsy_Metabolites/scripts/PRSice2/app"
params.plink     = "/vast/projects/Epilepsy_Metabolites/scripts/prsCS/prs_metab_GSA_Dec2025/data/imputted_GSA_SNP_Dec2025/EUR_only_gsa_QC3_maf_0.001"
params.outdir    = "/vast/projects/Epilepsy_Metabolites/scripts/PRSice2/prs_metab_GSA_Dec2025/2.output"

// ---------------------------------------------------------------------------
// Step 1: Reformat raw Karjalainen GWAS → PRSice2 column names
// ---------------------------------------------------------------------------
process REFORMAT_GWAS {
    tag "${trait}"

    input:
    val trait

    output:
    val trait

    script:
    def in_file  = "${params.gwas_dir}/${trait}.h.tsv.gz"
    def out_file = "${params.sumstats}/${trait}${params.suffix}"
    """
    module load R
    mkdir -p ${params.sumstats}
    Rscript ${projectDir}/bin/reformat_gwas.R "${in_file}" "${out_file}"
    """
}


// ---------------------------------------------------------------------------
// Step 2: Polygenic risk scores with PRSice2
// ---------------------------------------------------------------------------
process RUN_PRSICE {
    tag "${trait}"

    // All raw outputs → 1_raw/<trait>/; scores also copied to 2_score/
    publishDir "${params.outdir}/1_raw/${trait}", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/2_score",        mode: 'copy', overwrite: true, pattern: "*.all_score"

    input:
    val trait

    output:
    path "*"

    script:
    def file_name  = "${trait}${params.suffix}"
    def out_prefix = "${trait}_${params.cohort}"
    // Check published location so a second run / -resume can reuse a prior .valid file
    def valid_pub  = "${params.outdir}/1_raw/${trait}/${out_prefix}.valid"
    """
    module load R

    _run_prsice() {
        Rscript ${params.app}/PRSice.R \\
            --seed 2026 \\
            --prsice ${params.app}/PRSice_linux \\
            --base "${params.sumstats}/${file_name}" \\
            --target ${params.plink} \\
            --binary-target T \\
            --snp SNP --chr CHR --bp POS \\
            --A1 A1 --A2 A2 \\
            --pvalue P --stat BETA \\
            --out "${out_prefix}" \\
            --bar-levels 0.5,0.05,1e-5,5e-8 \\
            --fastscore \\
            --no-regress \\
            --no-full \\
            --model add \\
            --missing MEAN_IMPUTE \\
            --score avg \\
            --clump-kb 250 \\
            --clump-p 1.0 \\
            --clump-r2 0.1 \\
            --upper 1 \\
            --print-snp \\
            "\$@"
    }

    if [[ -f "${valid_pub}" ]]; then
        echo "Found published .valid file — running with --extract"
        _run_prsice --extract "${valid_pub}"
    else
        _run_prsice || {
            local_valid="${out_prefix}.valid"
            if [[ -f "\${local_valid}" ]]; then
                echo "Duplicate RSIDs detected — rerunning with --extract"
                _run_prsice --extract "\${local_valid}"
            else
                echo "PRSice2 failed and no .valid file was generated" >&2
                exit 1
            fi
        }
    fi
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    traits_ch = Channel
        .fromPath(params.traits_file)
        .splitText()
        .map  { it.trim() }
        .filter { it && !it.startsWith('#') }

    // Step 1: reformat (optional)
    if (params.run_reformat) {
        REFORMAT_GWAS(traits_ch)
        prsice_input_ch = REFORMAT_GWAS.out
    } else {
        prsice_input_ch = traits_ch
    }

    // Step 2: PRS (optional)
    if (params.run_prsice) {
        RUN_PRSICE(prsice_input_ch)
    }
}
