#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(vroom)
    library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: reformat_gwas.R <input_file> <output_file>")
}
in_file  <- args[1]
out_file <- args[2]

message("Reading: ", in_file)
gwas_raw <- vroom::vroom(in_file, show_col_types = FALSE)

gwas_formatted <- gwas_raw %>%
    select(rsid, chromosome, base_pair_location, effect_allele, other_allele, p_value, beta) %>%
    rename(
        SNP  = rsid,
        CHR  = chromosome,
        POS  = base_pair_location,
        A1   = effect_allele,
        A2   = other_allele,
        P    = p_value,
        BETA = beta
    ) %>%
    mutate(
        BETA = as.numeric(BETA),
        A1   = toupper(A1),
        A2   = toupper(A2),
        CHR  = as.integer(CHR),
        POS  = as.integer(POS)
    )

message("Writing: ", out_file)
vroom::vroom_write(gwas_formatted, out_file)
message("Done")
