#!/usr/bin/Rscript
source("main.R")
library(testthat)

csv <- paste0("data/verse_counts.tsv")

describe("Data Loading and Reduction", {
  it("should load and trim test data correctly", {
    test_df <- load_n_trim(csv)
    expect_equal(names(test_df), c("vP0_1", "vP0_2", "vAd_1", "vAd_2"))
    expect_equal(dim(test_df), c(55416, 4))
    expect_equal(class(test_df), "data.frame")
  })
})

describe("DESeq2 Functionality", {
  it("should return correct results", {
    load("data/mock_counts_df.RData") # loads the counts_df object into env
    coldata <- data.frame(condition = rep(c("day4", "day7"), each=2),
                          type="paired-end")
    row.names(coldata) <- c("vP4_1", "vP4_2", "vP7_1", "vP7_2")
    expect_warning(deseq <- run_deseq(counts_df, coldata, 10, "condition_day7_vs_day4"))
    expect_equal(dim(deseq), c(19127, 6))
    expect_equal(class(deseq)[1], "DESeqResults")
    expect_equal(c("pvalue", "padj") %in% names(deseq), c(TRUE, TRUE))
  })
})

describe("edgeR Functionality", {
  it("should return correct results", {
    load("data/mock_counts_df.RData")
    group <- factor(rep(c(1,2), each=2))
    edger_res <- run_edger(counts_df, group)
    expect_equal(names(edger_res), c("logFC", "logCPM", "PValue"))
    expect_equal(dim(edger_res), c(15026, 3))
  })
})

describe("Limma + Voom Integration", {
  it("should work correctly together", {
    load("data/mock_counts_df.RData")
    group <- factor(rep(c(1,2), each=2))
    design <- data.frame(day4=1, day4vsday7=c(0, 0, 1, 1))
    row.names(design) <- c("vP4_1", "vP4_2", "vP7_1", "vP7_2")
    expect_warning(voom_res <- run_limma(counts_df, design, group))
    expect_equal(dim(voom_res), c(15026, 6))
    expect_equal(names(voom_res), c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
  })
})
