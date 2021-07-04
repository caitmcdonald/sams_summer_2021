Sam’s summer research plan
================

  - [Overview](#overview)
  - [Goals](#goals)
      - [Learning goals](#learning-goals)
          - [1. Bioinformatiics](#bioinformatiics)
          - [2. Lab techniques](#lab-techniques)
      - [Output goals](#output-goals)
          - [1. Collaborate on RNAseq analysis of exogenous/endogenous
            FeLV dataset,
            including:](#collaborate-on-rnaseq-analysis-of-exogenousendogenous-felv-dataset-including)
          - [2. Fully proficient in (and able to teach me\!) mammalian
            fibroblast cell culture
            protocols.](#fully-proficient-in-and-able-to-teach-me-mammalian-fibroblast-cell-culture-protocols.)
          - [3. Help develop project design for larger FeLV cell culture
            and genome sequencing
            project.](#help-develop-project-design-for-larger-felv-cell-culture-and-genome-sequencing-project.)
  - [Expectations](#expectations)
      - [My expectations](#my-expectations)
      - [Your expectations](#your-expectations)
  - [Timeline](#timeline)
  - [FeLV RNAseq dataset](#felv-rnaseq-dataset)
      - [Background](#background)
      - [Dataset](#dataset)
      - [Pre-processing that I have already
        done:](#pre-processing-that-i-have-already-done)
      - [Your job](#your-job)

**last updated:** 2021-07-04

-----

# Overview

This repository will (hopefully) serve as a central location for Sam’s
projects this summer. We’ll set up some mutual expectations, goals, and
a timeline here, and amend as necessary.

-----

# Goals

## Learning goals

### 1\. Bioinformatiics

  - Rstudio + RMarkdown
  - Version control: git and github integration
  - Basic shell programming
  - Basic cluster computing
  - \+/- Snakemake

### 2\. Lab techniques

  - Cell culture
  - Electronic lab notebook documentation

## Output goals

### 1\. Collaborate on RNAseq analysis of exogenous/endogenous FeLV dataset, including:

  - Troubleshooting analytical pipeline
  - Generating and writing up preliminary results
  - Generating associated figures

### 2\. Fully proficient in (and able to teach me\!) mammalian fibroblast cell culture protocols.

  - Obtain training from Raegan before she leaves in July.
  - Additional training as necessary from Laura, other lab mates.

### 3\. Help develop project design for larger FeLV cell culture and genome sequencing project.

-----

# Expectations

## My expectations

Much of what you’ll be doing is computer-based and independent. This is
great, but it can also be daunting and/or hard to keep to a schedule and
self-motivate. To help with this, let’s set up the following
expectations:

1.  Weekly in-person (or virtual) meetings to go over progress.
2.  Monthly group meetings with Sue to keep her updated.
3.  Email communication as necessary for troubleshooting.
      - Since we’re all working independently on a bunch of different
        projects, it can be hard to keep track of what’s going on. I’ll
        check in with you minimally weekly, and more frequently if you
        prefer. If you need help with anything, please feel free to ask
        at any point\!
4.  Time commitment: 40 hrs/week.
5.  That you learn some useful tools that will be beneficial to you for
    future projects and career goals.
6.  That you have fun\!

## Your expectations

\<insert your expectations here once you figure out git/github 👍 \>

# Timeline

| Due date | Goal                                                                                                                                                  | Date completed | Output notes/links |
| -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- | ------------------ |
| 5/14/21  | Meet with Raegan to go over summer cell culturing plan                                                                                                |                |                    |
| 5/14/21  | [Classes 1-5 of online data science class](https://nt246.github.io/NTRES6940-data-science/syllabus.html)                                              |                |                    |
| 5/14/21  | [Read background literature on FeLV, generate annotated bibliography](FeLV_references.Rmd)                                                            |                |                    |
| 5/14/21  | Research+choose best electronic lab notebook, whether github integration is possible                                                                  |                |                    |
| 5/21/21  | Start cell culture with Raegan                                                                                                                        |                |                    |
| 5/21/21  | Obtain Summit account                                                                                                                                 |                |                    |
| 5/21/21  | [Chapters 1-4 of Eric Anderson’s Bioinformatics Handbook](https://eriqande.github.io/eca-bioinf-handbook/essential-unixlinux-terminal-knowledge.html) |                |                    |
| 5/21/21  | [Classes 6-9 of online data science class](https://nt246.github.io/NTRES6940-data-science/syllabus.html)                                              |                |                    |

# FeLV RNAseq dataset

## Background

Following Elliott’s work, endogenous FeLV (enFeLV) can potentially
decrease the effects of exogenous FeLV (exFeLV) infection in domestic
cats. The goal of this RNAseq project is to identify whether enFeLV LTR
insertion site significantly alters transcription of genes involved in
response to exogenous FeLV infection.

Part of answering this question comes from comparing baseline expression
in different cell types. So, the questions you’ll be working on are:

**1. How does baseline expression differ between fibroblasts and
polymorphonuclear blood cells (PMBCs)?**

**2. How does this relate to enFeLV LTR insertion site?**

## Dataset

You will be working with the **uninfected** samples in [this larger
dataset](data/felv_metadata.tsv)

## Pre-processing that I have already done:

1.  Quality control (`fastqc` and `multiqc`)
    
    To assess sequencing quality, read quantity, etc.

2.  Read trimming (`trimgalore`)
    
    To remove low-quality reads and contamination. See this
    [script](scripts/trimgalore_test.sh).

3.  Read alignment and quantification (`STAR`)
    
    `STAR` will align reads to a genome of interest (in this case the
    cat genome), then quantify read counts to gene level. See this
    [script](scripts/starquant_test.sh).

4.  Restrict reads to only cat-puma orthologs (`OrthoFinder`)
    
    Because we have both cat and puma samples, we need to restrict our
    analysis to only cat genes with putative puma orthologs. I did so
    using this [script](scripts/orthofinder_test.sh).

## Your job

1.  First, run this [data exploration](scripts/sam_DGE_dataexplore.R)
    script and make sure you understand what it’s doing.
    
    Check out [this
    manual](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
    if you are unsure what a line of code is doing.

2.  Next, run this [differential gene expression
    analysis](scripts/sam_DGE_edgeR.R) script.
