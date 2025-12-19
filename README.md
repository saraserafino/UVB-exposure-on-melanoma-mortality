# UVB exposure on melanoma mortality

This repository contains the exam project of the course "Bayesian Statistics" @ UniTS, which consisted in the reproduction of some results of the paper [Multi-level modelling of geographically aggregated health data: a case study on malignant melanoma mortality and UV exposure in the European Community](https://pubmed.ncbi.nlm.nih.gov/9463848/) using the R package mlmRev that contains the [Mmmec](https://rdrr.io/cran/R2MLwiN/man/mmmec.html) dataset on malignant melanoma mortality in 9 nations of the European Community during 1971-1980.

The main goals of the project are to perform some explorative analysis also with a Bayesian approach, and to build a model for the number of male deaths taking into account the hierarchical data structure (counties within regions within nations) — both using stan_glmer relying on Hamiltonian Monter Carlo and writing a stan model as model A of the paper —. Moveover, results are compared through a self-contained PDF presentation.

The compiled R script is presented as an R markdown in html format.
