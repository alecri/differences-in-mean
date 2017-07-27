###  Dose-response meta-analysis of differences in means

**[Alessio Crippa](http://alecri.github.io)<sup>1</sup>, and [Nicola Orsini](http://nicolaorsini.altervista.org)<sup>1</sup>**

_<sup>1</sup>Department of Public Health Sciences, Karolinska Institutet, Stockholm, Sweden_

---

The repository contains material to reproduce results (numbers in text, tables, and figures) of the article "Dose-response meta-analysis of differences in means" (BMC medical research methodology 2016).

Note: version 2.0.0 of the R package `dosresmeta` is required. It can be downloaded from http://github.com/alecri/dosresmeta or installed from R using the following command lines:

    install.packages("devtools")
    devtools::install_github("alecri/dosresmeta")

This repository contains:
- `diff_in_means.R`, an R script to reproduce the main findings.
- `no_lapply.R`, same as diff_in_means.R but without the use of lapply and Map functions (use of loops)

---

Reference: Crippa A, Orsini N. _Dose-response meta-analysis of differences in means_. BMC Medical Research Methodology. 2016 Aug 2;16(1):91. [download](https://www.researchgate.net/publication/305804878_Dose-response_meta-analysis_of_differences_in_means)
