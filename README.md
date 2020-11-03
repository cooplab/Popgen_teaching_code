# Popgen teaching code

Code to illustrate various ideas in population genetics. This code accompanies the [Population and Quantitative Genetics notes](https://github.com/cooplab/popgen-notes/releases).

Feel free to reuse/repurpose this code.

```R
pkgs <- c("RColorBrewer","shiny","kinship2")
dl_pkgs <- subset(pkgs,!pkgs %in% rownames(installed.packages()))
if(length(dl_pkgs)!=0){
  for(i in dl_pkgs) install.packages(i)
}
library(shiny)

runGitHub("cooplab/Popgen_teaching_code/",subdir = "Pheno_selection/")
runGitHub("cooplab/Popgen_teaching_code/",subdir = "Simulate_drift")
runGitHub("cooplab/Popgen_teaching_code/",subdir ="Simple_coalescent")
runGitHub("cooplab/Popgen_teaching_code/",subdir ="IBD_pedigree")
```
