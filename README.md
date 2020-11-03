# Popgen teaching code

Code to illustrate various ideas in population genetics. This code accompanies the [Population and Quantitative Genetics notes](https://github.com/cooplab/popgen-notes/releases).

Feel free to reuse/repurpose this code.

To run these shinys locally using the following code. 
```R
##download the required packages. Note kinship2 is needed only for drawing the pedigrees in IBD_pedigree.
pkgs <- c("RColorBrewer","shiny","kinship2")
dl_pkgs <- subset(pkgs,!pkgs %in% rownames(installed.packages()))
if(length(dl_pkgs)!=0){
  for(i in dl_pkgs) install.packages(i)
}
library(shiny)

##choose one of the following
runGitHub("cooplab/Popgen_teaching_code/",subdir = "Pheno_selection/")
runGitHub("cooplab/Popgen_teaching_code/",subdir = "Simulate_drift")
runGitHub("cooplab/Popgen_teaching_code/",subdir ="Simple_coalescent")
runGitHub("cooplab/Popgen_teaching_code/",subdir ="IBD_pedigree")
```
This software is free for use or modification under the GNU General Public License 3.0 (https://opensource.org/licenses/GPL-3.0)
