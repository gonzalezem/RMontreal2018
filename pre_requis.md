# Pré-requis
Pour cette session, vous aurez besoin de:

* Un ordinateur portable (windows, mac ou linux)
* Rstudio Desktop 1.1.453 - [Download RStudio](https://www.rstudio.com/products/rstudio/download/#download)
* R 3.5.0 - [Download R](https://cran.rstudio.com/)
* 3 libraries R: [dada2](https://benjjneb.github.io/dada2/dada-installation.html), [Phyloseq](http://joey711.github.io/phyloseq/install.html) et [ggplot2](https://cran.r-project.org/web/packages/ggplot2/README.html) (voir ci-dessous pour l'installation)
* Un dossier zippé contenant les sequences ADN  - [Download folder](https://drive.google.com/a/computationalgenomics.ca/file/d/1yw8O67HExy4mHuJo3EJRlHyTKUU_QTXe/view?usp=sharing)

## Pour installer les 3 libraries R
1. Ouvrir RStudio
2. Dans la console (partie en bas gauche), taper les commandes:
```
source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
biocLite('phyloseq')
biocLite("ggplot2")
```
