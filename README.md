## MBM: Multifaceted Biodiversity Modelling

R Package for simultaneously modelling multiple facets of biodiversity (α- and β-diversity, taxonomic, functional, and phylogenetic).

### Installation

Before installing MBM, it is necessary to have installations of [R](https://cran.r-project.org/), [Python](https://www.python.org/) and [GPy](https://sheffieldml.github.io/GPy/). Please install these packages and verify they are working before installing MBM.

Within R, MBM can be installed directly from github using the `devtools` package.

    # uncomment the line below if devtools is not yet installed
    # install.packages('devtools') 
    
    library('devtools')
    install_github('mtalluto/mbm')
    
    library('mbm')
    ?mbm 

Alternatively, if you would like additional features (principally, vignettes), you can clone or download the project and build and install it (from the command line):

    git clone https://github.com/mtalluto/mbm.git
    Rscript -e "devtools::build('mbm')"
    R CMD INSTALL mbm_i.j.k.tar.gz
    
Replacing the `i.j.k` with the version number you have downloaded; this file will have been created by devtools in the second step.