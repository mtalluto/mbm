## MBM: Multifaceted Biodiversity Modelling

R Package for simultaneously modelling multiple facets of biodiversity (α- and β-diversity, taxonomic, functional, and phylogenetic).

### Installation

Before installing MBM, it is necessary to have installations of [R](https://cran.r-project.org/), [Python](https://www.python.org/) and [GPy](https://sheffieldml.github.io/GPy/). Please install these packages and verify they are working before installing MBM. 

Additionally, some functionality (namely, `svgp` for fitting large models) requires the `climin` package in Python. Thus, it is highly recommended to install this package when instally GPy. Note that `climin pre-0.1` [is the only supported version](https://github.com/SheffieldML/GPy/issues/327). Thus it is necessary to [follow these instructions](https://climin.readthedocs.io/en/latest/installation.html) to install climin from Github.



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

### Macos Installation Notes
The following procedure has been tested on a clean mac running macos 10.13 using the system Python (version 2.7). From the command line:

    sudo easy_install pip
    sudo pip install -U --ignore-installed numpy scipy
    sudo pip install GPy tornado climin
    
From there, assuming everything proceeds without errors, you can install and test MBM.