FROM rocker/verse:3.5.1
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('flowCore'); biocLite('pamr'); install.packages('glmnet', repos = 'http://cran.us.r-project.org'); devtools::install_github('RcppCore/RcppEigen'); devtools::install_github('nolanlab/Rclusterpp'); biocLite('impute'); devtools::install_url('https://cran.r-project.org/src/contrib/Archive/samr/samr_2.0.tar.gz'); devtools::install_github('nolanlab/citrus'); install.packages('RANN');
devtools::install_github('JinmiaoChenLab/Rphenograph'); install.packages('Rtsne'); install.packages('hashmap'); install.packages('effsize'); devtools::install_github('Mamie/mwang')"
