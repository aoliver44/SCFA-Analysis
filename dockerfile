## Author: Andrew Oliver
## Version: aoliver44/scfa_rsutdio:1.1,latest
## Date: Mar 18, 2024

## base image to start with
FROM rocker/rstudio:4.2

## RENV version
ENV RENV_VERSION=1.0.5

## install some things that R needs (for intel machines)
## uncomment below if running intel machine
#RUN apt-get update && apt-get install -y libz-dev build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libxml2-dev libglpk-dev libnode-dev libv8-dev

RUN install2.r --error remotes
RUN install2.r --error igraph networkD3

## install RENV, which will then install all R project packages
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_version('renv', version = '1.0.5', repos = 'http://cran.us.r-project.org')"

## should be in the same directory as this file
COPY renv.lock ./
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'
