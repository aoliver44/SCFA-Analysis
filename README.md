# Diet, microbiome, and inflammation predictors of fecal and plasma short-chain fatty acids in humans.
 Scripts and docker images used in paper: Diet, microbiome, and inflammation predictors of fecal and plasma short-chain fatty acids in humans. A paper on this work has been recently submitted.

______________
### **Data Availibility**
___________________
- Metagenomes are deposited in NCBI Sequence Read Archive (SRA) under the [study accession SRP354271](https://dataview.ncbi.nlm.nih.gov/object/PRJNA795985) and BioProject PRJNA1090850
 (release pending publication). Requests for non-metagenomic data from the USDA ARS WHNRC Nutritional Phenotyping Study used in this analysis should be made via an email to the senior WHNRC author on the publication of interest. Requests will be reviewed quarterly by a committee consisting of the study investigators.

________________________________________
### **Containers for reproducibility**
_________________________________________________
- A docker container for analysis in R. Docker must be installed to run. These images were built and run on Docker v4.26.1. To build these images, make sure you have cloned the git repo, and are inside that directory: 
```
git clone aoliver44/SCFA-Analysis
cd SCFA-Analysis
docker build -t scfa_analysis:1.0 .
```
- EVEN EASIER: you can just pull the image from docker hub (though i still recommend you clone the enviornment to have the scripts):
  ```
  docker pull aoliver44/scfa_analysis:latest
  ```


- **In order to run the container:**

    1. (assuming you are still in the cloned repository folder) ```docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `pwd`:/home/docker/ scfa_analysis:1.0``` (or whatever you named the container)
    2. navigate to http://localhost:8787/ in a browser window
    3. log into the Rstudio local server
        - username: rstudio
        - password: yourpasswordhere (if you didnt set one in docker run command)
    4. change to the scripts working directory inside R.
        - setwd("/home/docker")
    5. navigate the filesystem to the working directory.
    
        ![plot showing changing working directory in the file pane](https://github.com/aoliver44/SCFA-Analysis/blob/main/utilities/readme_picture.png)
_____________
### **Figures**
_____________
- to generate most of the main/supp figures and tables in the manuscript, source each the script that is named that for the figure
____________
### **Computing Environments**
_______________
- all code used in this analysis is provided
- R analyses were run locally on a intel-based macbook pro
- ML analyses were run on Ceres, a supercomputer for use by USDA researchers
- Sequence preprocessing was run on Spitfire, a slurm-based HPC cluster managed by the UC Davis Genome Center
    - code for sequencing preprocessing can be found [here](https://github.com/dglemay/ARG_metagenome)

