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
docker build -t scfa_analysis:rstudio .
```
- EVEN EASIER: you can just pull the image from docker hub (though i still recommend you clone the enviornment to have the scripts):
  ```
  docker pull aoliver44/scfa_analysis:rstudio
  ```
________________________________________
### **To generate the figures from paper**
_________________________________________________

IMPORTANT: Assuming you cloned the GitHub repositories to your downloads foder. If this is not where you downloaded the GitHub repositories, change this

These commands will generate figures in a folder called ```~/Downloads/SCFA-Analysis/figure_scripts/output_figures/```

**Prior to generating figures**
```
cd ~/Downloads # if you download elsewhere, change throughout
git clone SCFA-Analysis-DATA ## PRIVATE REPO - For access, see instructions above
git clone aoliver44/SCFA-Analysis
```
**Generate a figure**
```
docker run --rm -it \
-v ~/Downloads/SCFA-Analysis/figure_scripts/:/home/scripts \
-v ~/Downloads/SCFA-Analysis-DATA/data/:/home/data \
-w /home \
scfa_analysis:rstudio bash -c "Rscript Figure1.R"
```
Note: you can generate any figure above by changing Figure1.R to another figure name. Options include:
1. Figure1.R
2. Figure2.R
3. Figure3.R
4. Figure4.R
5. Figure5.R
6. Table1.R
7. supplemental_figure2.R
8. supplemental_figure3.R
9. supplemental_figure4.R
10. supplemental_figure5.R
11. supplemental_figure6.R


**To work inside the container (if you want to look closer at the code and run it interactively):**

1. (assuming you cloned the repository to your Downloads folder) 

    ```docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v ~/Downloads/SCFA-Analysis/figure_scripts/:/home/scripts -v ~/Downloads/SCFA-Analysis-DATA/data/:/home/data -w /home scfa_analysis:rstudio```

2. navigate to http://localhost:8787/ in a browser window
3. log into the Rstudio local server
    - username: rstudio
    - password: yourpasswordhere (if you didnt set one in docker run command)
4. change to the scripts working directory inside R.
    - setwd("/home")
5. navigate the filesystem to the working directory.
   
    ![plot showing changing working directory in the file pane](https://github.com/aoliver44/SCFA-Analysis/blob/main/utilities/readme_picture.png)
____________
### **Computing Environments**
_______________
- all code used in this analysis is provided
- R analyses were run locally on a intel-based macbook pro
- Sequence preprocessing & ML analyses was run on Spitfire, a slurm-based HPC cluster managed by the UC Davis Genome Center
    - code for sequencing preprocessing can be found [here](https://github.com/dglemay/ARG_metagenome)

