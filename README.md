# Novel Classification of Human Ulcers with Host Gene Expression
---------------------------

This markdown will describe the data, scripts, programs, etc. required for reproducing the 
analysis and figures or the article: *Novel Classification of Human Ulcers with Host Gene Expression*.


# Scripts

---------------------------

Scripts for running tha analysis are located in the `./scripts` directory.

`human_alignment.py` - Python script that performs QC, trimming, rRNA depletion, alignment, and feature 
counting on the raw, demultiplexed read files.   
`analysis_utils.R` - Contains functions for running the analysis in R.
Needs to be sourced before running any of the other R analysis scripts.                      
`FINAL_DataAnalysis.R` - This script takes the featureCounts output files and metadata. 
It then performs all of the analysis, including normalization, PCA, k-means clustering, etc.  
`FINAL_Makeplots.R` - This script generates all of the figures included in the manuscript. 
This excludes Figure 1a, which was created manually in Biorender.      
`ML_linearclassifier_mainCLEAN.py` - This script takes the metadata and count file output from 
the *FINAL_DataAnalysis.R* script and uses them as input to train a linear classifier. This 
classifier is then used to identify gene features which are good for classifying each level
of the metadata variable of interest. 

# Generating read counts from raw sequence files (fastq.gz)

---------------------------

It is possible to recreate the analysis from the raw sequence files, however may be difficult to adapt to
different high-performance computing environments. 

### Programs Required

-----------------
The following programs should be available in the envronment. 

Prereqs for Processing of Sequence Data:
```
cutadapt 2.4
anaconda3 v4.4.0 (Python version 3.6)
ncbi-blast v2.10.0
kraken v2.0.8
sortmerna v2.1
seqtk v1.3
bwa v0.7.16a
samtools v1.10
subread v1.6.2
pigz v2.3.4
```

R libraries required for Analysis and Figure Generation*: 
```
dplyr v1.0.5  
stringr v1.4.0  
tibble v3.1.0  
factoextra v1.0.7  
ggplot2 v3.3.3  
gridExtra v2.3
DESeq2 v1.28.1
tidyverse v1.3.0
```
*Information about the R session and programs used to create the data 
for the manuscript is given in: `./data/Example_data/sessionInfo.txt`

### Inputting raw sequence files and running the analysis

------------------

The input is based on a tsv file, where the first row is a header containing 
the following metadata variables (see 
`./data/Example_data/Example_metadata_allsamples.tsv` for an example): 
  
`Sample_name`- The name of the sample       
`R1_name`- The name of the R1 file (forward read)  
`R2_name`- The name of the R2 file (reverse read). If single end, then R2_name
should be NA  
`data_dir`- The full path to the directory containing the R1 and R2 files   
`Endedness`- either SE for single-ended samples, or PE for paired-end samples  

This tsv file is then read, line-by-line, by the job submission script (which is
used by the TORQUE Resource
Manager to run the jobs) to define the jobs making up the entire batch job
submission. An example of the script is 
located at `./data/Example_data/Submission_script.pbs`. 
Several variables, which are required to run the analysis are also defined in the
job submission script and will likely need to be redefined when running the analysis
on a new system.
The job submission scrips then passes all necessary variables to
`human_alignment.py`, which then runs the analysis 
and writes all generated files to the output directory defined in the job submission script. 

Alternatively, the `human_alignment.py` script could also be run by manually passing the 
required arguments to it. 

### Analysis of Count Data and Generation of Figures

---------------------------------------

All analyses in the manuscript can be reproduced with the included R scripts. 
To reproduce the analysis, be sure that all necessary libraries are installed
with the correct versions. Also, be sure that the following folders exist and
that the WoundClassification2021 repository is set as the R source directory:  

    ./analysis/Figures/  
    ./analysis/DESeq2/      
    ./analysis/GO_analysis/    
    ./analysis/kmeans/    
    ./analysis/linearSVM_out/      
    ./analysis/normalized_counts/   


The count data used for the analysis is provided in `./data/Example_data/Example_counts.csv`. 
The corresponding metadata file used for the analysis is 
`./data/Example_data/Example_metadata_Ranalysis.csv`.

**Note** All analysis should be run from the main repo folder `./WoundClassification2021` as the working directory. 

To run the analysis, do the following: 
1. Download this git repository and unzip if necessary.
2. Download and unzip the 
[annotation file](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/)
from NCBI (GCA_000001405.15_GRCh38_full_analysis_set.fna.gz). Place it into `./data/annotation_gff/`. 
2. Open an R or R studio session and set the working directory as the repository folder. 
3. Open and run all of FINAL_DataAnalysis.R. This should run without errors and write output to the folders defined above.
Some warnings may appear, especially regarding conversion of characters to factors.  
    
    **Warning**: The analysis will overwrite files without asking, if there are already files with the same name
    as the output files in the output directories.

4. Rerun the linear classifier script by running the command `python3 ./scripts/ML_linearclassifier_mainCLEAN.py`
   from the command line while in the main directory (`WoundClassification2021`). This will generate the coefficient figures
   required for Figure 3.  
   
    **Warning**: The classifier script has only been tested with `pandas v1.0.3`. Other versions do not work properly.
    I have tried with pandas version 1.2.4, but it does seem not handle NA/NaNs correctly and includes them as levels
    in the model. 
   
5. Run the FINAL_Makeplots.R script to generate Figures 1-4. 
