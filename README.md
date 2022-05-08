# MetTracer

## Introduction

`MetTracer` is an R package for tracing isotopic labeling experiments. For convenience, we build the docker image contains `MetTracer` package ([zhulab/mettracer-r](https://hub.docker.com/r/zhulab/mettracer-r)). Users can pull the image and run MetTracer without installing further packages.

## Pulling image
Users can pull the MetTracer-r image with the following script
``` bash
docker pull zhulab/mettracer-r
```

## Data preparation

 `MetTracer` requires the unlabelled metabolomics data to be pre-processed by xcms (see section 3.1.2 in [http://metdna.zhulab.cn/metdna/help](http://metdna.zhulab.cn/metdna/help)) and identified by MetDNA (see [http://metdna.zhulab.cn](http://metdna.zhulab.cn/)) first.

The following files are required by `MetTracer` for labeled metabolites extraction:

- A folder named &quot;unlabelled&quot; containing unlabelled data files (.mzXML), annotation result (.csv) and peak information (.Rda) for targe isotopologue generation
- A folder named &quot;labelled&quot; containing data files (.mzXML) for isotope labelling data extraction
- A R script named &quot;run.R&quot; for your data processing

 ### Preparing the labelled folder

The &quot;labelled&quot; folder contains all the data files requiring extraction of labeled metabolites. We permit multiple subfolders in it, and each subfolder is an independent biological group. The data files in the folder should be in mzXML format.

###  Preparing the unlabelled folder

The unlabeled folder contains unlabeled data files, annotation result and peak information.

- unlabeled data files
  - unlabeled data files are those used for peak detection and metabolite annotation. These data files should be in mzXML format and put into the subfolder.
- annotation information
  - annotation information is a .csv format table named &quot;MRN.annotation.result&quot;. The table should be put into the subfolder named  "MRN_annotation_result" under subfolder "MetProcess-Result&quot"
- peak information
  - peak information is a .Rda file which is generated after peak detection with XCMS.

### Preparing the R script

File run.R is an R source code file paralleled with "labelled" folder and "unlabelled" folder. We provide a template here. Users only need to change the folder name in general. Other parameters are recommended parameters.

- Parameters
  - `adj.label` All subfolder names in &quot;labelled&quot; folder, the data files in these folder will do isotope contamination correction.
  - `adj.unlabel` A subfolder in &quot;labelled&quot; folder. The data files in this subfolder are unlabeled, and they will be used as reference files to correct the isotope contamination for data files in all &quot;labelled&quot; subfolder.
 
``` R
wd <- getwd()
setwd(wd)
library(MetTracer)

isotopologueParam <- IsotopologueParam(rt.extend=15, value="maxo")
experimentParam = ExperimentParam(wd = wd, nSlaves = 6)
extractParam <- ExtractParam(d.extract = "labelled",
                             adj.label = c("Glu_Gln_Ac","Glu_Gln","Gln_Ac","Unlabel"),
                             adj.unlabel = c(rep("Unlabel", 4)), 
                             adj.contaminate = TRUE)
iso.targets <- GenerateIsotopologue(isotopologueParam, experimentParam)
pdParam <- PeakdetectionParam(peakwidth = c(5,30))
iso.peaks <- ExtractIsotopicPeaks(extractParam, pdParam,
                                  iso.targets,check.peaks = TRUE)
ExtractIsotopologue(extractParam, iso.peaks, iso.targets, correct.iso = TRUE)
```
### Overview of the data preparation
<img src="https://raw.githubusercontent.com/ZhuMetLab/MetTracer/main/data_structure.png" alt="Overview of the data preparation" title="Overview of the data preparation">

## Run data processing work with mettracer-r image
* go to your data folder (e.g., data) 
``` base
 cd data
```
* run docker using following code
``` bash
# MUST keep the code exactly as it is!
docker run -it --rm -v "$PWD":/data -u $(id -u ${USER}):$(id -g ${USER}) zhulab/mettracer-r Rscript run.R
```
* wait till data processing work done

* Explaining `docker run` arguments
* `-v "$PWD":/home/${USER}`: mapping current dirctory as home directory in docker container
*  `-u $(id -u ${USER}):$(id -g ${USER})`: using current user to run the container
* `Rscript ~/run.R`: run run.R in container home directory with `Rscript`  command

## License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a> 
This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)
