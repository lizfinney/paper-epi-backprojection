
# Overview

This repository contains and data and scripts for reproducing the results accompanying the manuscript  

### Back-projection improves inference from sparsely sampled genomic surveillance data
Elizabeth Finney<sup>1</sup>, Brian Lee<sup>1</sup>, Syed Faraz Ahmed<sup>2,3</sup>, Muhammad Saqib Sohail<sup>4</sup>, Ahmed Abdul Quadeer<sup>2</sup>, Matthew R. McKay<sup>2,3,5</sup>, and John P. Barton<sup>1,6,#</sup>  

<sup>1</sup> Department of Physics and Astronomy, University of California, Riverside  
<sup>2</sup> Department of Electrical and Electronic Engineering, University of Melbourne  
<sup>3</sup> Department of Microbiology and Immunology, University of Melbourne, at The Peter Doherty Institute for Infection and Immunity  
<sup>4</sup> Department of Computer Science, Bahria University  
<sup>5</sup> Victorian Infectious Diseases Reference Laboratory, Royal Melbourne Hospital, at the Peter Doherty Institute for Infection and Immunity  
<sup>6</sup> Department of Computational and Systems Biology, University of Pittsburgh School of Medicine  
<sup>#</sup> correspondence to [jpbarton@pitt.edu](mailto:jpbarton@pitt.edu)  

The preprint is available at __INSERT LINK HERE__.

# Contents

1. Branching process simulations: The directory `src/simulations/` contains script to generate and analyze simulations. Generating simulations is done in the notebook `simulation-running.ipynb`, and `simulation-visualization.ipynb` contains analysis and paper figures.

2. Back-projection implementation: The directory `src/data-processing-pipeline/deconvolution` contains our implementation of back-projection for SARS-CoV-2, originally described by Becker et al. 1991. The notebook `back-projection-example.ipynb` can be run with example inputs found in this directory. This directory contains production scripts which we run with an incubation period distribution and binomial weights found in `/deconvolution/distribution-files`. This folder also contains a file describing the back-projected output of a single collection time, which aids in computation time.

3. Data processing: The notebook `SC2-back-projection-pipeline.ipynb` is used for processing and analyzing the SARS-CoV-2 sequence data. Back-projection is applied to pre-processed sequences, and back-projected sequences are further processed (i.e. time series trimming) before inferring mutational effects. Scripts used for both the back-projected method and comparison method (see [Lee et al. 2025](https://www.nature.com/articles/s41467-024-55593-0)) are found in `src/data-processing-pipeline/scripts`. The notebook also contains code that produces necessary job files to run much of the data analysis with cluster computing. The original sequence data and metadata used here were obtained from [GISAID](https://gisaid.org/)

4. Figures: A notebook for visualizing allele trajectories and their inferred coefficients is done by `SC2-back-projection-visualization.ipynb`. This notebook allows for comparison between back-projected results and comparison results obtained with a naive moving-average smoothing. An additional notebook for generating figures found in the paper and supplementary material is done by `SC2-back-projection-paper-figures.ipynb`. Processed trajectories and selection coefficients can be found in `data/analysis-visualization`.

### Software dependencies

Parts of the analysis are implemented in C++11 and the [GNU Scientific Library](https://www.gnu.org/software/gsl/).

# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).
