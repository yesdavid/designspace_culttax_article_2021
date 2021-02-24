## Research compendium for 'Design space constraints and the cultural taxonomy of European Final Palaeolithic large tanged points -- a comparison of typological, landmark-based and whole-outline geometric morphometric approaches' 

### Compendium DOI:

<!-- [![DOI](https://img.shields.io/badge/DOI-XXX-blue)](https://doi.org/XXX) -->

The files at the URL above will generate the results as found in the publication. The files hosted at <https://github.com/yesdavid/designspace_culttax_article_2021> are the development versions and may have changed since the paper was published.

### Maintainer of this repository:

[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--7349--5401-green.svg)](http://orcid.org/0000-0001-7349-5401) David N. Matzig (<david.matzig@cas.au.dk>) 

### Published in:

TBA

### Abstract:

The identification of material culture variability remains an important goal in archaeology, as such variability is commonly coupled to interpretations of cultural transmission and adaptation. While most such archaeological cultures are defined on the basis of typology and research tradition, cultural evolutionary reasoning combined with computer-aided methods such as geometric morphometrics (GMM) sheds new light on the validity of many such entrenched groupings, especially in regard to Upper Palaeolithic projectile points in Europe. Little methodological consistency, however, makes it difficult to compare their conclusions. Here, we present an effort towards a benchmarked, case-transferrable toolkit for such approaches, comparatively exploring relevant techniques centered around outline-based GMM. First, we re-analyze and compare two previously conducted landmark-based analyses of stone artefacts against our whole-outline approach, finding that outlines offer an efficient and reliable alternative. We then show how a theory-driven application of clustering algorithms to GMM outline data is able to successfully discriminate between different distinctive tool shapes, and suggests that such data can even infer phylogenies matching typo-chronological patterns of cultural evolution. Finally, we apply these methods to a dataset of large tanged points from the North and Eastern European Final Palaeolithic (ca. 15,000-11,000 cal BP). Exploratively comparing the overall design spaces across all our datasets and highlighting the lack of internal structure amongst the Final Palaeolithic points, our results indicate that these point shapes do not fall into meaningful regional or cultural evolutionary groupings but exhibit an internal outline variance comparable to that of spatio-temporally much closer confined types of post-Palaeolithic assemblages.

### Keywords:

Final Palaeolithic; geometric morphometrics; outline analysis; design space; cultural taxonomy; cultural evolution

### Overview of contents and how to reproduce:

This repository contains data (`1_data`) and code (`2_code`) for the paper. After downloading, the results can be reproduced using `matzig_et_al_2021.Rproj` and the existing folder structure. The required packages and their versions which have been used in this study are listed below and in the `DESCRIPTION`-file. All analyses and visualisations presented in this paper were prepared in R 3.6.3 under Ubuntu 18.04.5 LTS (64-bit).

### Required R-packages and their versions:

`ape` (>= 5.4-1), `asbio` (>= 1.6-5), `BiocManager` (>= 1.30.10), `caret` (>= 6.0-86), `cluster` (>= 2.1.0), `coin` (>= 1.3-1), `cowplot` (>= 1.0.0), `dispRity` (>= 1.5.1), `dplyr` (>= 1.0.0.9000), `fpc` (>= 2.2-7), `ggimage` (>= 0.2.8), `ggplot2` (>= 3.3.2), `ggpubr` (>= 0.4.0), `ggthemes` (>= 4.2.0), `ggtree` (>= 2.0.4), `maptree` (>= 1.4-7), `Momocs` (>= 1.3.0), `MVN` (>= 5.8), `NbClust` (>= 3.0), `parallel` (>= 3.6.3), `phangorn` (>= 2.5.5), `psych` (>= 1.9.12.31), `randomcoloR` (>= 1.1.0.1), `raster` (>= 3.3-7), `readr` (>= 1.3.1), `rgeos` (>= 0.5-3), `rworldmap` (>= 1.3-6), `splitstackshape` (>= 1.4.8), `vegan` (>= 2.5-6).




### Licenses:

Code: MIT <http://opensource.org/licenses/MIT> year: 2021, copyright holder: David Nicolas Matzig

Data: The data has been compiled from different resources. Please see the article's references this repository under `./1_data/Readme.md`.
