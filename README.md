## Research compendium for 'Design space constraints and the cultural taxonomy of European Final Palaeolithic large tanged points - a comparison of typological, landmark-based and whole-outline geometric morphometric approaches' 

### Compendium DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4560743.svg)](https://doi.org/10.5281/zenodo.4560743)

The files at the URL above will generate the results as found in the publication. The files hosted at <https://github.com/yesdavid/designspace_culttax_article_2021> are the development versions and may have changed since the paper was published.

### Maintainer of this repository:

[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--7349--5401-green.svg)](http://orcid.org/0000-0001-7349-5401) David N. Matzig (<david.matzig@cas.au.dk>) 

### Published in:

Matzig, D. N., Hussain, S. T., & Riede, F. (2021). Design Space Constraints and the Cultural Taxonomy of European Final Palaeolithic Large Tanged Points: A Comparison of Typological, Landmark-Based and Whole-Outline Geometric Morphometric Approaches. Journal of Paleolithic Archaeology, 4(4), 27. https://doi.org/10.1007/s41982-021-00097-2 [![DOI](https://zenodo.org/badge/DOI/10.1007/s41982-021-00097-2.svg)](https://doi.org/10.1007/s41982-021-00097-2) 

### Abstract:

The identification of material culture variability remains an important goal in archaeology, as such variability is commonly coupled with interpretations of cultural transmission and adaptation. While most archaeological cultures are defined on the basis of typology and research tradition, cultural evolutionary reasoning combined with computer-aided methods such as geometric morphometrics (GMM) can shed new light on the validity of many such entrenched groupings, especially in regard to European Upper Palaeolithic projectile points and their classification. Little methodological consistency, however, makes it difficult to compare the conclusion of such studies. Here, we present an effort towards a benchmarked, case-transferrable toolkit for such approaches, comparatively exploring relevant techniques centered on outline-based GMM. First, we re-analyze and compare two previously conducted landmark-based analyses of stone artefacts with our whole-outline approach, demonstrating that outlines can offer an efficient and reliable alternative. We then show how a careful application of clustering algorithms to GMM outline data is able to successfully discriminate between distinctive tool shapes, and suggest that such data can also be used to infer cultural evolutionary histories matching already observed typo-chronological patterns. Building on this baseline work, we apply the same methods to a dataset of large tanged points from the European Final Palaeolithic (ca. 15,000-11,000 cal BP). Exploratively comparing the structure of design space within and between the in this way analyzed datasets, our results indicate that Final Palaeolithic tanged point shapes do not fall into meaningful regional or cultural evolutionary groupings but exhibit an internal outline variance comparable to spatio-temporally much closer confined artefact groups of post-Palaeolithic age. We discuss these contrasting results in relation to the architecture of lithic tool design spaces and technological differences in blank production and tool manufacture.

### Keywords:

Final Palaeolithic; lithic technology; geometric morphometrics; outline analysis; design space; cultural evolution

### Overview of contents and how to reproduce:

This repository contains data (`1_data`) and code (`2_code`) for the paper. After downloading, the results can be reproduced using `matzig_et_al_2021.Rproj` and the existing folder structure. The required packages and their versions which have been used in this study are listed below and in the `DESCRIPTION`-file. All analyses and visualisations presented in this paper were prepared in R 3.6.3 under Ubuntu 18.04.5 LTS (64-bit).

### Required R-packages and their versions:

`ape` (>= 5.4-1), `asbio` (>= 1.6-5), `BiocManager` (>= 1.30.10), `caret` (>= 6.0-86), `cluster` (>= 2.1.0), `coin` (>= 1.3-1), `cowplot` (>= 1.0.0), `dispRity` (>= 1.5.1), `dplyr` (>= 1.0.0.9000), `fpc` (>= 2.2-7), `ggimage` (>= 0.2.8), `ggplot2` (>= 3.3.2), `ggpubr` (>= 0.4.0), `ggthemes` (>= 4.2.0), `ggtree` (>= 2.0.4), `maptree` (>= 1.4-7), `Momocs` (>= 1.3.0), `MVN` (>= 5.8), `NbClust` (>= 3.0), `parallel` (>= 3.6.3), `phangorn` (>= 2.5.5), `psych` (>= 1.9.12.31), `randomcoloR` (>= 1.1.0.1), `raster` (>= 3.3-7), `readr` (>= 1.3.1), `rgeos` (>= 0.5-3), `rworldmap` (>= 1.3-6), `splitstackshape` (>= 1.4.8), `vegan` (>= 2.5-6).




### Licenses:

Code: MIT <http://opensource.org/licenses/MIT> year: 2021, copyright holder: David Nicolas Matzig

Data: The data has been compiled from different resources. Please see the article's references or this repository under `./1_data/Readme.md`.
