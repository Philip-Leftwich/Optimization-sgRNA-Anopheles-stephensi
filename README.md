# Optimization of sgRNA expression with RNA Pol III regulatory elements in Anopheles stephensi

# Authors 

Estela Gonzalez1,2†, Michelle A. E. Anderson1,3†, Joshua X. D. Ang1,3, Katherine Nevard1, Lewis Shackleford1,3, Mireia Larrosa-Godall1, Philip T. Leftwich4, Luke Alphey1,3*

†these authors contributed equally

*Email address: luke.alphey@york.ac.uk

1 Arthropod Genetics, The Pirbright Institute, Pirbright GU24 0NF, U.K. 

2 Animal and Plant Health Agency, Woodham Lane, Addlestone, Surrey, KT15 3NB

3 Department of Biology, University of York, Wentworth Way, York, YO10 5DD

4 School of Biological Sciences, University of East Anglia, Norwich Research Park, Norwich, NR4 7TJ, U.K.


# Abstract

*Anopheles stephensi*, a major Asian malaria vector, is invading Africa and has been implicated in recent outbreaks of urban malaria. Control of this species is key to eliminating malaria in Africa. Genetic control strategies, and CRISPR/Cas9-based gene drives are emerging as promising species-specific, environmentally friendly, scalable, affordable methods for pest control. To implement these strategies, a key parameter to optimize for high efficiency is the spatiotemporal control of Cas9 and the gRNA. Here, we assessed the ability of four RNA Pol III promoters to bias the inheritance of a gene drive element inserted into the cd gene of *An. stephensi*. We determined the homing efficiency and examined eye phenotype as a proxy for non-homologous end joining (NHEJ) events in somatic tissue. We found all four promoters to be active, with mean inheritance rates up to 99.8%. We found a strong effect of the Cas9-bearing grandparent (grandparent genotype), likely due to maternally deposited Cas9.

# Repository

Contains raw data and analyses for inheritance and cleavage rates for Figures 2 and 3 from published manuscript. 

- homing_mosaic_analysis.R:  produces models and figures for Figure 2
- cutting_analysis.R: produces models and figures for Figure 3
- model_summaries.R: produces model summary tables found in Supplementals

## R packages used

Versions can be initialised from renv.lock file included with the repository

|Package  |Version |
|:--------|:-------|
|base     |4.3.3   |
|DHARMa   |0.4.7   |
|emmeans  |1.10.7  |
|renv     |1.1.1   |
|scales   |1.3.0   |
|showtext |0.9.7   |
|sjPlot   |2.8.17  |
|glmmTMB   |1.1.10  |
|tidyverse |2.0.0  |


We used R version 4.3.3 (R Core Team 2024) and the following R packages: DHARMa v. 0.4.7 (Hartig 2024), emmeans v. 1.10.7 (Lenth 2025), renv v. 1.1.1 (Ushey and Wickham 2025), scales v. 1.3.0 (Wickham, Pedersen, and Seidel 2023), showtext v. 0.9.7 (Qiu and See file AUTHORS for details. 2024), sjPlot v. 2.8.17 (Lüdecke 2024), glmmTMB v 1.1.10 (Brooks et al. 2017), tidyverse v 2.0.0 (Wickham et al. 2019).

Package citations:

Brooks, Mollie E., Kasper Kristensen, Koen J. van Benthem, Arni Magnusson, Casper W. Berg, Anders Nielsen, Hans J. Skaug, Martin Maechler, and Benjamin M. Bolker. 2017. “glmmTMB Balances Speed and Flexibility Among Packages for Zero-Inflated Generalized Linear Mixed Modeling.” The R Journal 9 (2): 378–400. https://doi.org/10.32614/RJ-2017-066.

Hartig, Florian. 2024. DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. http://florianhartig.github.io/DHARMa/.

Lenth, Russell V. 2025. emmeans: Estimated Marginal Means, Aka Least-Squares Means. https://rvlenth.github.io/emmeans/.

Lüdecke, Daniel. 2024. sjPlot: Data Visualization for Statistics in Social Science. https://CRAN.R-project.org/package=sjPlot.

Qiu, Yixuan, and authors/contributors of the included software. See file AUTHORS for details. 2024. showtext: Using Fonts More Easily in r Graphs. https://github.com/yixuan/showtext.

R Core Team. 2024. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing. https://www.R-project.org/.

Ushey, Kevin, and Hadley Wickham. 2025. renv: Project Environments. https://rstudio.github.io/renv/.

Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy D’Agostino McGowan, Romain François, Garrett 
Grolemund, et al. 2019. “Welcome to the tidyverse.” Journal of Open Source Software 4 (43): 1686. https://doi.org/10.21105/joss.01686.

Wickham, Hadley, Thomas Lin Pedersen, and Dana Seidel. 2023. scales: Scale Functions for Visualization. https://scales.r-lib.org.