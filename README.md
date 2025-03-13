# cpp3d
Microbial metagenome functional profiling at the resolution of individual compounds (... and in 3-d!)

## Background
Microbiomes are everywhere and must be better harnessed to address significant challenges in climate and environmental change, food security, and sustainable management of soil, land and water. But we need to understand what they are doing, in different environments, in order to harness them.

Metagenomes are used to examine functions performed by microbiomes (i.e., microbial communities, their genes and metabolites). But, existing approaches require specialised knowledge of complex metabolic pathways and biological systems, and do not provide resolution at the deeper fundamental level of the chemical compounds involved (e.g., microbial substrates, feedstocks and metabolites such as carbohydrates, lipids, proteins, amino-sugars, vitamins and phytochemicals).

## A compound-oriented view can improve interpretability
Examining the role of microbiomes from a metabolite- or compound-oriented perspective offers promise for enhanced insights and interpretability. Chemical compounds are very often inputs, intermediates, or products of interest from microbiomes. In agriculture, ecosystem management, and microbiome-mediated health we typically have greater knowledge of the functional importance of compounds (e.g., chemical reactivity, health associations, quality indicators) than can be interpreted from overarching metabolic pathways reported by standard functional profiling tools. 

Therefore, stakeholders can benefit from improved insights and understanding of the roles of microbiomes through the innovative lens of their compound processing potential.

There are existing metabolite prediction frameworks, but these are typically limited to particular organisms and environments and are poorly suited to revealing potential metabolism in microbially diverse and fluctuating systems, such as soils.

## Compound processing potential (CPP) method overview

This updated CPP method uses outputs from existing [SUPER-FOCUS](https://github.com/metageni/SUPER-FOCUS) metagenome functional profiling software, which estimates the relative distribution of genetic capacity across various functional pathways. I link up existing look-up tables from the [ModelSEED biochemistry database](https://github.com/ModelSEED/ModelSEEDDatabase/tree/master/Biochemistry) to connect each functional pathway (or 'function') to one or more relevant chemical reactions, and the many compounds in each reaction. Functional relative abundances are then divided among all linked reactions and compounds with weighting to account for molar ratios (stoichiometry) of reactions.  Reaction inputs and products are all treated equally, because often microbes can facilitate a string of reactions (so products become inputs, and so on). For all carbon (C)-containing compounds, elemental ratios of oxygen (O):C, hydrogen (H):C, and nitrogen (N):C are calculated to enable visualisation and and analysis of energetically and chemically similar compounds. 

CPP values represent the functional capacity (%) allocated to each compound, reflecting their potential metabolism in a given metagenome.

## Code example

Bioinformatics for generating SUPER-FOCUS functional profiling data are provided via my earlier method, see [https://github.com/liddic/compound_potential] > [sunbad-resto](https://github.com/liddic/compound_potential/tree/main/sunbad-resto) and [forslund-t2d](https://github.com/liddic/compound_potential/tree/main/forslund-t2d) case studies which use published data from [Sun and Badgley 2019](https://doi.org/10.1016/j.soilbio.2019.05.004) and [Forslund et al 2015](https://doi.org/10.1038/nature15766).

Example R code for the case study analysis titled 'Overlapping soil and gut microbiome compound processing potential in a gradient of ecosystem quality and subjects with type 2 diabetes' is provided in the attached .R file. Outputs from this case study are also available in the attached .xlsx file.
