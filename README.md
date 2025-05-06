# TRIP

This repository contains the code for replicating the analysis and visuals presented in .... This paper is a distillation of my master's thesis, which I completed under the mentorship of Danny Krizanc, with lots of really valuable input from Professor Valerie Nazzaro. The thesis should be under an embargo until May 2026, but hopefully there will be a paper on the internet you can check out! If you don't see anything linked here it probably means that we haven't published anything yet :(. Below is a table with file names and a description of the code that is in each.


| File | Description |
| --- | --- |
| ch3.ipynb | The content of this chapter didn't make it to the paper. I consider the increase in efficiency for the algorithm that results from eliminating repeat computations. Under some assumptions we can get very clean recurrences, which perform decently in practice as estimates. |
| curseOfDimensionality.R |The test breaks down in higher dimensions, which we demonstrate through simulation in this file.|
| pairNoPair.R |It should be considered whether the permutations treat the same row with different features as the unit to be exchanged or if the exchanging should ignore this. We argue that since the similarity of a point to it's corresponding row-pair is almost always greater than the similarity to other observations, we should exchange within rows (resembling the Wilcoxon signed-rank test). We explore this question through simulation.|
| singh.R | The application to the gene expression dataset used in Singh et al. I'm not going to provide the data in this repository, but you should be able to download it for yourself [here](https://www.stats.uwo.ca/faculty/aim/2015/9850/microarrays/FitMArray/data/).|
| spcaSim.R| We demonstrate the effectiveness of [SPCA](https://www.jstor.org/stable/27594179?seq=1) for producing interpretable components that allow the extension of TRIP to higher dimension.|
| wineData.R| The application to the wine classification dataset. Again, you can find the data and easily download it from the [UCI Machine Learning Repository](https://archive.ics.uci.edu/dataset/109/wine).|
