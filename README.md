# Joint Mean and Correlation Regression Modelling for Multivariate Data

This repository code contain template code associated with the manuscript "Joint Mean and Correlation Regression Modelling for Multivariate Data" by [Tho](https://rsfas.anu.edu.au/about/staff-directory/zhi-yang-tho), [Hui](https://francishui.netlify.app/), and [Zou](https://cbe.anu.edu.au/about/staff-directory/dr-tao-zou).

# Getting started

There are currently three directories in this repository:

-   `Codes`, which contains `R` functions implementing the proposed joint mean and correlation regression model and the iterative estimation method along with other functions for performing inferences on the model parameters.

-   `Simulations`, which contains template scripts to implement the simulation study in the manuscript based on Bernoulli/Poisson/Negative Binomial/Gaussian responses.

-   `Applications`, which contains scripts applying the proposed model to the Scotland Carabidae ground beetle dataset (originally from [Ribera et al. 2001](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/0012-9658%282001%29082%5B1112%3AEOLDAS%5D2.0.CO%3B2)) and the France Alpine plant dataset (originally from the [aravo](https://cran.r-project.org/web/packages/lori/vignettes/aravo_data_analysis.html) package). For the ground beetle dataset, we only provide the environmental covariates and species traits that are used in the manuscript. **Users are recommended to start here by examining `carabidae_ground_beetle.R` or `aravo.R` to understand how to fit the proposed joint mean and correlation regression model and conduct statistical inference on its parameters.**

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem, then that is also much appreciated;
3.  Required data files etc...

Alternatively, please contact the corresponding author at [ZhiYang.Tho\@anu.edu.au](mailto:ZhiYang.Tho@anu.edu.au).
