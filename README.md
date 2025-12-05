# landpred <img src="hex_landpred.png" align="right" height="220" alt="landpred hex logo" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/landpred)](https://CRAN.R-project.org/package=landpred )
<!-- badges: end -->

`landpred` is an `R` package that provides nonparametric models for landmark prediction of long-term survival outcomes, incorporating covariate and short-term event information. The package supports the construction of flexible varying-coefficient models that use discrete covariates, as well as multiple continuous covariates. The goal is to improve prediction accuracy when censored short-term events are available as predictors, using robust nonparametric procedures that don't require correct model specification and avoid restrictive parametric assumptions found in existing multistate survival methods. 

Earlier versions of this package (1.0, 1.1, 1.2) did not include the varying coefficient method. Version 2.0, developed by [Dylan Huynh](https://baolong281.github.io/) does include the varying coefficient method and other major improvements in package implementation. 

More information on these methods can be found in Parast et al. (2012, Journal of the American Statistical Association, [doi:10.1080/01621459.2012.721281](https://doi.org/10.1080/01621459.2012.721281), and 
Parast et al. (2011, Biometrical Journal, [doi:10.1002/bimj.201000150](https://doi.org/10.1002/bimj.201000150)). 

Go here to view a tutorial for this package: [landpred Tutorial](https://htmlpreview.github.io/?https://github.com/laylaparast/landpred/blob/main/landpredtutorial.html). 
