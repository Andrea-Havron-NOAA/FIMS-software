
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FIMS-software

<!-- badges: start -->
<!-- badges: end -->

The goal of FIMS software platform project is to compare statistical
computing environments for use in the
[FIMS](https://www.fisheries.noaa.gov/national/population-assessments/fisheries-integrated-modeling-system)
intiative. This is ongoing research. For more details, see the [Working
Document](https://andrea-havron-noaa.github.io/FIMS-software/docs/FIMS-software-exploratory-analysis.html)

## Resources

-   [Julia GitHub](https://github.com/JuliaLang/julia)
-   [Stan User’s Guide](https://mc-stan.org/docs/2_19/stan-users-guide/)
-   [Stan GitHub](https://github.com/stan-dev/stan)
-   [cmdstan GitHub](https://github.com/stan-dev/stan)
-   [TMB
    Documentation](https://kaskr.github.io/adcomp/_book/Introduction.html)
-   [TMB GitHub](https://github.com/kaskr/adcomp)

## FIMS Software Organization

The FIMS software platform project is organized as follows:

<table style="width:68%;">
<colgroup>
<col style="width: 22%" />
<col style="width: 45%" />
</colgroup>
<thead>
<tr class="header">
<th>Directory</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>data</td>
<td>data simulation script, simulated data files, and initial values files</td>
</tr>
<tr class="even">
<td>docs</td>
<td>working documents detailing project and results</td>
</tr>
<tr class="odd">
<td>R</td>
<td><ul>
<li>model_setup.R: code used to setup TMB and Stan models<br />
</li>
<li>modular_demo.R: code runs modular example of passing C++ code to different platforms<br />
</li>
<li>run_models.R: runs different platforms for benchmark study excluding julia<br />
</li>
<li>utils.R: base functions for R code</li>
</ul></td>
</tr>
<tr class="even">
<td>results</td>
<td>results files from benchmark study</td>
</tr>
<tr class="odd">
<td>src</td>
<td>source files for project</td>
</tr>
<tr class="even">
<td>  src/julia</td>
<td><ul>
<li>stateSpace.jl: julia file that specifies and run gompertz and logistic models</li>
</ul></td>
</tr>
<tr class="odd">
<td>  src/Rcpp</td>
<td><ul>
<li>Common.hpp:<br />
</li>
<li>likelihoods.hpp: likelihood definitions<br />
</li>
<li>logisticGrowth.cpp: TMB model specification<br />
</li>
<li>model.hpp: C++ model definitions</li>
</ul></td>
</tr>
<tr class="even">
<td>  src/stan</td>
<td><ul>
<li>gompertz.stan: gompertz model<br />
</li>
<li>logistic.stan: logistic growth model</li>
</ul></td>
</tr>
<tr class="odd">
<td>  src/tmb</td>
<td><ul>
<li>sptail_poisson.cpp: spatially explicit Poisson model<br />
</li>
<li>stateSpace.cpp: gompertz and logistic growth model</li>
</ul></td>
</tr>
</tbody>
</table>

[Example](https://andrea-havron-noaa.github.io/FIMS-software/docs/ar1xar1.html)
incorporating TMB functions in a modular C++ framework
