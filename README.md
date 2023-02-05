# Code examples for doctoral thesis of Anna-Lena Wölwer

This project contains code for chosen parts of the doctoral thesis titled *Model-Based Prediction and Estimation Using Incomplete Survey Data*, submitted in partial fulfillment of the requirements for the degree *Doctor rerum politicarum* (Dr. rer. pol.) to the Faculty IV at Trier University by Anna-Lena Wölwer. The thesis is freely available for download under the following [DOI: 10.25353/ubtr-xxxx-25a6-5f2c](https://doi.org/10.25353/ubtr-xxxx-25a6-5f2c).

The doctoral thesis was submitted in June 2022, defended in December 2022 and published in January 2023.

## Structure

### Multivariate Fay-Herriot Models under Missing Direct Estimates (Chapter 6)

-   *MMFH_example.qmd*: Shows a commented example of generating data according to a MMFH model (using *MMFH_gen_dat_m3.R*) and estimating its parameters using algorithm *MMFH_fitting.R*. The presented example can be used as the basis for a model-based Monte Carlo simulation study similar to the ones shown in Chapter 6.
-   *MMFH_gen_dat_m3*: Code for generating the different data parts used in model-based Monte Carlo simulation studies of multivariate Fay-Herriot Models without and with randomly missing direct estimates. The code is written for 3 dependent variables.
-   *MMFH_fitting.R*: Code for fitting a multivariate Fay-Herriot model under missing direct estimates (MMFH). For fitting, we use a Fisher-Scoring algorithm, either with ML or REML.

### Empirical Best Prediction in Multivariate Fay-Herriot Models (Chapter 5)

-   XXXX: Under construction
