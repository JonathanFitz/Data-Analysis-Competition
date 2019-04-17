This project is the submission by Jonathan Fitz and Dong Liang for the
2017 Statistical Society of Canada Data Analysis Competition. 

Supervised by Shaun Sun of the University of the Fraser Valley, this project
marks the first time our university has been represented at this annual
competition.


PURPOSE

The Statistical Society of Canada uses a national survey to collect various
health-related measurements. We were provided an example of such a dataset
consisting of 3060 observations and 509 variables, of which 500 represented
bootstrap weights. The purpose of these bootstrap weights was to approximately
take into account the effect of the multi-level stratification used in the
survey design.

Generally, the aim of the analysis was two-fold:

    1) Understand the risk factors associated with hypertension (high blood
       pressure) among Canadian adults. Do these risk factors change across
       different segments of the population?

    2) Explore how the bootstrap weights affect the above analysis.

Although the approach and conclusions are described in our poster (see
Final_Poster.pdf), I briefly summarize both of these below.


APPROACH

The programming language R was used to perform all parts of the analysis,
including data exploration, cleansing, visualization, and statistical
analysis. Notably, the R package ggplot2 was used for all the visualizations.
A Generalized Linear Model was used to perform the hypothesis testing.


CONCLUSIONS

For both models (with and without bootstrap weights) there is very strong
evidence that age, sex, and bmi affect mean blood pressure, even when other
risk factors are accounted for. For the unweighted model, there is also some
evidence that mercury levels affect mean blood pressure.

There is also some evidence of interaction effects. For example, in the
unweighted model there is some evidence that the effect of bmi on mean blood
pressure changes with sex.

For the weighted model, there is some evidence that the effect of sex and bmi
on mean blood pressure changes with age.

I again refer the reader to our poster (Final_Poster.pdf) for a more in-depth
summary of the project.
