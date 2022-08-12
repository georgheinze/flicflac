# SAS macros for logistic regression

* fl.sas: SAS macro to fit logistic regression with Firth's penalized likelihood 

This macro is somewhat outdated by FIRTH option in PROC LOGISTIC. However, unlike PROC LOGISTIC, it computes penalized likelihood ratio tests for all parameters which are compatible with the profile penalized likelihood confidence intervals.
Documented in fl.pdf.

* flicflac.sas: SAS macros to fit FLIC or FLAC models in logistic regression, as described by Puhr et al, StatMed 2017.

* augGEE1.sas: SAS macro to fit a GEE model with a binary outcome using the single-step augmented GEE (augGEE1) approach described by Geroldinger et al, BMC Med Res Meth 2022.

The macro uses an exchangeable working correlation matrix and on exit estimates the empiricial (sandwich) variance estimate.
