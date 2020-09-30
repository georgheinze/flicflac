# flicflac
SAS-macros for Firth's corrected logistic, conditional logistic and Poisson regression, FLIC and FLAC methods

## Description

Here we provide our SAS-macros to fit Firth-corrected regression models, in particular logistic, conditional logistic and Poisson regression models. Special macros are available to implement the FLIC and FLAC methods of Puhr et al (2017) <doi:10.1002/sim.7273>.

## LogisticRegression/FL.SAS

The 'old' SAS macro to fit logistic regression models using SAS/PROC IML code. Unlike the implementation of Firth's correction in SAS/PROC LOGISTIC, this macro is also able to provide p-values based on penalized likelihood ratio tests for each regression coefficient. It also computes profile penalized likelihood confidence intervals as described by Heinze and Schemper (2002), Heinze and Ploner (2003),  Heinze (2006), and Mansournia, Geroldinger, Greenland and Heinze (2018).

## ConditionalLogisticRegression/CFL.SAS

Implements the conditional Firth-corrected logistic regression methods described in Heinze and Puhr (2010). It uses SAS/PROC LOGISTIC to compute the conditional distribution of sufficient statistics which can computationally burdensome.

## CoxRegression/FC06.ZIP

Implements Firth's correction for Cox regression as described by Heinze and Schemper (2001). It is based on FORTRAN code and an external routine (either included as EXE or DLL) which must be made invokable from the SAS macro. The accompanying technical report contains instructions on how to implement analyses for 1:m matched case-control studies with Firth-correction using the macro.

## LogisticRegression/FLICFLAC.SAS

This macro implements the FLIC and FLAC methods as described by Puhr, Heinze, Nold, Lusa and Geroldinger (2017). These methods are particularly interesting for predicting with penalized logistic regression. Unlike the default Firth correction, with FLIC and FLAC it is guaranteed that the average predicted probability is equal to the observed event rate. With rare events, Firth correction can lead to inflated average predicted proabilities such that predictions are biased high. Recently, van Calster, van Smeden, de Cock and Steyerberg (2020) showed that the FLIC method can yield calibration slopes which have mean squared error smaller than competing methods that use cross-validation to tune penalty parameters such as the Lasso or ridge regression .

## PoissonRegression/FLACPOISSON.SAS

With this macro, the Firth and FLIC/FLAC methods can be used with Poisson and Negative Binomial regression. The macro builds on iterated calls of PROC GENMOD. Multiple, equally-structured data sets can be processed with very efficient use of BY-processing.

## Acknowledgment

This work was supported by the Austrian Science Fund (FWF), award I-2276 and by the European Commission's MSCA programme, grant agreement number 795292.

## References

Heinze, G., Schemper, M. (2001): "A Solution to the Problem of Monotone Likelihood in Cox Regression", Biometrics 57(1):114-119 <doi:10.1111/j.0006-341X.2001.00114.x>.

Heinze, G., Schemper, M. (2002): "A Solution to the Problem of Separation in logistic regression", Statistics in Medicine 21:2409-2419 <doi:10.1002/sim.1047>.

Heinze, G., Ploner, M. (2002): "SAS and SPLUS programs to perform Cox regression without convergence problems", Computer Methods and Programs in Biomedicine 67:217-223 <10.1016/S0169-2607(01)00149-3>.

Heinze, G., Ploner, M. (2003): "Fixing the nonconvergence bug in logistic regression with SPLUS and SAS", Computer Methods and Programs in Biomedicine 71:181-187  <doi:10.1016/S0169-2607(02)00088-3>.

Heinze, G. (2006): "A comparative investigation of methods for logistic regression with separated or nearly separated data", Statistics in Medicine 25:4216-4226 <doi:10.1002/sim.2687>.

Heinze, G., Dunkler, D. (2008): "Avoiding infinite estimates of time-dependent effects in small-sample survival studies", Statistics in Medicine 27:6455-6469 <doi:10.1002/sim.3418>.

Heinze, G., Puhr, R. (2010): "Bias-reduced and separation-proof conditional logistic regression with small or sparse data sets", Statistics in Medicine 29:770-777 <doi:10.1002/sim.3794>.

Puhr R, Heinze G, Nold M, Lusa L, Geroldinger A (2017): "Firth's logistic regression with rare events: accurate effect estimates and predictions?", Statistics in Medicine 36:2302-2317 <doi:10.1002/sim.7273>.

Mansournia MA, Geroldinger A, Greenland S, Heinze G (2018): "Separation in Logistic Regression: Causes, Consequences, and Control", Am J Epidemiol 187:864-870 <doi:10.1093/aje/kwx299>.

Van Calster, B., van Smeden, M., De Cock, B., Steyerberg, S. (2020). "Regression shrinkage methods for clinical prediction models do not guarantee improved performance: Simulation study", Statistical Methods in Medical Research, to appear <doi:10.1177/0962280220921415>.




