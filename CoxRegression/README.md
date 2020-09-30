# Firth's correction for Cox regression

As of SAS Version 9, Firth's correction has been implemented for Cox regression in PROC PHREG also, including profile penalized likelihood confidence intervals. However, in PROC PHREG p-values are always based on Wald statistics rather than penalized likelihood ratio tests, and hence are not fully compatible with the confidence intervals.

The SAS macro FC06.SAS allows to compute p-values based on penalized likelihood ratio tests which are compatible with the profile penalized likelihood confidence intervals. Furthermore, the macro can handle time-dependent effects and counting-process style to represent survival times with left-truncation.

The core routines reside in a DLL or EXE file which must be made callable from the macro. This requires some installation steps which are outlined in the Technical Report. You can download the full suite of files needed to run the macro as a zip file.
