/* Adapted from LogisticRegression/augGEE1.sas (georgheinze/flicflac) */
/* %augGEE1 is reproduced verbatim from the repo. Below it is a small   */
/* caller with mock clustered binary-outcome data (the repo ships no    */
/* example caller for its macros) exercising the single-step augmented */
/* GEE approach of Geroldinger, Blagus, Ogden and Heinze (2022). */

%macro augGEE1(data=_last_, y=, varlist=, id=);

* Version: 20220812;
* implements single-step augmented GEE (augGEE1) approach of Geroldinger et al (BMC Med Res Meth 2022);
* uses a compound symmetry working correlation matrix;
* empirical (sandwich) variance estimation without small sample correction term.;

data _tmp;
set &data;
_id_ = &id;
keep y &varlist _id_;
run;

proc logistic descending data=_tmp;
title3 "Firth (FL) model";
model &y = &varlist /firth;
output out=_flactmp0_ h=_hdiag_;
run;

data _flactmp_;
set _flactmp0_;
_added_covariate_=0;
_weight_=  1;
_ind_=1;
output;
_added_covariate_=1;
_weight_= _hdiag_ * 0.5;
_ind_=2;
output;
y=1-y;
_ind_=3;
output;
run;


data _flactmp2_;
set _flactmp_;
run;

proc sort data=_flactmp2_;
by _ind_ _id_;
run;

data _flactmp2_;
set _flactmp2_;
by _ind_ _id_;
retain _id2_ 0;
if first._id_ then _id2_=_id2_+1;
run;


*pseudo observations in extra clusters depending on y or 1-y;
proc gee data=_flactmp2_ descending;
ods output GEEEmpPEst=_AUGEST2_ GEEExchCorr=_correst2_;
title3 "augmented GEE model";
class _id2_;
model &y = &varlist/ dist=bin;
repeated subject=_id2_ / corrw  type=cs;
weight _weight_;
run;


%mend;

/* Mock clustered rare-event data: binary outcome, one id per row (each */
/* its own cluster), matching what %augGEE1's &id argument expects. */
data rare_event;
input id y x1 x2;
datalines;
1 1 5.2 1.1
2 0 1.4 0.3
3 0 1.1 0.2
4 0 0.9 0.4
5 0 1.6 0.5
6 0 1.2 0.3
7 0 0.8 0.2
8 0 1.5 0.6
9 0 1.0 0.3
10 0 0.7 0.1
11 0 1.3 0.4
12 0 0.9 0.2
13 0 1.1 0.3
14 0 0.6 0.1
15 0 1.4 0.5
16 0 1.0 0.2
17 0 0.8 0.3
18 0 1.2 0.4
19 0 0.9 0.1
20 0 1.1 0.2
;
run;

%augGEE1(data=rare_event, y=y, varlist=x1 x2, id=id);
