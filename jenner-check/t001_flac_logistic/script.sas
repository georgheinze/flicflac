/* Adapted from LogisticRegression/flicflac.sas (georgheinze/flicflac) */
/* %flac is reproduced verbatim from the repo. Below it is a small      */
/* caller with a mock rare-event dataset (this repo ships no example    */
/* caller for its macros) exercising %flac exactly as documented in the */
/* repo README: "the FLAC method of Puhr, Heinze, Nold, Lusa and        */
/* Geroldinger (2017)". */

%macro flac(data=_last_, y=, varlist=, weight=, tau=0.5);


proc logistic descending data=&data;
title3 "Firth (FL) model";
model &y = &varlist /firth;
output out=_flactmp1_ h=_hdiag_;
%if &weight ne %then %do; weight &weight; %end;
run;

data _flactmp2_;
set _flactmp1_;
_added_covariate_=0;
_weight_= %if &weight = %then %do; 1 %end; %else %do; &weight %end; ;
output;
_added_covariate_=1;
_weight_= _hdiag_ * &tau * %if &weight = %then %do; 1 %end; %else %do; &weight %end; ;
output;
y=1-y;
output;
run;

*** to fix: c-index not correct, as it involves the pseudo observations. Perhaps also do a 'dummy' call to proc logistic as in flic?;

proc logistic descending data=_flactmp2_ outest=_FLACEST_ outmodel=_flacmodel_;
title3 "FLAC model";
model &y = &varlist _added_covariate_;
weight _weight_;
run;


%mend;

/* Mock rare-event data: binary outcome with a low event rate, the      */
/* exact scenario FLIC/FLAC exists for (Firth correction on ordinary    */
/* /firth alone can bias predicted probabilities high with rare events).*/
data rare_event;
input y x1 x2;
datalines;
1 5.2 1.1
0 1.4 0.3
0 1.1 0.2
0 0.9 0.4
0 1.6 0.5
0 1.2 0.3
0 0.8 0.2
0 1.5 0.6
0 1.0 0.3
0 0.7 0.1
0 1.3 0.4
0 0.9 0.2
0 1.1 0.3
0 0.6 0.1
0 1.4 0.5
0 1.0 0.2
0 0.8 0.3
0 1.2 0.4
0 0.9 0.1
0 1.1 0.2
;
run;

%flac(data=rare_event, y=y, varlist=x1 x2);

proc print data=_FLACEST_ noobs label;
run;
