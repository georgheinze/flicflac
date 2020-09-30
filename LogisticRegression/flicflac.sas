%macro flic(data=_last_, y=, varlist=, weight=);

proc logistic descending data=&data outest=_FLMODTMP_ outmodel=_flmodel_;
ods output  CLoddsPL = _cloddsPL_;
title3 "Firth (FL) model";
model &y = &varlist /firth clodds=pl;
%if &weight ne %then %do; weight &weight; %end; 
output out=_FLTMP_ xbeta=_xbeta_ p=_pred_;
run;

proc logistic descending data=_fltmp_ outest=_FLICmodTMP_;
title3 "FLIC intercept modifier model";
model &y = /offset=_xbeta_;
%if &weight ne %then %do; weight &weight; %end;
run;

data _flicmodTMP_;
set _flicmodTMP_;
call symput("delta", intercept);
run;

data _FLmodTMP_;
set _FLmodTMP_;
intercept = intercept + &delta;
run;

proc logistic descending data=_fltmp_ inest=_FLmodTMP_ outest=_FLICest_ outmodel=_FLICmodel_;
title3 "FLIC model";
model &y = &varlist / firth maxiter=0;
%if &weight ne %then %do; weight &weight; %end;
output out=_FLICTMP_ p=_predflic_;
store out=_FLICstore_;
run;

proc print data=_cloddspl_ noobs label;
run;

%mend;

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




