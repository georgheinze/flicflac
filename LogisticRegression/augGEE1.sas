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
