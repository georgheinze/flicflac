%macro flacpoisson(data=, varlist=, y=, maxiter=15, dist=poisson, outmodelfirth=_FIRTHMODEL, outmodelFLAC=_FLACMODEL, betastart=, offset=, Xconv=0.0001, by=, print=1, nonotes=0, iterinfo=1, pl=1,
iterhist=_iterhist, finalparms=_finalparms, odsselect=all, deletework=1, keep=);

* by: by-processing variable (e.g. for simulation). must not be missing (.).
* odsselect: none= suppresses all output and leaves ods select none on exit,   all= sets ods select all on exit;
* nonotes: 1 if options nonotes is active, 0 else (to ensure nonotes is activated on exit);
* keep: any extra variables that should be carried over to the data set with predictions;

%let iter=1;
%let conv=0;
%let maxdbeta=.;

%let nvar=0;
%do %while(%scan(&varlist,&nvar+1)~=);
 %let nvar=%eval(&nvar+1);
 %let var&nvar=%scan(&varlist,&nvar);
%end;


data _work;
set &data;
_hdiag_=0.01;
run;

data &iterhist;
run;
data _workfinished;
run;

ods select none;    * turns off output temporarily for the iterations;

%do %while (&conv = 0 & &iter <= &maxiter);
	data _work (drop=_hdiag_);
	set _work;
	y_mod = &y + _hdiag_/2;
	run;

	proc genmod data=_work;
	ods output ParameterEstimates=_parms;
	model y_mod = &varlist /dist=&dist 
		%if &offset ne %then %do; offset=&offset %end;  ;
	output out=_work h=_hdiag_;
	%if &by ne %then %do; 
		by &by;
	%end;
	run;
	quit;

	data _parms;
	set _parms;
	keep estimate &by;
	run;

	proc transpose data=_parms out=_beta_t;
	%if &by ne %then %do; 
		by &by;
	%end;
	run;

	data _beta_t;
	set _beta_t;
	_iter_ = &iter;
	run;


	data &iterhist;
	%if &iter=1 %then %do; 
		set _beta_t;
	%end;
	%else %do;
		set &iterhist _beta_t;
	%end;
	run;

	proc sort data=&iterhist out=&iterhist;
	by &by _iter_;
	run;

	data &iterhist;
	set &iterhist;
	_dbeta_= abs(Col1 - lag(col1)) 
			%do j=2 %to &nvar+1;
				+ abs(col&j - lag(col&j))
			%end; ;
	if _iter_=1 then _dbeta_=9999;
	run;

	%if &iter > 1 %then %do;
		data _convcheck;
		set &iterhist;
		if _iter_ = &iter;
		run;
		proc means data=_convcheck;
		var _dbeta_;
		output out=_convsum max=_maxdbeta_ N=_ndata_;
		run;
	 	data _convsum;
		set _convsum;
		call symput("MAXDBETA", _maxdbeta_);
		call symput("NDATA", _ndata_);
		run;
		%if &by ne %then %do;
			data _finished;
			merge _work _convcheck;
			by &by;
			run;
			data _finished _work(keep=&by &varlist &keep &y &offset _hdiag_);
			set _finished;
			if _dbeta_<&xconv then output _finished;
			else output _work;
			run;
			
			data _workfinished;
			set _workfinished _finished;
			if &by ne .;
			run;

		%end;
		%if %sysevalf(&maxdbeta < &XCONV) = 1 %then %do; 
			%let conv=1; 
			%if &by = %then %do;
				data _workfinished;
				set _work;
				run;
			%end;
		%end;
		%if &nonotes=1 %then %do; 
			options notes;             * temporarily turn on notes to put iteration info into log;
		%end;
		%if &iterinfo=1 %then %do;
			%put NOTE: ITERATION: &iter  DATASETS: &ndata    MAXDBETA: &maxdbeta;  
		%end;
		%if &nonotes=1 %then %do;
			options nonotes;           * turn off again;
		%end;
	%end;

	%let iter= %sysevalf(&iter + 1);

%end;

%if &conv = 0 %then %do;
	%if &nonotes=1 %then %do; 
		options notes;             * temporarily turn on notes to put iteration info into log;
	%end;
	%if &iterinfo=1 %then %do;
		%put WARNING: Iterations not converged for &ndata datasets. MAXDBETA: &maxdbeta;  
	%end;
	%if &nonotes=1 %then %do;
		options nonotes;           * turn off again;
	%end;
	data _work (drop=_hdiag_);
	set _work;
	y_mod = &y + _hdiag_/2;
	run;

	data _workfinished;
	set _workfinished _work;
	run;
%end;


%if &by ne %then %do;
		proc sort data=_workfinished out=_workfinished;
		by &by;
		run;
%end;

data &iterhist;
set &iterhist;
rename col1=Intercept; 
	%do j=1 %to &nvar; 
		%let jj=%sysevalf(&j+1);
		rename col&jj = &&var&j;
	%end;
%let jj=%sysevalf(&nvar + 2);
rename col&jj = Dispersion;
run;

proc means data=&iterhist noprint;
var _iter_;
output out=_maxit max=_maxiter_;
%if &by ne %then %do;
	by &by;
%end;
run;

data &finalparms;
merge &iterhist _maxit;
%if &by ne %then %do;
	by &by;
%end;
if _maxiter_=_iter_;
drop _maxiter_ _type_ _freq_;
run;



%if &print=1 %then %do;
	ods select all;
%end;

proc genmod data=_workfinished;
ods output ParameterEstimates=_FirthParms;
title3 "Firth model";
model y_mod = &varlist  /dist=&dist 
	%if &pl=1 %then %do; lrci %end;
	%if &offset ne %then %do; offset=&offset %end; ;
output out=_workfinished h=_hdiag2_ Predicted=_FIRTHPred;
%if &by ne %then %do; by &by; %end;
run;

data _FIRTHPredictions;
set _workfinished;
run;


data _workfinished;
set _workfinished;
drop _FIRTHPred;
_added_covariate_=0;
output;
&y = _hdiag2_/2;
_added_covariate_=1;
output;
run;

proc genmod data=_workfinished;
ods output ParameterEstimates=_FLACParms;
title3 "FLAC model";
model &y = &varlist _added_covariate_ /dist=&dist 
	%if &pl=1 %then %do; lrci %end;
	%if &offset ne %then %do; offset=&offset %end; ;
output out=_FLACPredictions Predicted=_FLACPred;
%if &by ne %then %do; by &by; %end;
run;

data _FLACPredictions;
set _FLACPredictions;
if _added_covariate_=0;
run;

%if &deletework = 1 %then %do;
	data _workfinished;
	run;
%end;

title3;

%if %upcase(&odsselect) = NONE %then %do;
	ods select none;
%end;
%if %upcase(&odsselect) = ALL %then %do;
	ods select all;
%end;

%mend;
