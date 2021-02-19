%macro dataugpoisson(data=, y=, offset=, varlist=, dist=poisson, priorinterval=1000, outmodel=_DATAUGMODEL, finalparms=_DATAUGParms, by=, S=10000,
procoptions=, modeloptions=, printprior=0, id=, print=1, odsselect=all);

* implements Bayesian data augmentation (ridge regression, or Bayesian analysis with zero-centered normal prior) for Poisson regression with fixed penalty parameter by data augmentation and approximation;
* pseudo data is added which represents the Bayesian weakly informative priors centered around 0 for each parameter;
* penalty parameter can be controlled by specifying the upper limit of a 95% prior interval for the IRR for each variable (option priorinterval);

* data 			... data set to be analyzed. If no data set is specified, macro will just compute prior variances;
* y 			... outcome variable;
* offset		... offset variable (log of rate multiplier) - if not specified the macro will create a variable _offset_ which is constantly 0;
* varlist		... explanatory variables (blank-separated);
* dist			... distribution, normally poisson or negbin;
* priorinterval ... upper limit of prior interval for exp(beta) (see above), one value for each variable. Must be integer.
					Macro will use it to specify prior variances v_j and penalty parameters 1/v_j.
					If one value is given, it applies to all variables. If several values are given, they correspond to the variables. If fewer values aregiven than there are variables,
 					the program will repeat the last given one for all further variables;
* outmodel		... model file to keep the output model (not yet implemented);
* finalparms	...	data set with parameter estimates and 95%CI (as produced by ods output ParameterEstimates;
* by			... grouping variable to efficiently analyze multiple similar-structured data sets;
* S				... approximation constant, should be >= 25 according to Sullivan and Greenland or better >=10000 according to own experience. can be checked by inspecting
					symmetry of prior intervals by option printprior;
* procoptions	...	further options to be passed to PROC GENMOD statement;
* modeloptions	... further options to be passed to MODEL statement;
* printprior 	... =1 will print the analysis of the prior (pseudo) data;
* id 			... variables that should be transferred into the output data set _DATAUGPredictions with predicted values;

%let nvar=0;
%do %while(%scan(&varlist,&nvar+1)~=);
 %let nvar=%eval(&nvar+1);
 %let var&nvar=%scan(&varlist,&nvar);
 %if %scan(&priorinterval, &nvar)~= %then %do; %let pi_&nvar=%scan(&priorinterval, &nvar); %let pi_work=&nvar; %end;
 %else %do; %let pi_&nvar=%scan(&priorinterval, &pi_work); %end;
%end;


data _prior;
length variable $ 30;
%do j=1 %to &nvar; 
	Variable="&&var&j";
	prior_beta_up_&j=log(&&pi_&j);
	prior_se_&j=prior_beta_up_&j/1.96;
	prior_v_&j=prior_se_&j**2;
	call symput("prior_v_&j", prior_v_&j);
	UpperLimit_IRR = &&pi_&j;
	UpperLimit_beta = log(UpperLimit_IRR);
	StdErr = prior_se_&j;
	Variance = prior_v_&j;
	Penalty = 1 / prior_v_&j;
	output;
%end;
run;


%if &printprior=1 %then %do;
	proc print data=_prior;
	var Variable UpperLimit_IRR UpperLimit_beta StdErr Variance Penalty;
	run;
%end;

ods select none;

%if &data ne %then %do;
	data _work;
	set &data;
	_const_=1;
	%if &offset= %then %do;
		_offset_=0;
		%let offset=_offset_;
	%end;
	%if &by = %then %do;
		%let byproc=0;
		%let by=_by_;
		_by_=1;
	%end;
	%else %let byproc=1;
	run;

	proc sort data=_work out=_work;
	by &by;
	run;


	data _aug;
	set _work;
	by &by;
	if first.&by;
	keep &by;
	run;

	data _aug;
	set _aug;
	_const_=0;
	%do j=1 %to &nvar;
		&&var&j=0;
	%end;
	%do j=1 %to &nvar;
		&&var&j = 1/&S;
		&y = &S **2/&&prior_v_&j;
		&offset = log(&y);
		output;
		&&var&j = 0;
	%end;
	run;

	data _work2;
	set _work _aug;
	keep &id _const_ &varlist &y &offset &by;
	run;

	proc sort data=_work2;
	by &by descending _const_;
	run;

	%if &print=1 %then %do; 
		ods select all; 
	%end;

	%if &printprior=1 %then %do;
		proc genmod data=_work2 &procoptions;
		title3 "Prior data";
		model &y =  &varlist /dist=&dist noint lrci offset=&offset &modeloptions;
		where _const_=0;
		%if &byproc = 1 %then %do; by &by; %end;
		run;
	%end;

	ods output ParameterEstimates=&finalparms;
	proc genmod data=_work2 &procoptions;
	title3 "Augmented data";
	model &y = _const_ &varlist /dist=&dist noint lrci offset=&offset &modeloptions;
	output out=_DATAUGPredictions Predicted=_DATAUGPred;
	%if &byproc = 1 %then %do; by &by; %end;
	run;

	data _dataugPredictions;
	set _dataugPredictions;
	if _const_ ne 0;
	run;

	%if %upcase(&odsselect) = NONE %then %do;
		ods select none;
	%end;
	%if %upcase(&odsselect) = ALL %then %do;
		ods select all;
	%end;
%end;
%else %do;
	%put NOTE: no data set specified - macro stops;
%end;



%mend;

