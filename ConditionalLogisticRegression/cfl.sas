%macro cfl(data=, varlist=, covar=, strata=, y=, maxit=50, epsilon=0.001, conflev=0.95,  by=,
            outtab=_res_, print=1, indist=, insuff=, rf=0.001, optn=0, lownumber=1e-7);

*** version 120716;
*** last changes:
	p-value computation by PLR tests implemented (120716);

*** computes Firth-corrected conditional logistic regression model;
*** written by Georg Heinze, 2008-12;		

*** outputs (log) odds ratio estimates and 95% profile likelihood confidence intervals ;

*** uses PROC IML/NLPNRA and NLPLM routines;

*** Input:

	data ...		data set for analysis
	varlist ...		list of variables for which estimates are computed
	covar ...		list of	continuous covariates for which no estimates are needed. These covariates will be conditioned out of likelihood.
    strata ...		any stratifying variables (nominal, e. g. clusters or matched sets). Will be conditioned out of likelihood.
    y ...			dependent binary (0/1) variable
	conflev ...		confidence level for confidence intervals
	by ...			BY variables (for simulation)
	indist ...		input data set of distribution of sufficient statistic, as output by outdist option of proc logistic/exact statement
					(specify to save computing time, if such a distribution has already been computed)
					if this distribution has been computed by a previous call of the macro (with the same data set but different modeling options)
					then you can access it by specifying indist=__indist_
	insuff ...		input data set with observed values of sufficient statistics 
			 		(specify to save computing time, if it has already been computed)
					if this data set has been computed by a previous call of the macro (with the same data set, but with different modeling options)
			        then you can access it by specifying insuff=__insuff_

*** Options controlling iteration:

	maxit ...		maximum number of iterations
	epsilon ...		maximum allowed change in parameter estimates to stop iteration
	rf ...			rounding factor for distribution of sufficient statistic
	lownumber ...	limit for computation of log
	
*** Output:
	optn ...		control output of optimization routines (see SAS/IML manual). Values > 0 will allow to retrace iteration.
	print ...		print output (1=yes, 0=no)
	outtab ...		SAS data set with results
;	



%let nvar=0;
%do %while(%scan(&varlist,&nvar+1)~=);
 %let nvar=%eval(&nvar+1);
 %let var&nvar=%scan(&varlist,&nvar);
%end;


data _work;
set &data;
%if &by= %then %do;
 _rby_=1;
 %let maxrby=1;
 run;
%end;
%else %do;
 run;
 proc freq data=_work;
 tables &by/noprint out=_rby_;
 run;
 data _rby_;
 set _rby_;
 _rby_=_n_;
 call symput("maxrby",_rby_);
 run;
 data _work;
 merge _work _rby_;
 by &by;
 run;
%end;




%if &indist = %then %do;
 ods listing close;
 ods output suffstats=__insuff_ exactparmest=__xparmest;
 proc logistic descending data=_work;
 %if &strata ne %then %do; strata &strata; %end;
 model &y=&varlist &covar;
 exact &varlist/outdist=__indist_ jointonly;
 by _rby_;
 run;
 ods listing;
%end;
%else %do;
 %if &by ne %then %do;
  data __insuff_;
  merge &insuff _rby_;
  by &by;
  run;
  data __indist_;
  merge &indist _rby_;
  by &by;
  run;
 %end;
 %else %do;
  data __insuff_;
  set &insuff;
  _rby_=1;
  run;
  data __indist_;
  set &indist;
  _rby_=1;
  run;
 %end;
%end;

*proc print data=__insuff_;
*run;

data _dist_;
set __indist_;
%do j=1 %to &nvar;
 _t&j=&&var&j;
%end;
_match_=1;
run;


proc means noprint data=_dist_;
var %do j=1 %to &nvar; _t&j %end; ;
by _rby_ _match_;
output out=_msd_ mean= std=/autoname;
run;


data _dist_;
merge _dist_ _msd_;
by _rby_ _match_;
%do j=1 %to &nvar;
 _t&j=((_t&j-_t&j._mean));
%end;
run;




data __insuff_;
set __insuff_;
if 1=0 %do j=1 %to &nvar; or parameter="&&var&j" or parameter="%upcase(&&var&j)" %end;;
run;



proc iml;
start like (beta) global(tobs, tc, grad, hess);
*    cobs=tc[loc(tc[,1:&nvar]=t(tobs)),&nvar+1];
*	tc1=tc[loc(tc[,1:&nvar]=t(tobs)),];;
*	print cobs tc1;
*    put "start like";
	S0 = 0;
	SR = repeat(0,&nvar,1);
	SRS = repeat(0, &nvar, &nvar);
	%do j=1 %to &nvar;
	 SRS&j=repeat(0,&nvar,&nvar);
	 DABL&j=repeat(0,&nvar,&nvar);
	%end;
	*** compute sums;
	like= tobs`*beta`;
*	put "noch da (vor i schleife";
	do i=1 to nrow(tc);
*	 if tc[i,1:&nvar]=t(tobs) then cobs=tc[i,&nvar+1];
	 score=tc[i,ncol(tc)]*exp(tc[i,1:&nvar]*beta`);
     *print i score; 
	 S0 = s0 + score;
	 do j=1 to &nvar;
	  SR[j] = SR[j] + tc[i,j]* score;
	  do jj=1 to &nvar;
	   SRS[j,jj] = SRS[j,jj] + tc[i,j]*tc[i,jj] * score;
	   %do j=1 %to &nvar;
	    SRS&j[j,jj] = SRS&j[j,jj] +  tc[i,j]*tc[i,jj]*tc[i,&j] * score;
	   %end;
	  end;
	 end;
	end;
*	put "noch da (nach i schleife)";
*	if cobs=0 then cobs=1;
*	like =log(cobs)+ tobs`*beta` - log(S0); * cobs does not depend on beta, so it can be omitted as a constant;
	if s0<=0 then s0=&lownumber;
    like =tobs`*beta` - log(S0);
	grad=t(tobs - SR/S0);
	hess=-(SRS/S0-SR*SR`/(S0**2));
	dethess=det(-hess);
	if dethess<=0 then dethess=&lownumber;
    like = like + 0.5 * log(dethess);
	do j=1 to &nvar;
	 do jj=1 to &nvar;
	  %do j=1 %to &nvar;
	   DABL&j[j,jj] = SRS&j[j,jj]/S0 - SRS[j,jj]*SR[&j]/S0/S0 - SR[jj]/s0*(SRS[j,&j]/S0 - SR[j]*SR[&j]/S0/S0) - SR[j]/S0*(SRS[jj,&j]/s0 - SR[jj]*SR[&j]/S0/S0);
	  %end;
	 end;
	end;
	dinfoi=inv(-hess);
	%DO j=1 %to &nvar;
     TRACE=0;
     DO jj=1 to &nvar;
      TRACE=TRACE+DINFOI[jj,]*DABL&j[,jj];
	 end;
     grad[&j]=grad[&j]+trace/2;
    %end;
	return(like);
finish;

start f_grad(beta) global(tobs, tc, grad, hess);;
 return(grad);
finish;

start f_hess(beta) global(tobs, tc, grad, hess);;
 return(hess);
finish;


  
start f_pl(x) global(tobs, tc, ipar,lstar); 
       /* x[1]=sigma, x[2]=c */ 
       lik = like(x); 
       grad = f_grad(x); 
       grad[ipar] = lik - lstar; 
       return(grad`); 
finish f_pl;





use _dist_;
read all var (%do j=1 %to &nvar; "_t&j" || %end;  "count" || "_rby_") into tc_all;
close _dist_;

use _msd_;
read all var ("_t1_mean" %if &nvar>1 %then %do; %do j=2 %to &nvar; || "_t&j._mean" %end; %end; || "_rby_") into tmeans_all;
close _msd_;

use __insuff_;
read all var("value" || "_rby_") into tobs_all;
close __insuff_;


do iby = 1 to &maxrby;

tc=tc_all[loc(tc_all[,ncol(tc_all)]=iby),1:(ncol(tc_all)-1)];
tmeans=tmeans_all[loc(tmeans_all[,ncol(tmeans_all)]=iby),1:(ncol(tmeans_all)-1)];
tobs=tobs_all[loc(tobs_all[,ncol(tobs_all)]=iby),1:(ncol(tobs_all)-1)];


tobs=tobs-t(tmeans);

tobs=round(tobs,&rf);
tc[,1:&nvar]=round(tc[,1:&nvar],&rf);
tc[,&nvar+1]=tc[,&nvar+1]/tc[><,&nvar+1];

npar=&nvar;

file log;

put "NOTE: DATA set: " iby;


beta=repeat(0,1,&nvar);
optn={1, &optn};


call nlpNRA(rc,betares,"LIKE",beta, optn) grd="F_GRAD" hes="F_HESS"; 

maxlike=like(betares);
betahat=betares;
pval=repeat(-1,&nvar,1);
chi2=repeat(-1,&nvar,1);

*** calculate p-values;
con=repeat(.,3,&nvar+2);
con2=con;
x0=betahat;
do j=1 to &nvar;
* x0=betahat;
* x0[j]=0;
 x0=repeat(0,1,&nvar);
 con2[3,]=0;		*reset all constraints;
 con2[3,j]=1;		*specifies the parameter to be subjected to constraints;	
 con2[3,&nvar+1]=0; *specifies equality constraint;
 con2[3,&nvar+2]=0; *specifies the value at which the parameter is constraint;


* print j,x0,con2,betares,optn,rc;

 
 call nlpNRA(rc,betares,"LIKE",x0,optn,con2) grd="F_GRAD" hes="F_HESS"; 
 *print j, betares;
 l0=like(betares);
 Chi2[j]=2*(maxlike-l0);
 pval[j]=1-probchi(chi2[j],1);

end;


betares=betahat;
prob= 0.05;

chqua = cinv(1-prob,1); 
lstar = maxlike - .5 * chqua; 
*print chqua lstar;
optn = {2 0}; 
hes2=hess;
xopt=betares`;
xlb=repeat(.,&nvar,1);
xub=xlb;
con=repeat(.,2,&nvar);
optn={&nvar, &optn};
   do ipar = 1 to npar;
     if npar>1 then do;
	  ind=setdif(1:npar,ipar);
    /* Compute initial step: */ 
    /* Choose (alfa,delt) to go in right direction */ 
    /* Venzon & Moolgavkar (1988), p.89 */ 
*       if ipar=1 then ind = 2;
*		else ind = 1; 
       delt = - inv(hes2[ind,ind]) * hes2[ind,ipar]; 
       alfa = - (hes2[ipar,ipar] - delt` * hes2[ind,ipar]); 
       if alfa > 0 then alfa = .5 * sqrt(chqua / alfa); 
       else do; 
          print "Bad alpha"; 
          alfa = .1 * xopt[ipar]; 
       end; 
       delt = 1 || delt`; 
       indices= ipar || ind;
	   ord = rank(indices);
       delt[ord] = delt[1:npar];
   	  end; 
	  else do;
	   delt=1;
	   alfa = .1* xopt[ipar];
	  end;
	  alfadelt=(alfa * delt)`;
	  if alfadelt[ipar]<=1e4 then alfadelt[ipar]=0.1;
  
    /* Get upper end of interval */ 
       x0 = xopt + alfadelt; 
	   
    /* set lower bound to optimal value */ 
       con2 = con; con2[1,ipar] = xopt[ipar]; 
	   f = f_pl(xopt`);
	   *print f;
       call nlplm(rc,betares,"f_pl",x0,optn,con2); 
       f = f_pl(betares); s = ssq(f); 
       if (s < 1.e-6) then xub[ipar] = betares[ipar]; 
          else xub[ipar] = .; 
  
    /* Get lower end of interval */ 
       x0 = xopt - alfadelt; 
    /* reset lower bound and set upper bound to optimal value */ 
       con2[1,ipar] = con[1,ipar]; con2[2,ipar] = xopt[ipar]; 
       call nlplm(rc,betares,"f_pl",x0,optn,con2); 
       f = f_pl(betares); s = ssq(f); 
       if (s < 1.e-6) then xlb[ipar] = betares[ipar]; 
          else xlb[ipar] = .; 
    end; 
	if iby=1 then res=repeat(1,&nvar,1)||t(1:&nvar)||xopt||xlb||xub||chi2||pval;
	else res=res//(repeat(iby,&nvar,1)||t(1:&nvar)||xopt||xlb||xub||chi2||pval);
    *print "Profile-Likelihood Confidence Interval"; 
    *print xlb xopt xub;
end;

create _res_ from res[colname={"_rby_" "_var_" "Estimate" "LowerCL" "UpperCL" "Chi2" "Pvalue"}];
append from res;
close _res_;

quit;

data &outtab;
%if &by ne %then %do;
 merge _res_ _rby_;
 by _rby_;
%end;
%else %do;
 set _res_;
%end;
OddsRatio=exp(Estimate);
ORLowerCL=exp(LowerCL);
ORUpperCL=exp(UpperCL);
%do j=1 %to &nvar;
 if _var_=&j then VariableName="&&var&j";
%end;
run;

%if &print=1 %then %do;
 proc print;
 title2 "Conditional logistic regression";
 title3 "Firth-corrected Parameter estimates and Profile likelihood confidence intervals";
 var VariableName Estimate LowerCL UpperCL OddsRatio ORLowerCL ORUpperCL Chi2 pvalue;
 %if &by ne %then %do;
  by &by;
  pageby &by;
 %end;
 run;
%end;


title2;

%mend;
