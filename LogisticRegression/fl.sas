%macro fl(data=, y=, varlist=, maxit=50,
 epsilon=0.0001, noint=0, by=' ', outtab=_outtab, outest=_outest, print=1,
 offset=, maxstep=5, maxhs=5, pl=1, plint=0, plmaxit=50, alpha=0.05,
 odds=0, test=, outlr=_test, out=, pred=_pred_, lower=_lo_, upper=_up_, h=_h_, 
 outprof=_prof, profile=, profsel=, profser=, profn=100, notes=1, global=_globtest, standard=1, pl1=0, ortho=1, huge=1000, where=%str(1=1));

%let version=2009.06;
%let build=231505;

%let critpi=0.00000001;

* Georg Heinze Nov 1999 - Jan 2009                                ;

* fits a logistic regression model using Firth`s bias reduction method (modified score function);

* Please find the complete documentation of this program in a TechRep at
  http://www.meduniwien.ac.at/user/georg.heinze/techreps/tr2_2004.pdf;

* log of changes: 						;
*	2009-01-28:		in computation of pl confidence limits, iteration stops if absolute values of 
					standardized estimates exceed 100. The program will write a NOTE into the log, but proceed.;


* additional macro options not contained in the TechRep:
    where       ... define a subset of the data set (e.g., where=%str(age > 49 & sex = 0)) (added 2009-06-23)

	pl1 		... Profile penalized likelihood confidence interval only for first variable in varlist (no=0=default, yes=1)
					you can use this option to save computing time if a confidence interval is needed for one variable only.

	ortho		... Orthogonalization of covariates before estimation of profile penalized likelihood
					confidence intervals (yes=1=default value, no=0)
					this option offers more robustness in numerical computations at the cost of a potential increase in computing time
                    in most cases total computing time will be less

    huge        ... controls the memory allocation in the PROC IML parts of the macro. The options bypasses the
			 		allocation of n x n matrices which is usually done by the macro to compute the H matrix by
					computing the H matrix diagonal elements elementwise (slow but safe). There are three possible values:
					1  			 ... bypass allocation of n x n matrices
					0  			 ... always allocate n x n matrices (faster with small data sets)
				    any value >1 ... automatic decision (bypass turned on if N>&huge)
					default value... 1000 (meaning that bypass will be turned on if N>1000);


%let maxdiff=%sysevalf(&maxstep);
%if &noint=1 %then %do;
	%let plint=%eval(&pl);
	%if &standard=1 %then %do;
	 %put NOTE: Standardization turned off because of noint=1 option.;
	 %let standard=0;
	%end;
%end;

%let nvar=0;
%do %while(%scan(&varlist,&nvar+1)~=);
 %let nvar=%eval(&nvar+1);
 %let var&nvar=%scan(&varlist,&nvar);
%end;

%let ntest=0;
%do %while(%scan(&test,&ntest+1)~=);
 %let ntest=%eval(&ntest+1);
 %let test&ntest=%scan(&test,&ntest);
%end;


%let conflev=%sysevalf(100-100*&alpha);

%let by=%upcase(&by);
%let nby=0;
%if &by ne ' ' %then %do;
 %do %while(%scan(&by,&nby+1)~=);
  %let nby=%eval(&nby+1);
  %let by&nby=%scan(&by,&nby);
 %end;
%end;


%let nprof=0;
%do %while(%scan(&profile,&nprof+1)~=);
 %let nprof=%eval(&nprof+1);
 %let prof&nprof=%scan(&profile,&nprof);
%end;

%if &profn <3 %then %let profn=3;


%if &by=' ' %then %let byproc=0;
%else %let byproc=1;

%if &out~= %then %do;
 %if &outest=' ' %then %do;
  %let outest=_outest;
 %end;
%end;


data _work;
set &data;
if &where;
intercep=1;
%let var0=INTERCEP;
if &y ne .;
%do k=1 %to &nvar;
 if &&var&k ne .;
%end;
%if &by=' ' %then %do;
 %let by=_by_;
 _by_=1;
 %let nby=1;
 %let by1=_by_;
%end;
run;

proc sort;
by &by;
run;

data _fby;
set _work;
by &by;
if first.&&by&nby;
keep &by;
run;

data _fby;
set _fby;
_rby_=_n_;
call symput("maxby", _rby_);
run;


data _work;
merge _work _fby;
by &by;
run;


* rby is now column with numbers from 1 to #by-groups ;

%if &offset~= %then %do;
 %if &byproc=1 %then %do;
  proc sort data=&offset out=&offset;
  by &by;
  run;
 %end;
 data _offset_;
 merge &offset _fby;
 %if &byproc=1 %then %do;
  by &by;
 %end;
 run;
%end;
%else %do;
 data _offset_;
 set _fby;
  %do j=0 %to &nvar;
   &&var&j=.;
  %end;
 run;
%end;

data _estwot_;
set _offset_;
%do j=&noint %to &nvar;
 if &&var&j = . then do;
  _hilf_=1;
 end;
 else do;
  _hilf_=0;
 end;
 &&var&j=_hilf_;
%end;
drop _hilf_;
run;

* avoid missings in _offset_ data set ;

data _offset_;
set _offset_;
%do j=&noint %to &nvar;
 if &&var&j=. then &&var&j=0;
%end;
run;


%if &noint=0 %then %do;
 %let nvar=%eval(&nvar+1);
%end;


proc iml;

*******************************************************************;
start pl_comp;
 errcode=0;
 linpred=x*beta;
 pi=1/(1+exp(-linpred));
 logpi=pi;
 wvec=pi#(1-pi);
 loglike=0;
 rtwvec=sqrt(wvec);
 rtwx=x#rtwvec;
 Fisher=t(rtwx)*rtwx;
 detfish=det(Fisher);
 if detfish=0 & comefromprofile=1 then do;
  errcode=1;
  put "ERROR (numerical): Profile likelihood cannot be determined.";
 end;
 else do;
  if detfish=0 & comefromprofile=0 then do;
   errcode=2;
   put "ERROR (numerical): Likelihood cannot be determined.";
  end;
  else do;
   if huge=0 | (huge>1 & n<=huge) then do;
    * for huge data sets, the term in brackets can be too large (nxn Matrix);
    hvec=vecdiag(rtwx*inv(Fisher)*t(rtwx));
	hvec_pi=diag(hvec)*pi;
   end;
   else do;
    * workaround for huge data sets;
    help1=rtwx*inv(Fisher);
    do i=1 to n;
     hvec[i]=help1[i,]*t(rtwx[i,]);
 	 hvec_pi[i]=hvec[i]*pi[i];
    end; 
   end;
   logpi[loc(y=1)]=log(pi[loc(y=1)]);
   logpi[loc(y=0)]=log(1-pi[loc(y=0)]);
   loglike=logpi[+];

   if detfish<&critpi then do;
    file log;
    put "WARNING: det(Fisher)=0, switched to slow but safe mode";
    loglike=0;
    do i=1 to n;
     pi[i]=1/(1+exp(-x[i,]*beta));
     wvec[i]=pi[i]*(1-pi[i]);
     rtwvec[i]=sqrt(wvec[i]);
    end;
    rtwx=repeat(0,n,&nvar);
    do i=1 to n;
     rtwx[i,]=x[i,]*rtwvec[i];
    end;
    Fisher=t(rtwx)*rtwx;


    ifish=inv(Fisher);
    help1=rtwx*ifish;
    do i=1 to n;
     hvec[i]=help1[i,]*t(rtwx[i,]);
	 hvec_pi[i]=hvec[i]*pi[i];
    end; 

    do i=1 to n;
     if y[i]=1 then loglike=loglike+log(pi[i]);
     else loglike=loglike+log(1-pi[i]);
    end;
    penlike=loglike+0.5*log(det(Fisher));
    g=t(x)*(y+0.5*hvec-pi-hvec_pi);
    v=-Fisher;
    Cov=ifish;
   end;
   else do;
    penlike=loglike+0.5*log(detfish);
	* huge!;
    g=t(x)*(y+0.5*hvec-pi-hvec_pi);
    v=-Fisher;
    Cov=inv(Fisher);
   end;
  end;
 end;
finish pl_comp;


*******************************************************************;


**********************************************************************;
start fl03;

 * uses:    x, y, ew, os, n, comefromprofile, &nvar, &noint, and the subroutine pl_comp;
 * output:  beta, loglike, penlike, vm, U, it;

 nvar=&nvar;
 xsave=x;
 if &standard=1 then do;
  xmean=repeat(x[:,],n,1);
  xstd=repeat(((x#x)[:,]-x[:,]#x[:,])##0.5,n,1);
 end;
 else do;
  xmean=repeat(0,n,nvar);
  xstd=repeat(1,n,nvar);
 end;
 if &noint=0 then do;
  xmean[,1]=repeat(0,n,1);
  xstd[,1]=repeat(1,n,1);
 end;
* xmean=repeat(xmean,n,1);
 file log;

 



 xs=x;
 xs=(x-xmean)/xstd;
 x=xs;
 if (ew[1]=0 & &noint=0) then
  x=xs+xmean/xstd;
  beta=os;
 if &noint=1 then
  beta=os#t(xstd[1,]);
 else do;
  beta[2:&nvar]=os[2:&nvar]#t(xstd[1,2:&nvar]);
  beta[1]=os[1]+xmean[1,2:&nvar]*os[2:&nvar];
 end;

 dimx=ew[+];

 if dimx=0 then do;
  run pl_comp;
 end;
 else do;
  xwot=repeat(0,n,dimx);
  map=repeat(0,dimx,1);   * to map x to xwot ;

  jj=0;
  do j=1 to nvar;
   if ew[j]=1 then do;
    jj=jj+1;
    map[jj]=j;
   end;
  end;


  do j=1 to dimx;
   xwot[,j]=x[,map[j]];
  end;



  if &noint=0 & ew[1]=1 & comefromprofile=0 then do;
   eta=x*beta;
   eta_bar=eta[+]/n;
   y_bar=y[+]/n;
   beta[1]=log(y_bar/(1-y_bar))-eta_bar;
  end;


  pi=repeat(0.5,n,1);
  wvec=repeat(0.25,n,1);
  rtwvec=repeat(0.5,n,1);
  rtpi=repeat(0.5,n,1);
  rtwxwot=repeat(0,n,ncol(xwot));
  di=repeat(1,&nvar,1);
  hvec=repeat(0,n,1);
  hvec_pi=repeat(0,n,1);
  diff=1;
  it=0;

  run pl_comp;
  if errcode=2 then put "ERROR: Numerical problem in computation of likelihood.";
  else if errcode =0 then do;
   do while(it<&maxit & diff>abs(&epsilon));
    it=it+1;
    oldbeta=beta;
    rtwxwot=xwot#rtwvec;
    vmwot=ginv(t(rtwxwot)*rtwxwot);
    vm=repeat(0,&nvar,&nvar);
    do j=1 to dimx;
     do jj=1 to dimx;
      vm[map[j],map[jj]]=vmwot[j,jj];
     end;
    end;
	*** huge ! ;
    fd=t(x)*(y+0.5*hvec-pi-hvec_pi);
    incr=0;
    incr=vm*(fd#ew);
    if any(incr>&maxstep) then incr=sign(incr)*&maxstep;
    beta=beta+incr#ew;

    hs=0;
    oldlike=penlike;
    run pl_comp;
    do while(hs<&maxhs & penlike<oldlike);
     hs=hs+1;
     beta=(beta+oldbeta)/2;
     run pl_comp;
    end;

    U=t(g);
    diff=abs(beta-oldbeta)[+];
  
   end;

   run pl_comp;

   lp=x*beta;
   betastar=beta;
   Covstar=Cov;
   if &noint=0 then do;
    redo=1//t(1/xstd[1,2:&nvar]);
*    beta[2:&nvar]=beta[2:&nvar]#t(1/xstd[1,2:&nvar]);
    beta=beta#redo;
	Cov=Covstar#(redo*redo`);
    beta[1]=beta[1]-xmean[1,2:&nvar]*beta[2:&nvar];
   end;
   else do;
    redo=t(1/xstd[1,]);
    beta=beta#redo;
	Cov=Covstar#(redo*redo`);
   end;

   x=xsave;
   wx1=(x-repeat(0||xmean[1,2:&nvar],n,1))`*sqrt(pi#(1-pi));   *p*n x n*1 = p*1;
*   wx2=(1-pi)`*x; * 1*n x n*p = 1*p;
*   print Covstar;
*   Cov=inv(wx1*wx1`);
*    Cov={1 0,0 1};
*   Cov=inv(x`*diag(pi#(1-pi))*x);
  end;
 end;
finish fl03;

****************************************************************************;

start fl03pl;

 xsave=x;
 Covsave=Cov;
 Covstarsave=Covstar;
 betastarsave=betastar;

 if &standard=1 then do;
  xmean=repeat(x[:,],n,1);
  xstd=repeat(((x#x)[:,]-x[:,]#x[:,])##0.5,n,1);
 end;
 else do;
  xmean=repeat(0,n,nvar);
  xstd=repeat(1,n,nvar);
 end;

 if &noint=0 then do;
  xmean[,1]=repeat(0,n,1);
  xstd[,1]=repeat(1,n,1);
 end;
 file log;


 xs=x;
 xs=(x-xmean)/xstd;
 %if &notes=1 %then %do;
  put "NOTE: Likelihood evaluations at computation of pl confidence intervals:";
  put "      Variable     Limit Iteration Maxstep = &maxdiff exceeded";
 %end;
 %if &pl1=1 %then %do;
  do jjj=2-&plint to 2-&plint;
 %end;
 %else %do;
  do jjj=2-&plint to &nvar;
 %end;
   x=xs;
   if &noint=0 & jjj=1 then do;
    x=xs+xmean/xstd;
   end;
   %if &ortho=1 and jjj>1-&noint %then %do;
	xest=x[,jjj];
	xleft=repeat(1,nrow(xest),1);
	xright=xleft;
	if jjj>1 then do;
	 xleft=xleft||x[,1:jjj-1];
	end;
	if jjj<&nvar then do;
	 xright=x[,jjj+1:&nvar]||xright;
	end;
	if ncol(xleft)>1 then do;
     if ncol(xright)>1 then do;
      xdes=xleft[,2:ncol(xleft)]||xright[,1:ncol(xright)-1];
	 end;
	 else xdes=xleft[,2:ncol(xleft)];
    end;
	else do;
     if ncol(xright)>1 then do;
 	  xdes=xright[,1:ncol(xright)-1];
 	 end;
	end;
	help1=xdes*inv(t(xdes)*xdes);  * n*p*((p*n)*(n*p)) = n*p;
	help2=t(xdes)*xest; * p*n * n*p = p*p;
	x[,jjj]=xest-help1*help2;
   %end;
   do l=-1 to 1 by 2;           /*  -1=lower limit, 1=upper limit  */

*	print l jjj x;
    penlike=pl_save;
    beta=betastarsave;
    ej=repeat(0,&nvar,1);
    ej[jjj]=1;
	run pl_comp;
	it=0;
    rest=penlike-l0;
    do while(it<&plmaxit & abs(penlike-l0)>&epsilon);
	 iv=inv(v);
     wurzel=2*(l0-penlike+0.5*t(g)*iv*g)/(t(ej)*iv*ej);
     if wurzel <0 then do;
      gvg=t(g)*iv*g;
      eve=t(ej)*iv*ej;
      put "ERROR: Numerical problem in computation of profile likelihood c.i.";
     end;
     lambda=l*sqrt(wurzel);
     delta=-iv*(g+lambda*ej);
     oldbeta=beta;
     beta=beta+delta;
     diff=abs(delta)[+];
     hs=0;
     do while(diff>&maxdiff & hs<&maxhs);
      beta=beta-delta;
      delta=delta/(diff/&maxdiff);
     diff=abs(delta)[+];
      beta=beta+delta;
	  %if &notes=1 %then %do;
	   bjds=beta[jjj]/xstd[1,jjj];
	   bjs=beta[jjj];
	   put "     " jjj " " l " " it " yes" " current value (std): " bjs " original: " bjds;
	  %end;
	  run pl_comp;
      hs=hs+1;
     end;
 	 if hs=0 then do;
	  %if &notes=1 %then %do;
	   bjds=beta[jjj]/xstd[1,jjj];
	   bjs=beta[jjj];
	   put "     " jjj " " l " " it " no " " current value (std): " bjs " original: " bjds;
 	  %end;
      run pl_comp;
	 end;
	  if beta[jjj]>100 then do; 
        beta[jjj]=100;
		it=&plmaxit *2;
		put "NOTE: iterations stopped (overflow).";
	  end;
	  if beta[jjj]<-100 then do;
       beta[jjj]=-100;
	   it=&plmaxit *2;
		put "NOTE: iterations stopped (overflow).";
	  end;

     rest=penlike-l0;
     it=it+1;
     bj=beta[jjj];
    end;
	if (&noint=0 & jjj > 1) | (&noint=1) then
     limits[jjj,(l+1)/2+1]=bj/xstd[1,jjj];
	else do;
     limits[jjj,(l+1)/2+1]=bj;
	end;
    achl[jjj,(l+1)/2+1]=penlike;
   end;
  end;
 x=xsave;
 penlike=pl_save;
 beta=b_save;
 Cov=Covsave;
 Covstar=Covstarsave;
 betastar=betastarsave;
finish fl03pl;


**************************************************************************;

use _work;
 read all var("_rby_" || "&by") into by;


%if &noint=0 %then %do;
 read all var("intercep"
  %do k=1 %to &nvar-1;
   ||"&&var&k"
  %end;
  ) into x_all;
%end;
%else %do;
 read all var("&var1"
  %if &nvar>=2 %then %do;
   %do k=2 %to &nvar;
    ||"&&var&k"
   %end;
  %end;
  ) into x_all;
%end;
read all var("&y") into y_all;
close _work;

use _offset_;
read all var (
 %if &noint=1 %then %do;
  "&var1"
 %end;
 %else %do;
  "&var0"
 %end;
 %if &nvar-1+&noint ge &noint+1 %then %do;
  %do j=&noint+1 %to &nvar-1+&noint;
    || "&&var&j"
  %end;
 %end;
 ) into offset_a;
close _offset_;

use _estwot_;
read all var (
 %if &noint=1 %then %do;
  "&var1"
 %end;
 %else %do;
  "&var0"
 %end;
 %if &nvar-1+&noint ge &noint+1 %then %do;
  %do j=&noint+1 %to &nvar-1+&noint;
    || "&&var&j"
  %end;
 %end;
 ) into estwot_a;
close _estwot_;


huge=&huge;

comefromprofile=0;
n_all=nrow(x_all);
maxrby=by[n_all,1];
outtab=repeat(0,maxrby*&nvar,10);
%if &nprof ne 0 %then %do;
 outprof=repeat(0,maxrby*&nprof*&profn,10);
%end;
hvec_all=repeat(0,n_all,1);
%if &test~= %then %do;
 p_spec=repeat(0,maxrby,3);
%end;

%if &outest ne ' ' %then %do;
 outest=repeat(0,(&nvar+1)*maxrby,&nvar+7);
%end;

stop=0;
h_index=0;


do iby=1 to maxrby;

 x=x_all[loc(by[,1]=iby),];
 y=y_all[loc(by[,1]=iby)];
 
 n=nrow(y);

 os=t(offset_a[iby,]);

 Cov=repeat(0,&nvar,&nvar);

 ew=repeat(0,&nvar,1);
 %if &noint=0 %then %do;
  ew[1]=1;
 %end;


  run fl03;
  pl0_save=penlike;


  

 ew=t(estwot_a[iby,]);
 os=t(offset_a[iby,]);

 run fl03;

 %if &noint=0 %then %do;
  ChiSqWald=t(beta[2:&nvar])*inv(Cov[2:&nvar,2:&nvar])*beta[2:&nvar];
  df=&nvar-1;
 %end;
 %else %do;
  ChiSqWald=t(beta)*inv(Cov)*beta;
  df=&nvar;
 %end;
 pvalWald=1-ProbChi(max(0,ChiSqWald),df);
 ChiSqLR=2*(penlike-pl0_save);
 pvalLR=1-ProbChi(max(0,ChiSqLR),df);

 file log;


 *print Cov;
 stderr=repeat(0,&nvar,1);
 ci=repeat(0,&nvar,2);
 p_value=repeat(0,&nvar,1);
 do j=1 to &nvar;
  stderr[j]=sqrt(Cov[j,j]);
  p_value[j]=2*(1-probnorm(abs(beta[j]/stderr[j])));
 end;
 limits=repeat(0,&nvar,2);
 limits[,1]=beta-1.96*stderr;
 limits[,2]=beta+1.96*stderr;

 b_save=beta;
 h_save=hvec;
 Cov_save=Cov;
 pl_save=penlike;
 ll_save=loglike;
 achl=repeat(0,&nvar,2);
 it_save=it;
 pl0=repeat(pl_save,&nvar,1);
 ll0=repeat(ll_save,&nvar,1);


 %if &pl=1 %then %do;

  l0=penlike-0.5*cinv(1-&alpha,1);

  run fl03pl;

   * LR-Tests;

  maxpl=pl_save;
  *os=b_save;
  do k=1 to &nvar;
   os[k]=0;
   ew[k]=0;
   run fl03;
*  print k beta maxpl penlike;
   if maxpl<penlike then do;
    p_value[k]=1.0;
   end;
   else do;
    p_value[k]=1-probchi(max(0,2*(maxpl-penlike)),1);
   end;
   if penlike-maxpl > &epsilon then do;
     print "WARNING: negative log likelihood ratio!", k;
   end;
   os[k]=offset_a[iby,k];
   ew[k]=estwot_a[iby,k];
   beta=b_save;
   pl0[k]=penlike;
   ll0[k]=loglike;
   penlike=pl_save;
  end;

 %end;  * pl;

* Special LR test;

 %if &ntest ne 0 %then %do;

  b_save=beta;
  os=t(offset_a[iby,]);
  %let match=0;
  %do k=0 %to %eval(&nvar-1);
   ew[&k+1]=estwot_a[iby,&k+1];
   %do kk=1 %to %eval(&ntest);
    %if %upcase(&&test&kk)=%upcase(&&var&k) %then %do;
     %let match=%eval(&match+1);
     os[&k+1]=0;
     ew[&k+1]=0;
    %end;
   %end;
  %end;

  %if &match=&ntest %then %do;
   run fl03;
   p_spec[iby,1]=2*(pl_save-penlike);
   p_spec[iby,2]=&ntest;
   p_spec[iby,3]=1-probchi(max(0,2*(pl_save-penlike)),&ntest);

   beta=b_save;
   penlike=pl_save;
   loglike=ll_save;

  %end;
  %else %do;
   %put   ;
   %put ERROR: Some of the variables (%upcase(&test)) are not contained in varlist (%upcase(&varlist)).;
   %put WARNING: No tests were computed.;
   %put NOTE: The above message was produced by the macro FL.;
   %put   ;
   %let ntest=0;
   %let test=;
  %end;

 %end;

 %if &nprof ne 0 %then %do;
  %do j=1 %to &nprof;
   %do jj=1 %to &nvar;
    %let jjdata=%eval(&jj-1);
    %if &&prof&j=&&var&jjdata %then %do;
	 %if &profsel~= %then %do;
	  profsel=&profsel;
	  profser=&profser;
	 %end;
	 %else %do;
	  waldlo=b_save[&jj]-probit(1-&alpha/2)*stderr[&jj];
      waldup=b_save[&jj]+probit(1-&alpha/2)*stderr[&jj];
	  profsel=-(min(limits[&jj,1],waldlo)-b_save[&jj])/stderr[&jj]+0.5;
	  profser=(max(limits[&jj,2],waldup)-b_save[&jj])/stderr[&jj]+0.5;
	 %end;
	 do k=1 to &nvar;
	  os[k]=0;
	  ew[k]=1;
	  *os[k]=offset_a[iby,k];
	  ew[k]=estwot_a[iby,k];
	 end;

     ew[&jj]=0;
	 indp=1;
	 outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+1,7]=limits[&jj,1];
	 outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+2,7]=b_save[&jj];
	 outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+3,7]=limits[&jj,2];

     outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+1,9]=b_save[&jj]-probit(1-&alpha/2)*stderr[&jj];
	 outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+2,9]=b_save[&jj];
	 outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+3,9]=b_save[&jj]+probit(1-&alpha/2)*stderr[&jj];
     
	 minl=0;
	 os=b_save;
	 do bprof=b_save[&jj] to b_save[&jj]-profsel*stderr[&jj] by -(profser+profsel)*stderr[&jj]/(&profn-1);
      os[&jj]=bprof;
      put "NOTE: Computing profile likelihood at beta[&jj]=" bprof;	
	  comefromprofile=1;
	  run fl03;
	  if errcode=1 then put "ERROR: Profile likelihood cannot be determined.";
	  do k=1 to &nvar;
	   os[k]=beta[k];
	  end;
	  indexp=(iby-1)*&nprof*&profn+(&j-1)*&profn+indp;
	  outprof[indexp,1]=iby;
	  outprof[indexp,2]=&jjdata;
	  outprof[indexp,3]=bprof;
	  outprof[indexp,4]=penlike;
	  outprof[indexp,5]=pl_save-0.5*(((bprof-b_save[&jj])/stderr[&jj])**2);
	  if outprof[indexp,4]<minl then minl=outprof[indexp,4];
	  outprof[indexp,6]=pl_save-0.5*cinv(&conflev/100,1);
	  indp=indp+1;
	 end;
	 os=b_save;
	 do bprof=b_save[&jj]+(profser+profsel)*stderr[&jj]/(&profn-1) to b_save[&jj]+profser*stderr[&jj] by (profser+profsel)*stderr[&jj]/(&profn-1);
      os[&jj]=bprof;
	  %if &notes=1 %then %do;
       put "NOTE: Computing profile likelihood at beta[&jj]=" bprof;	
	  %end;
	  comefromprofile=1;
	  run fl03;
	  if errcode=1 then put "ERROR: Profile likelihood cannot be determined.";
	  do k=1 to &nvar;
	   os[k]=beta[k];
	  end;
	  indexp=(iby-1)*&nprof*&profn+(&j-1)*&profn+indp;
	  outprof[indexp,1]=iby;
	  outprof[indexp,2]=&jjdata;
	  outprof[indexp,3]=bprof;
	  outprof[indexp,4]=penlike;
	  outprof[indexp,5]=pl_save-0.5*(((bprof-b_save[&jj])/stderr[&jj])**2);
	  if outprof[indexp,4]<minl then minl=outprof[indexp,4];
	  outprof[indexp,6]=pl_save-0.5*cinv(&conflev/100,1);
	  indp=indp+1;
	 end;

	 outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+1,8]=minl-(pl_save-minl)/10;
	 outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+2,8]=minl-(pl_save-minl)/10;
	 outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+3,8]=minl-(pl_save-minl)/10;

     outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+1,10]=minl-1.5*(pl_save-minl)/10;
	 outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+2,10]=minl-1.5*(pl_save-minl)/10;
	 outprof[(iby-1)*&nprof*&profn+(&j-1)*&profn+3,10]=minl-1.5*(pl_save-minl)/10;
    %end;
    do k=1 to &nvar;
	 os[k]=offset_a[iby,k];
	 ew[k]=offset_a[iby,k];
	end;
   %end;
  %end;
 %end;


 do k=1 to &nvar;
  outtab[(iby-1)*&nvar+k,]=iby||k||b_save[k]||stderr[k]||limits[k,1]||limits[k,2]||
                            p_value[k]||it_save||pl_save||ll_save;
 end;

 do i=h_index+1 to h_index+n;
  hvec_all[loc(by[,1]=iby)]=h_save;
 end;

 if iby=1 then 
  global=(iby||1||ChiSqLR||df||pvalLR)//
  (iby||2||ChiSqWald||df||pvalWald);
 else
  global=global//
  (iby||1||ChiSqLR||df||pvalLR)//
  (iby||2||ChiSqWald||df||pvalWald);
 


 %if &outest ne ' ' %then %do;
  outest[(iby-1)*(&nvar+1)+1,&nvar+1]=-1;
  outest[(iby-1)*(&nvar+1)+1,&nvar+2]=iby;
  outest[(iby-1)*(&nvar+1)+1,&nvar+3]=pl_save;
  outest[(iby-1)*(&nvar+1)+1,&nvar+4]=ll_save;
  outest[(iby-1)*(&nvar+1)+1,&nvar+5]=it_save;
  outest[(iby-1)*(&nvar+1)+1,&nvar+6]=y[+];
  outest[(iby-1)*(&nvar+1)+1,&nvar+7]=n-y[+];

  do k=1 to &nvar;
   outest[(iby-1)*(&nvar+1)+1+k,&nvar+2]=iby;
   outest[(iby-1)*(&nvar+1)+1+k,&nvar+3]=pl0[k];
   outest[(iby-1)*(&nvar+1)+1+k,&nvar+4]=ll0[k];
   outest[(iby-1)*(&nvar+1)+1+k,&nvar+5]=it_save;
   outest[(iby-1)*(&nvar+1)+1+k,&nvar+6]=y[+];
   outest[(iby-1)*(&nvar+1)+1+k,&nvar+7]=n-y[+];
   outest[(iby-1)*(&nvar+1)+1,k]=b_save[k];
   outest[(iby-1)*(&nvar+1)+1+k,&nvar+1]=k-1;
   do kk=1 to &nvar;
    outest[(iby-1)*(&nvar+1)+1+kk,k]=Cov_save[k,kk];
   end;
  end;
 %end;

end;

create &global from global [colname={"_rby_" "Testtype" "ChiSq" "df" "P_value"}];
append from global;
close &global;

create &outtab from outtab [colname={"_rby_" "_VAR_" "BETA" "STDERR" "CI_LO" "CI_UP" "P_VALUE" "_ITER_" "_PENLIK_" "_LNLIKE_"}];
append from outtab;
close &outtab;

%if &nprof ne 0 %then %do;
 create &outprof from outprof[colname={"_rby_" "_var_" "_b_" "_profli_" "_normal_" "_refer_" "_bpl_" "_cipl_" "_bwald_" "_ciwald_"}];
 append from outprof;
 close &outprof;
%end;

%if &outest ne ' ' %then %do;
 create &outest from outest [colname={
  %if &noint=0 %then %do;
   "INTERCEP"
  %end;
  "&var1"
  %if &nvar>2 %then %do;
   %do k=2 %to &nvar-1;
     "&&var&k"
   %end;
  %end;
   "_code_" "_rby_" "_penlik_" "_lnlike_" "_IT_" "_RESP_" "_NORESP_"}];
 append from outest;
 close &outest;
%end;

%if &out~= %then %do;
 create __h from hvec_all[colname={"&h"}];
 append from hvec_all;
 close __h;
%end;

%if &test~= %then %do;
 create &outlr from p_spec[colname={"_CHI_" "_DF_" "p_value"}];
 append from p_spec;
 close &outlr;
%end;

quit;

%if &test~= %then %do;
 data &outlr;
 set &outlr;
 _rby_=_n_;
 run;

 data &outlr;
 merge &outlr _fby;
 by _rby_;
 _name_="%upcase(&test)";
 label _name_="Tested variable(s)" _chi_="Penalized log likelihood Chi-square"
  _df_="Degrees of freedom" p_value="Pr > Chi-square";
 run;
%end;

%if &nprof ne 0 %then %do;
 proc sort data=&outprof out=&outprof;
 by _rby_ _var_ _b_;
 run;

 data &outprof;
 merge &outprof _fby;
 by _rby_;
 length _name_ $ 8;
 %do j=0 %to &nvar-1;
  if _var_=&j then do;
   _name_="&&var&j";
  end;
 %end;
 if _ciwald_=0 and _bwald_=0 and _bpl_=0 and _cipl_=0 then do;
  _ciwald_=.;
  _bwald_=.;
  _bpl_=.;
  _cipl_=.;
 end;
 %let clevp=%sysevalf(100*(1-&alpha));
 label _name_="Variable" _b_="beta" _profli_="Profile penalized likelihood" _normal_="Wald" _refer_="&clevp.% reference line" _bpl_="beta" _bwald_="beta"
  _cipl_="&clevp.% Profile penalized likelihood c.i." _ciwald_="&clevp.% Wald c.i.";
 run;

 proc sort data=&outprof;
 by _rby_ _name_;
 run;
%end;

data &global;
merge &global _fby;
by _rby_;
if testtype=1 then      Type="Likelihood Ratio";
else if testtype=2 then Type="Wald            ";
label ChiSq="ChiSq" df="df" p_value="Pr>ChiSq" type="Type";
run;

data &outtab;
merge &outtab _fby;
by _rby_;
run;

%if &outest ne ' ' %then %do;
 data &outest;
 merge &outest _fby;
 by _rby_;
 length _type_ _name_ _link_ $ 8;
 _total_=_resp_+_noresp_;
 _LINK_="LOGIT";
 if _code_=-1 then do;
  _type_="PARMS";
  _name_="%upcase(&y)";
 end;
 if _code_=0 then do;
  _type_="COV";
  _name_="INTERCPT";
 end;
 %do co=1 %to &nvar-1;
  if _code_=&co then do;
   _type_="COV";
   _name_="%upcase(&&var&co)";
  end;
 %end;
 label _it_="Iterations" _penlik_="Penalized log likelihood"
 _lnlike_="Log likelihood" _resp_="Number of responses" _noresp_="Number of nonresponses" 
 _total_="Number of observations";
 run;
%end;


%if &out~= %then %do;
 data __oe;
 set &outest;
 if _code_=-1;
 run;


 proc iml;

 start getcov;
   do j=1 to &nvar;
    do jj=1 to &nvar;
     Cov[j,jj]=vcv[(&nvar+1)*(rby[i]-1)+1+j,jj];
    end;
   end;
 finish getcov;

 use _work;
 read all var("&&var&noint"
  %if &nvar-1 > 0 %then %do;
   %do j=&noint+1 %to &nvar-1+&noint;
    || "&&var&j"
   %end;
  %end;
  ) into x;
 read all var("_rby_") into rby;
 close _work;
 use __oe;
 read all var("&&var&noint"
  %if &nvar-1 > 0 %then %do;
   %do j=%eval(&noint+1) %to &nvar-1+&noint;
    || "&&var&j"
   %end;
  %end;
  ) into beta;
 close __oe;
 maxrby=nrow(beta);
 use &outest;
 read all var ("&&var&noint"
  %if &nvar-1 > 0 %then %do;
   %do j=%eval(&noint+1) %to &nvar-1+&noint;
    || "&&var&j"
   %end;
  %end;
   ) into vcv;
 close &outest;
 n_all=nrow(x);
 cov=repeat(0,&nvar,&nvar);

 pred=repeat(0,n_all,3);
 zalpha2=probit((&conflev+100)/200);
 do i=1 to n_all;
  eta_i=x[i,]*t(beta[rby[i],]);
  pred[i,1]=1/(1+exp(-eta_i));
  if (i=1) then run getcov;
  if (i>1) then do;
   if rby[i]>rby[i-1] then run getcov;
  end;
  sd_eta=sqrt(x[i,]*Cov*t(x[i,]));
  rbyi=rby[i];
*  print rbyi sd_eta Cov zalpha2;
  pred[i,2]=1/(1+exp(-eta_i+zalpha2*sd_eta));
  pred[i,3]=1/(1+exp(-eta_i-zalpha2*sd_eta));
 end;

 create &out from pred[colname={"&pred" "&lower" "&upper"}];
 append from pred;
 close &out;

 quit;

 data &out;
 merge _work __h &out;
 label &pred="Predicted probabilty" &lower="Lower &conflev.% c.l." &upper="Upper &conflev.% c.l."
  &h="Hat matrix diagonal";
 run;

%end;

data &outtab;
set &outtab;
length _name_ $ 8;
%do k=1 %to &nvar;
 %let km1=%eval(&k-1);
 %if &noint=0 %then %do;
  if _var_=&k then do;
   _name_="%upcase(&&var&km1)";
  end;
 %end;
 %else %do;
  if _var_=&k then do;
   _name_="%upcase(&&var&k)";
  end;
 %end;
%end;
%if &odds=1 %then %do;
 odds=exp(beta);
 or_lo=exp(ci_lo);
 or_up=exp(ci_up);
%end;
%if &noint=0 %then %do;
 _var_=_var_-1;
 %end;
label _var_="Variable No." _name_="Variable" beta="Parameter estimate" stderr="Standard Error" p_value="Pr > Chi-Square"
 ci_lo="Lower &conflev.% c.l." ci_up="Upper &conflev.% c.l." _iter_="Iterations" _penlik_="Penalized log likelihood"
 _lnlike_="Log likelihood"
 %do k=1 %to &nby;
  &&by&k="%upcase(&&by&k)"
 %end;
 %if &odds=1 %then %do;
  odds="Odds ratio" or_lo="Lower &conflev.% c.l." or_up="Upper &conflev.% c.l."
 %end;
;
run;

%if &print=1 %then %do;

 data ___a1;
 set &outtab;
 file print;
 if _n_=1 then do;
  put "   ";
*  put "*******************************************************************************";
  put "  ";
  put "  FFFFF L                Logistic regression";
  put "  F     L                with Firth`s bias reduction:";
  put "  FFF   L       ";
  put "  F     L                A solution to the problem of separation";
  put "  F     LLLLL            in logistic regression";
  put "  ";
  put "   ";
*  put "*******************************************************************************";
  put "  ";
  put "  Author:                Georg Heinze";
  put "  Version:               &version";
  put "   ";
  put "  Methods published in:  Heinze, G. & Schemper, M. (2002). A solution";
  put "                         to the problem of separation in logistic";
  put "                         regression. Statistics in Medicine 21(16)";
  put "                         2409-2419.";
  put "   ";
  put "  Data set:              %upcase(&data)";
  put "  Dependent variable:    %upcase(&y)";
  %if &nvar-1+&noint > 1 %then %do;
   put "  Independent variables: %upcase(&varlist)";
  %end;
  %else %do;
   put "  Independent variable:  %upcase(&varlist)";
  %end;
  put " ";
  put "  Table with parameter estimates saved as %upcase(&outtab).";
  %if &outest ne ' ' %then %do;
   put "  Estimates and covariance matrix saved as %upcase(&outest).";
  %end;
  %if &out~= %then %do;
   put "  Predicted probabilities, confidence limits";
   put "  and hat matrix diagonals saved as %upcase(&out).";
  %end;
 end;
 run;
 %do rby=1 %to &maxby;
  proc print data=&Outest noobs label;
  title4 "Model fitting information";
  var _it_ _penlik_ _resp_ _noresp_ _total_;
  where  _NAME_="%upcase(&y)"
  %if &byproc=1 %then %do;
   and _rby_=&rby;
   by &by;
  %end;
  %else %do; ; %end;
  %if &standard=1 %then %do;
   title5 "NOTE: Penalization computed from standardized data.";
  %end;
  run;

  proc print data=&global noobs;
  title4 "Testing global null hypothesis beta=0";
  var Type ChiSq df p_value;
  format p_value pvalue6.;
  %if &byproc=1 %then %do;
   where _rby_=&rby;
   by &by;
  %end;
  run;

  proc print data=&outtab noobs label;
  %if &pl=0 %then %do;
   title4 "FL estimates and Wald confidence limits and tests";
  %end;
  %else %do;
   title4 "FL estimates, profile penalized likelihood confidence limits";
   title5 "and penalized likelihood ratio tests";
  %end; 
  %if &plint.&noint=00 %then %do;
   title6 "NOTE: Confidence interval for Intercept based on Wald method.";
  %end;
  var _name_ beta stderr ci_lo ci_up p_value;
  %if &byproc=1 %then %do;
   by &by;
   where _rby_=&rby;
  %end;
  format p_value pvalue6.;
  run;
  title4;
  %if &odds=1 %then %do;
   proc print data=&outtab noobs label;
   %if &pl=0 %then %do;
    title4 "FL odds ratio estimates and Wald confidence limits and tests";
   %end;
   %else %do;
    title4 "FL odds ratio estimates, profile penalized likelihood confidence limits";
    title5 "and penalized likelihood ratio tests";
   %end; 
    var _name_ odds or_lo or_up p_value;
	where _var_ ne 0;
	%if &byproc=1 %then %do;
     by &by;
     where _rby_=&rby;
    %end;
	format p_value pvalue6.;
   run;
   title4;
  %end; 

  %if &test~= %then %do;
   proc print data=&outlr noobs label;
   var _name_ _chi_ _df_ p_value;
   %if &byproc=1 %then %do; 
    by &by; 
    where _rby_=&rby;
   %end; 
   format p_value pvalue6.;
   run;
  %end;

  %if &nprof ne 0 %then %do;
   symbol1 i=join v=none c=black line=1;
   symbol2 i=join v=none c=black line=2;
   symbol3 i=join v=plus c=black line=1;
   symbol4 i=join v=plus c=black line=3;
   legend1 down=4 label=none
        position=(bottom center outside)
        mode=reserve
        shape=symbol(10,1) value=("profile penalized likelihood" "&clevp.% reference line" 
         "point estimate and &clevp.% profile penalized likelihood confidence interval" 
         "point estimate and &clevp.% Wald confidence interval");
   axis1 label=(angle=90);
   proc gplot data=&outprof;
   where _rby_=&rby;
   plot (_profli_ _refer_)*_b_  _cipl_*_bpl_ _ciwald_*_bwald_/ overlay legend=legend1 vaxis=axis1;
   by %if &byproc=1 %then %do; &by %end; _name_;
   run;
  %end;


 %end;

%end;


%mend;
