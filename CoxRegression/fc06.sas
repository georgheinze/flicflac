%macro fc06(data=_last_, time1=0, time2=t, time=, cens=ic, censval=0, varlist=, by=, print=1, offset=,
          outmod=_mod, outest=_est, outtab=_tab, outtest=_test, outprof=_prof, outs0=, s0method=EMP,
          maxit=50, maxhs=5, maxstep=2.5, epsilon=0.000001,
          alpha=0.05, risk=0, out=, test=, pl=1, profile=, profsel=, profser=, profn=100, firth=1, notes=0,
          testtype=SCORES, ft=, ftmap=, tdenames=,
          global=_global,
          call=dll, path=%str(S:\sasdll\));

%let version=2006.01 (beta);
%let build=11;

%let mynameis=fc06;

%macro picture;
  put "  FFFFF  CCC  000   666       Cox regression";
  put "  F     C    0   0 6          using &penalize maximum likelihood estimation";
  put "  FFFF  C    0   0 6666  ";
  put "  F     C    0   0 6   6 ";
  put "  F      CCC  000   666  ";
%mend;

* FC06                                                                                                    ;

* estimates Cox model with Firths modified score function                                               ;
* Start-Stop-Syntax allowed ;
* time-dependent effects allowed;
* Georg Heinze April 2004 - January 2006                                                                ;



* uses the file fc06.dll that must be installed in a directory accessible by SAS                       ;

* THE USER TAKES RESPONSABILITY FOR CORRECT ENTRIES                                                     ;
* THERE ARE NO CHECKS FOR UNPLAUSIBLE ENTRIES (MISSING DATASETS, ...)                                   ;
* INCORRECT OPTION SETTINGS MAY LEAD TO CRASHES OF SAS                                                  ;

*************************************************************************                               ;
*              SO SAVE YOUR WORK BEFORE CALLING FC06!!!                 *                               ;
*************************************************************************                               ;

********************************************************************************************************;
* KNOWN PROBLEMS:                                                                                       ;
********************************************************************************************************;
* SAS crashes when calling FC06:                                                                        ;

*         Are all options correct (no missing data sets, no missing variables in these data sets, ...)? ;
*            No -> correct the entries                                                                  ;
*            Yes ->
*                FIRTH=1 -> try a lower value for maxstep (e.g. 0.3)                                    ;
*                FIRTH=0 -> there might be monotone likelihood and an overflow in the FORTRAN Routine.  ;
*                           try lower values of maxstep and maxit (0.5 and 8, respectively)             ;
*                                                                                                       ;
* Please report all other problems to Georg.Heinze@meduniwien.ac.at                                       ;
********************************************************************************************************;


%let call=%upcase(&call);

%if &call=EXE and (&pl=3 or &pl=2) %then %do;
 %put ERROR: call=EXE not available with pl>1;
 %goto ende;
%end;


%if &call=DLL %then %do;
 filename SASCBTBL "&path.&mynameis..def";
%end;
%if &call=EXE %then %do;
 data _ttt;
 file "copyfc.bat";
 put "copy &path.&mynameis.exe.exe &mynameis.exe.exe";
 put "copy &path.&mynameis..dll &mynameis..dll";
 put "exit";
 run;
 data _ttt;
 file "&mynameis..bat";
 put "&mynameis.exe";
 put "exit";
 run;
 data _ttt;
 file "delfcfiles.bat";
 put "copy ioarray.txt ioarrayoutput.txt";
 put "copy parms.txt parmsoutput.txt";
 put "copy cards.txt cardsoutput.txt";
 put "copy parmsinput.txt parmsinputoutput.txt";
 put "del ioarray.txt";
 put "del parms.txt";
 put "del cards.txt";
 put "del parmsinput.txt";
 put "exit";
 run;
 x copyfc.bat;
%end;


%macro toputc(a=, b=, f=);
 toput=cards[&a,&b];
 put toput %if &f ne %then %do; &f. %end; " " @;
%mend;

%macro toinputc(a=, b=);
 input toinput @;
 cards[&a,&b] = toinput;
%mend;

%macro toputl(a=, b=);
 toput=iopl[&a,&b];
 put toput " " @;
%mend;

%macro toinputl(a=, b=);
 input toinput @;
 iopl[&a,&b] = toinput;
%mend;

%macro toputi(a=, b=);
 toput=ioarray[&a,&b];
 put toput " " @;
%mend;

%macro toinputi(a=, b=);
 input toinput @;
 ioarray[&a,&b] = toinput;
%mend;

%macro toputp(a);
 toput=parmspl[&a];
 put toput " " @;
%mend;

%macro toinputp(a);
 input toinput @;
 parmspl[&a] = toinput;
%mend;

%macro toputh(a);
 toput=hugo[&a];
 put toput @;
%mend;

%macro toinputh(a);
 input toinput @;
 hugo[&a] = toinput;
%mend;



%macro thecall(what=fc);

%if &call=DLL %then %do;
 %if &what=fc %then %do;
   CALL modulei('*E','FIRTHCOX',cards,hugo,ioarray);
 %end;
 %else %do;
   CALL modulei('*E','PLCOMP',cards,parmspl,iopl);
 %end;
%end;
%else %if &call=EXE %then %do;
  file "cards.txt";
  do icw1=1 to nrow(cards);
   do icw2=1 to ncol(cards);
    if icw2=ncol(cards)-2 then do;
	 %toputc(a=icw1,b=icw2,f=20.17);
    end;
	else do;
     %toputc(a=icw1,b=icw2);
	end;
   end;
   put ;
  end;
  put " ";
  closefile "cards.txt";
  file "ioarray.txt";
  %if &what=fc %then %do;
   do icw1=1 to nrow(ioarray);
    do icw2=1 to ncol(ioarray);
	 %toputi(a=icw1,b=icw2);
	end;
    put;
   end;
  %end;
  %else %do;
   do icw1=1 to nrow(iopl);
    do icw2=1 to ncol(iopl);
	 %toputl(a=icw1,b=icw2);
	end;
    put;
   end;
  %end;
  put " ";
  closefile "ioarray.txt";
  file "parms.txt";
*  print hugo;
  %if &what=fc %then %do;
   put 0;
  %end;
  %else %do;
   put 1;
  %end;
  nrc=nrow(cards);
  put nrc;
  ncc=ncol(cards);
  put ncc;
  %if &what=fc %then %do;
   nri=ncol(ioarray);
   put nri;
   do icw1=1 to 14;
    %toputh(icw1);
	put;
   end;
  %end;
  %else %do;
   nri=ncol(iopl);
   put nri;
   do icw1=1 to 14;
    %toputp(icw1);
	put;
   end;
  %end;
  put " ";
  closefile "parms.txt";
  file "parmsinput.txt";
  %if &what=fc %then %do;
   put 0;
  %end;
  %else %do;
   put 1;
  %end;
  nrc=nrow(cards);
  put nrc;
  ncc=ncol(cards);
  put ncc;
  %if &what=fc %then %do;
   nri=ncol(ioarray);
   put nri;
   do icw1=1 to 14;
    %toputh(icw1);
	put;
   end;
  %end;
  %else %do;
   nri=ncol(iopl);
   put nri;
   do icw1=1 to 14;
    %toputp(icw1);
	put;
   end;
  %end;
  put " ";
  closefile "parmsinput.txt";
*  store _ALL_ module=_all_;
*  quit;
  file log;
  put "CALL: &mynameis.EXE.EXE";
  x "&mynameis..bat";
*  proc iml;
*  load _ALL_ module=_ALL_;
  file log;
  infile "parms.txt";
  %if &what=fc %then %do;
*   input junk;
   do icw1=1 to 14;
    %toinputh(icw1);
	input;
   end;
*   print hugo;
  %end;
  %else %do;
*   input junk;
   do icw1=1 to 14;
    %toinputp(icw1);
	input;
   end;
*   print parmspl;
  %end;
  closefile "parms.txt";
  infile "ioarray.txt" flowover;
  %if &what=fc %then %do;
   do icw1=1 to nrow(ioarray);
    do icw2=1 to ncol(ioarray);
	 %toinputi(a=icw1,b=icw2);
	end;
    input;
   end;
  %end;
  %else %do;
   do icw1=1 to nrow(iopl);
    do icw2=1 to ncol(iopl);
	 %toinputl(a=icw1,b=icw2);
	end;
    input;
   end;
  %end;
  closefile "ioarray.txt";
  x "delfcfiles.bat";
%end;

%mend;



%if &by~= %then %let byproc=1;
%else %let byproc=0;

%let nby=0;
%if &by ne ' ' %then %do;
 %do %while(%scan(&by,&nby+1)~=);
  %let nby=%eval(&nby+1);
  %let by&nby=%scan(&by,&nby);
 %end;
%end;

%let nvar=0;
%do %while(%scan(&varlist,&nvar+1)~=);
 %let nvar=%eval(&nvar+1);
 %let var&nvar=%scan(&varlist,&nvar);
 %let vg&nvar=&&var&nvar;
%end;

%let ntde=0;
%do %while(%scan(&ft,&ntde+1)~=);
 %let ntde=%eval(&ntde+1);
 %let ftmap&ntde=%scan(&ftmap,&ntde);
 %let tdename&ntde=%scan(&tdenames,&ntde);
 %let numm=%eval(&nvar+&ntde);
 %let vg&numm=&&tdename&ntde;
%end;

%let ntde=0;
%let stopit=0;
%do %while(&stopit=0);
 %let plusone=%eval(&ntde+1);
 %let ft&plusone=%scan(&ft,&plusone," ");
 %if "&&ft&plusone"="" %then %do;
  %let stopit=1;
 %end;
 %else %do;
  %let ntde=%eval(&ntde+1);
 %end;
%end;

%let genvar=;
%let gentype=;


%let ngv=0;
%do %while(%scan(&genvar,&ngv+1)~=);
 %let ngv=%eval(&ngv+1);
 %let gv&ngv=%scan(&genvar,&ngv);
%end;


%do j=1 %to &nvar;
 %let weit&j=0;
 %do jj=1 %to &ngv;
  %if %upcase(&&var&j)=%upcase(&&gv&jj) %then %let weit&j=1;
 %end;
%end;

%let ntest=0;
%do %while(%scan(&test,&ntest+1)~=);
 %let ntest=%eval(&ntest+1);
 %let test&ntest=%scan(&test,&ntest);
%end;

%let testtype=%upcase(&testtype);
%if %substr(&testtype,1,1)=S %then %let testtype=SCORES;
%if %substr(&testtype,1,1)=L %then %let testtype=LR;

%let nprof=0;
%do %while(%scan(&profile,&nprof+1)~=);
 %let nprof=%eval(&nprof+1);
 %let prof&nprof=%scan(&profile,&nprof);
%end;

%if &profn <3 %then %let profn=3;

%let conflev=%sysevalf(1-&alpha);


* names for output;

%if &firth=1 %then %do;
 %let estname=FC;
 %let penalize=penalized;
 %let penalcap=Penalized log;
 %let pentname=Penalized likelihood ratio;
%end;
%else %do;
 %let estname=ML;
 %let penalize=;
 %let penalcap=Log;
 %let pentname=Likelihood ratio;
%end;

%if &testtype=LR %then %do;
  %let penlik=&penalcap likelihood;
  %let penli0=Restricted &penalize log likelihood; 
  %let restrict=_penlik_ _penli0_;
%end;
%else %if &testtype=SCORES %then %do;
  %let penlik=Scores statistic;
  %let penli0=Restricted scores statistic;
  %let pentname=Partial scores;
  %let restrict=;
%end;

%if &pl=0 %then %do;
 %let clname=Wald;
 %let testname=Wald;
%end;
%if &pl=1 %then %do;
 %let clname=profile &penalize likelihood;
 %let claname=Profile &penalize likelihood;
 %let statval=profile &penalize likelihood;
 %let testname=&penalize likelihood ratio;
%end;
%if &pl=2 %then %do;
 %let clname=profile &penalize likelihood;
 %let claname=Profile &penalize likelihood;
 %let statval=profile &penalize likelihood;
 %let testname=&penalize likelihood ratio;
%end;
%if &pl=3 %then %do;
 %let clname=partial scores;
 %let claname=Partial scores;
 %let statval=partial scores statistic;
 %let testname=partial scores;
%end;


data _work;
set &data;
%if &time ne %then %do;
 %let time2=&time;
 %let time1=_start_;
 _start_=0;
%end;

*&time=&time+1;
if &time2>&time1;

if &cens ne . then do;
 if &cens=&censval then _ic=0;
 else _ic=1;
end;

if &cens ne .;

%do k=1 %to &nvar;
 if &&var&k ne .;
%end;

%if &by~= %then %do;
 %let byproc=1;
%end;
%else %do;
 %let byproc=0;
 %let by=_BY_;
 %let by0=_BY_;
 _By_=1;
%end;
run;


proc sort;
by &by &time2 descending _ic;
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
run;

data _work;
merge _work _fby;
by &by;
run;

%if &offset~= %then %do;
 %if &by~= %then %do;
  proc sort data=&offset out=&offset;
  by &by;
  run;

  data &offset;
  merge &offset _fby;
  by &by;
  run;
 %end;
%end;

* rby is now column with numbers from 1 to #by-groups ;


proc means noprint;
var _rby_;
output out=_maxby max=maxby;
run;

data _maxby;
set _maxby;
call symput('maxby',maxby);
run;

%let maxby=%cmpres(&maxby);

data _offset_;
%if &offset~= %then %do;
 set &offset;
 %if &ntest ne 0 %then %do;
  %do kt=1 %to &ntest;
    %do kt2=1 %to &nvar;
     %if %upcase(&&test&kt)=%upcase(&&var&kt2) %then %do;
      if _rby_<=&maxby/2 then _off&kt2._=0;
     %end;
    %end;
    %do kt2=1 %to &ntde;
	 %let numm=%eval(&nvar+&kt2);
     %if %upcase(&&test&kt)=%upcase(&&tdename&kt2) %then %do;
      if _rby_<=&maxby/2 then _off&numm._=0;
     %end;
    %end;
  %end;
 %end;
%end;
%else %do;
 do _rby_=1 to &maxby;
  %do j=1 %to &nvar+&ntde;
   _off&j._=.;
  %end;
  output;
 end;
%end;
run;

data _offset_;
set _offset_;
%do k=1 %to &nvar+&ntde;
 if _off&k._=. then do;
  _off&k._=0;
  _fl&k._=1;
 end;
 else do;
  _fl&k._=0;
 end;
%end;
run;


*proc rank data=_work out=_workr;
*var &time;
*ranks trank;
*by _rby_;
*run;


proc means noprint data=_work;
var &cens;
output out=_n n=_n sum=_sum;
by _rby_;
run;


data _n;
merge _n _offset_;
by _rby_;
run;


proc iml;
use _n;
read all var("_n") into nby;
close _n;

use _work;
read all var("&var1" %if &nvar>1 %then %do; %do k=2 %to &nvar; || "&&var&k" %end; %end;) into xall;
read all var("&time1" || "&time2") into timeall;
read all var("_ic") into censall;
read all var("_rby_") into bygroup;
close _work;

use _offset_;
read all var("_off1_" %if &nvar+&ntde>1 %then %do; %do j=2 %to &nvar+&ntde; || "_off&j._" %end; %end;) into offs;
read all var("_fl1_" %if &nvar+&ntde>1 %then %do; %do j=2 %to &nvar+&ntde; || "_fl&j._" %end; %end;) into flags;
close _offset_;

k=&nvar+&ntde;
maxrby=max(bygroup);
%if &ntde ne 0 %then %do;
 ftmap=t({&ftmap});
%end;

outtab=repeat(0,(k)*maxrby,8);
outest=repeat(0,(k+1)*maxrby,3+k);
outmod=repeat(0,maxrby,12);
%if &ntest ne 0 %then %do;
 outtest=repeat(0,maxrby,6);
%end;
%if &nprof ne 0 %then %do;
 outprof=repeat(0,maxrby*&nprof*&profn,10);
%end;
index=0;
%do rby=1 %to &maxby;
 rby=&rby;
 %if &notes=1 %then %do;
  %put NOTE: &mynameis. is processing data set &rby of &maxby..;
 %end;
 *cards=repeat(0,nby[rby],&nvar+3+&nvar);
* n=nby[rby];
 cards=xall[loc(bygroup=rby),]||timeall[loc(bygroup=rby),]||censall[loc(bygroup=rby)];
 n=nrow(cards);
 %if &ngv>0 %then %do;
  %if %upcase(&gentype)=N %then %do;
   sp=(n+1)-t(1:n);
   sp=sp#cards[,&nvar+3];
  %end;
  %else %if %upcase(&gentype)=SQRTN %then %do;
   sp=sqrt((n+1)-t(1:n));
   sp=sp#cards[,&nvar+3];
  %end;
  %else %if %upcase(&gentype)=1 %then %do;
   sp=repeat(1,n,1);
   sp=sp#cards[,&nvar+3];
  %end;
  * normalize such that sum of weights equals #risksets;
  sp=sp/sp[+]*cards[+,&nvar+3];
*  print "Used data and normalized weights", rby,cards sp;
 %end;
 %else %do;
  sp=repeat(1,n,1);
 %end;
 %if &ntde ne 0 %then %do;
  t2=timeall[loc(bygroup=rby),2];
  ft=repeat(0,n,&ntde);
  do i=1 to n;
   _time_=t2[i];
   %do j=1 %to &ntde;
    ft[i,&j]= &&ft&j ;
   %end;
  end;
 %end;
 xmean=cards[:,1:&nvar];
 cc=cards[,1:&nvar]#cards[,1:&nvar];
 xstd=sqrt( cc[:,] - xmean##2 );
 %if &ntde > 0 %then %do;
  xstd=xstd||t(xstd[ftmap]);
 %end;
* xstd=repeat(1,&nvar,1);
 do j=1 to &nvar;
*  cards[,j]=(cards[,j]-xmean[j])/xstd[j];
   cards[,j]=(cards[,j])/xstd[j];
 end;

 index=index+n;
 ioarray=repeat(0,3+k,k);
 hugo=repeat(0,14,1);
 cards=cards %do j=1 %to &nvar; 
              %if &&weit&j=1 %then %do; || sp %end;
			  %else %do; || repeat(1,n,1) %end;
			 %end;
             %do j=1 %to &ntde;
			  %let numm=%eval(&&ftmap&j);
			  %if &&weit&numm=1 %then %do; || sp %end;
			  %else %do; || repeat(1,n,1) %end;
             %end;;

 %if &ntde ne 0 %then %do;
  cards=cards||ft;
 %end;
 hugo=n//&nvar//&firth//&maxit//&maxhs//&maxstep//&epsilon//0//0//0//0//0//&ngv//&ntde;
* do j=1 to k;
*  ioarray[1,j]=flags[rby,j];
*  ioarray[2,j]=offs[rby,j];
*  if j<=&nvar then do;
*  ioarray[2,j]=offs[rby,j]*xstd[j];
*  end;
* end;
 ioarray[1,]=flags[rby,];
 ioarray[2,]=offs[rby,]#xstd;
 %if &ntde ne 0 %then %do;
  ioarray[4,&nvar+1:k]={&ftmap};
 %end;

*print rby, cards, ioarray, hugo;
 file log;
* put "bis hierher";
* print cards, hugo, ioarray;
* show space;
 %thecall(what=fc);
* print hugo, ioarray;
 scorechi=hugo[7];
 hugo[7]=&epsilon;
 if hugo[8]>0 then put "ERROR: Problem in FIRTHCOX";
 iosave=ioarray;
 penlike=hugo[11];
 penlike0=hugo[12];
 if k=1 & flags[rby,1]=0 then
  penlike0=penlike;
 
 iter=hugo[10];
 isep=hugo[9];
 jcode=hugo[8];

 b=repeat(0,k,1);
 vm=repeat(0,k,k);
 se=repeat(0,k,1);

 vm=ioarray[4:(k+3),];
 do j=1 to k;
  b[j]=ioarray[3,j]/xstd[j];
  do jj=1 to k;
   vm[j,jj]=ioarray[3+j,jj]/xstd[j]/xstd[jj];
  end;
  se[j]=sqrt(vm[j,j]);
 end;
* do j=&nvar+1 to k;
*  b[j]=ioarray[3,j]/xstd[ftmap[j-&nvar]];
*  do jj=1 to &nvar;
*   vm[j,jj]=vm[j,jj]/xstd[jj]/xstd[ftmap[j-&nvar]];
*   vm[jj,j]=vm[j,jj];
*  end;
*  do jj=&nvar+1 to k;
*   vm[j,jj]=vm[j,jj]/xstd[ftmap[jj-&nvar]]/xstd[ftmap[j-&nvar]];
*   vm[jj,j]=vm[j,jj];
*  end;
*  se[j]=sqrt(vm[j,j]);
* end;
 outtab[(rby-1)*k+(1:k),8]=t(ioarray[3,(1:k)]);
 modchi=(penlike-penlike0)*2;
 pmodchi=-1;
 if modchi>0 then
  pmodchi=1-probchi(modchi,k);

 %if &pl=1 %then %do;
  iopl=repeat(0,8,k);
  do jj=1 to k;
   iopl[1,jj]=1;
   iopl[3,jj]=b[jj]*xstd[jj];
  end;
*  do jj=&nvar+1 to k;
*   iopl[1,jj]=1;
*   iopl[3,jj]=b[jj]*xstd[ftmap[jj-&nvar]];
*  end;
  parmspl=repeat(0,14,1);
  do j=1 to 6;
   parmspl[j]=hugo[j];
  end;
  parmspl[7]=&epsilon;
  parmspl[8]=&alpha;
  parmspl[9]=0;
  parmspl[13]=&ngv; 
  parmspl[14]=&ntde;
  %if &ntde ne 0 %then %do;
   iopl[4,&nvar+1:k]={&ftmap};
  %end;
  %if &ntest=0 %then %do;
   %thecall(what=pl); 
  %end;
  %else %do;
   if rby>maxrby/2 then do;
    parmspl[7]=&epsilon;
    %thecall(what=pl); 
   end;
  %end;
*  print iopl;
  ci=repeat(0,k,2);
  pvalue=repeat(0,k,1);
  do j=1 to k;
   do l=1 to 2;
    ci[j,l]=iopl[3+l,j]/xstd[j];
   end;
   pvalue[j]=iopl[6,j];
  end;
 %end;
 %else %if &pl=2 %then %do;
  * emergency estimation of pl c. l. by binary search;
  ci=repeat(0,k,2);
  pvalue=repeat(0,k,1);

  if &ntest=0 | rby>maxrby/2 then do;
   do j=1 to k;
    do l=1 to 2;
	 estimate=b[j];
	 li1=b[j];
	 * first calibrate "outer" limit;
	 target=penlike-0.5*cinv(&conflev,1);
	 hugo[11]=target+1;
	 it_cal=1;
	 do while(hugo[11]>target & it_cal<&maxit);
	   it_cal=it_cal+0.5;
       li2=b[j]-it_cal*se[j]*(-1)**(l=2);
	   llike1=hugo[11];
	   hugo[(1:14)]=n//&nvar//&firth//&maxit//&maxhs//&maxstep//&epsilon//0//0//0//0//0//&ngv//&ntde;
       do jj=1 to k;
        ioarray[1,jj]=flags[rby,jj];
        ioarray[2,jj]=offs[rby,jj]*xstd[jj];
       end;
       ioarray[1,j]=0;
	   ioarray[2,j]=li2*xstd[j];
       %if &ntde ne 0 %then %do;
        ioarray[4,&nvar+1:k]={&ftmap};
	   %end;
*	   print j l it_cal cards hugo ioarray;
       %thecall(what=fc);
*	   print j l it_cal "done.";
	   llike2=hugo[11];
       s_chi=hugo[7];
       hugo[7]=&epsilon;
     end;
	 li1=b[j]-(it_cal-1)*se[j]*(-1)**(l=2);

	 it_bin=0;
	 li1=li1*xstd[j];
	 li2=li2*xstd[j];
	 do while(abs(li1-li2)>&epsilon & it_bin<&maxit);
*	   print j l it_bin li1 li2;
	   it_bin=it_bin+1;
	   hugo[(1:14)]=n//&nvar//&firth//&maxit//&maxhs//&maxstep//&epsilon//0//0//0//0//0//&ngv//&ntde;
       do jj=1 to k;
        ioarray[1,jj]=flags[rby,jj];
        ioarray[2,jj]=offs[rby,jj]*xstd[jj];
       end;
       ioarray[1,j]=0;
*	   ioarray[2,j]=(li1+li2)/2;
	   ioarray[2,j]=(target-llike1)/(llike2-llike1)*(li2-li1)+li1;
	   %if &ntde ne 0 %then %do;
	    ioarray[4,&nvar+1:k]={&ftmap};
	   %end;
       %thecall(what=fc);
       s_chi=hugo[7];
       hugo[7]=&epsilon;
       
/**       if hugo[11]>target then do;
*        li1=ioarray[2,j];
*		llike1=hugo[11];
*	   end;
       else do;*/
	    li1=li2;
		llike1=llike2;
        li2=ioarray[2,j];
	    llike2=hugo[11];
		
*	   end;
	 end;
	 ci[j,l]=(li1+li2)/2/xstd[j];
	end;
	hugo[(7:14)]=&epsilon//0//0//0//0//0//&ngv//&ntde;
    do jj=1 to k;
     ioarray[1,jj]=flags[rby,jj];
     ioarray[2,jj]=offs[rby,jj]*xstd[jj];
    end;
	if k>1 then do;
     ioarray[1,j]=0; 
 	 ioarray[2,j]=0;
	 %if &ntde ne 0 %then %do;
	  ioarray[4,&nvar+1:k]={&ftmap};
	 %end;
     %thecall(what=fc);
	 pvalue[j]=1-probchi(2*(penlike-hugo[11]),1);
	end;
	else do;
	 pvalue[j]=1-probchi(2*(penlike-penlike0),1);
	end;
   end;
  end;
 %end;
 %else %if &pl=3 %then %do;
  * estimation of partial scores c. l. by binary search;
  ci=repeat(0,k,2);
  pvalue=repeat(0,k,1);

  if &ntest=0 | rby>maxrby/2 then do;
   do j=1 to k;
    * calculate target value of score statistic;
	hugo[(1:14)]=n//&nvar//1//&maxit//&maxhs//&maxstep//&epsilon//0//0//0//0//0//&ngv//&ntde;
    do jj=1 to k;
     ioarray[1,jj]=flags[rby,jj];
     ioarray[2,jj]=offs[rby,jj]*xstd[jj];
    end;
    ioarray[1,j]=0;
   	ioarray[2,j]=b[j]*xstd[j];
	%if &ntde ne 0 %then %do;
	 ioarray[4,&nvar+1:k]={&ftmap};
	%end;
    %thecall(what=fc);
    scoest=hugo[9];
	hugo[7]=&epsilon;
    do l=1 to 2;
	 estimate=b[j];
	 li1=b[j];
	 * first calibrate "outer" limit;
	 target=cinv(&conflev,1);
	 hugo[9]=target-1;
	 it_cal=0;
	 do while(hugo[9]<target & it_cal<&maxit);
	   it_cal=it_cal+1;
       li2=b[j]-it_cal*se[j]*(-1)**(l=2);
  	   hugo[(1:14)]=n//&nvar//1//&maxit//&maxhs//&maxstep//&epsilon//0//0//0//0//0//&ngv//&ntde;
       do jj=1 to k;
        ioarray[1,jj]=flags[rby,jj];
        ioarray[2,jj]=offs[rby,jj]*xstd[jj];
       end;
       ioarray[1,j]=0;
	   ioarray[2,j]=li2*xstd[j];
	   %if &ntde ne 0 %then %do;
	    ioarray[4,&nvar+1:k]={&ftmap}; 
	   %end;
       %thecall(what=fc);
       s_chi=hugo[9];
*	   print l it_cal s_chi;
*       hugo[7]=&epsilon;
     end;
	 li1=b[j]-(it_cal-1)*se[j]*(-1)**(l=2);
	 li2=b[j]-(it_cal+1)*se[j]*(-1)**(l=2);

	 it_bin=0;
	 li1=li1*xstd[j];
	 li2=li2*xstd[j];
	 do while(abs(li1-li2)>&epsilon & it_bin<&maxit);
	   *print j l it_bin li1 li2;
	   it_bin=it_bin+1;
	   hugo[(1:14)]=n//&nvar//&firth//&maxit//&maxhs//&maxstep//&epsilon//0//0//0//0//0//&ngv//&ntde;
       do jj=1 to k;
        ioarray[1,jj]=flags[rby,jj];
        ioarray[2,jj]=offs[rby,jj]*xstd[jj];
       end;
       ioarray[1,j]=0;
	   ioarray[2,j]=(li1+li2)/2;
	   %if &ntde ne 0 %then %do;
	    ioarray[4,&nvar+1:k]={&ftmap};
	   %end;
       %thecall(what=fc);
       s_chi=hugo[9];
 *      hugo[7]=&epsilon;
       
       if hugo[9]<target then li1=ioarray[2,j];
       else li2=ioarray[2,j];
	 end;
	 ci[j,l]=(li1+li2)/2/xstd[j];
	end;
	hugo[(7:14)]=&epsilon//0//0//0//0//0//&ngv//&ntde;
    do jj=1 to k;
     ioarray[1,jj]=flags[rby,jj];
     ioarray[2,jj]=offs[rby,jj]*xstd[jj];
    end;
*	if k>1 then do;
     ioarray[1,j]=0; 
 	 ioarray[2,j]=offs[rby,j]*xstd[j];
	 %if &ntde ne 0 %then %do;
	  ioarray[4,&nvar+1:k]={&ftmap}; 
	 %end;
	 hugo[7]=&epsilon;
*	 print j,cards,hugo,ioarray;
     %thecall(what=fc);
	 pvalue[j]=1-probchi(hugo[9],1);
*	end;
*	else do;
*	 pvalue[j]=1-probchi(2*(penlike-penlike0),1);
*	end;
   end;
  end;
 %end;
 %else %do;
  ci=repeat(0,k,2);
  pvalue=repeat(0,k,1);
  do j=1 to k;
   ci[j,1]=b[j]-probit(1-&alpha/2)*se[j];
   ci[j,2]=b[j]+probit(1-&alpha/2)*se[j];
   pvalue[j]=2*(1-probnorm(abs(b[j]/se[j])));
  end;
*  print b se ci pvalue;
 %end;

 wald=t(b)*inv(vm)*b;
 pwald=1-probchi(wald,k);

* print hugo modchi pmodchi wald pwald b se ci pvalue;
 outest[(rby-1)*(k+1)+1,1]=rby;
 outest[(rby-1)*(k+1)+1,k+2]=penlike;
 outest[(rby-1)*(k+1)+1,k+3]=0;
 do j=1 to k;
  outest[(rby-1)*(k+1)+1,j+1]=b[j];
  outest[(rby-1)*(k+1)+1+j,1]=rby;
  outest[(rby-1)*(k+1)+1+j,k+2]=penlike;
  outest[(rby-1)*(k+1)+1+j,k+3]=j;
  do jj=1 to k;
   outest[(rby-1)*(k+1)+1+j,jj+1]=vm[jj,j];
  end;
  outtab[(rby-1)*k+j,(1:7)]=rby||j||b[j]||se[j]||ci[j,1]||ci[j,2]||pvalue[j];
 end;
 pscore=-1;
 if scorechi>0 then pscore=(1-probchi(scorechi,k));
 outmod[rby,(1:12)]=rby||penlike||penlike0||modchi||pmodchi||
                    iter||cards[+,&nvar+3]||(n-cards[+,&nvar+3])||scorechi||pscore||wald||pwald;

 ioas=ioarray;
 %if &ntest ne 0 %then %do;
  ioarray[1,]=flags[rby,];
  %do j=1 %to &ntest;
   %do jj=1 %to &nvar+&ntde;
    %if %upcase(&&test&j)=%upcase(&&vg&jj) %then %do;
	 %if &testtype=LR %then %do;
	  ioarray[1,&jj]=1-flags[rby,&jj];
	 %end;
	 %else %do;
	  ioarray[1,&jj]=0;
	 %end;
   	 ioarray[2,&jj]=offs[&jj]*xstd[&jj];
    %end;
   %end;
  %end;
  hugo[7]=&epsilon;
  %if &ntde ne 0 %then %do;
   ioarray[4,&nvar+1:k]={&ftmap}; 
  %end;
*  print "TEST",cards,hugo,ioarray;
  %thecall(what=fc); 
  outtest[rby,1]=rby;
  %if &testtype=LR %then %do;
   %if &offset~= %then %do;
    outtest[rby,(2:3)]=penlike||hugo[11];
   %end;
   %else %do;
    outtest[rby,(2:3)]=hugo[11]||penlike;
   %end;
   outtest[rby,(4:6)]=(2*abs(penlike-hugo[11]))||&ntest||(1-probchi((2*abs(penlike-hugo[11])),&ntest));
  %end;
  %if &testtype=SCORES %then %do;
   outtest[rby,(2:6)]=0||0||hugo[9]||&ntest||(1-probchi(hugo[9],&ntest));
  %end;
  ioarray=ioas;
 %end;
 %if &nprof ne 0 %then %do;
  %do j=1 %to &nprof;
   %do jj=1 %to &nvar+&ntde;
    %if &&prof&j=&&vg&jj %then %do;
	 ioarray=iosave;
	 %if &profsel~= %then %do;
	  profsel=&profsel;
	  profser=&profser;
	 %end;
	 %else %do;
	  waldlo=b[&jj]-probit(1-&alpha/2)*se[&jj];
      waldup=b[&jj]+probit(1-&alpha/2)*se[&jj];
	  profsel=-(min(ci[&jj,1],waldlo)-b[&jj])/se[&jj]+0.5;
	  profser=(max(ci[&jj,2],waldup)-b[&jj])/se[&jj]+0.5;
	 %end;
	 indp=1;
	 outprof[(rby-1)*&nprof*&profn+(&j-1)*&profn+(1:3),7]=ci[&jj,1]//b[&jj]//ci[&jj,2];
     outprof[(rby-1)*&nprof*&profn+(&j-1)*&profn+(1:3),9]=(b[&jj]-probit(1-&alpha/2)*se[&jj])//b[&jj]//(b[&jj]+probit(1-&alpha/2)*se[&jj]);
     minl=0;
	 do bprof=b[&jj]-profsel*se[&jj] to b[&jj]+profser*se[&jj] by (profser+profsel)*se[&jj]/(&profn-1);
	  do jjj=1 to k;
	   ioarray[1,jjj]=2*iosave[1,jjj];
	   ioarray[2,jjj]=ioarray[3,jjj];
	  end;
      ioarray[1,&jj]=0;
      ioarray[2,&jj]=bprof*xstd[&jj];
*     print bprof ioarray iosave;
	  hugo[7]=&epsilon;
      %if &ntde ne 0 %then %do;
       ioarray[4,&nvar+1:k]={&ftmap}; 
      %end;
	  %thecall(what=fc);
	  indexp=(rby-1)*&nprof*&profn+(&j-1)*&profn+indp;
	  outprof[indexp,(1:3)]=rby||&jj||bprof;
	  %if &pl ne 3 %then %do;
 	   outprof[indexp,4]=hugo[11];
	   outprof[indexp,5]=penlike-0.5*(((bprof-b[&jj])/se[&jj])**2);
	   if min(outprof[indexp,4],outprof[indexp,5])<minl then minl=min(outprof[indexp,4],outprof[indexp,5]);
	   outprof[indexp,6]=penlike-0.5*cinv(&conflev,1);
	   indp=indp+1;
	  %end;
	  %else %do;
 	   outprof[indexp,4]=hugo[9];
	   outprof[indexp,5]=(((bprof-b[&jj])/se[&jj])**2);
	   if max(outprof[indexp,4],outprof[indexp,5])>minl then minl=max(outprof[indexp,4],outprof[indexp,5]);
	   outprof[indexp,6]=cinv(&conflev,1);
	   indp=indp+1;
	  %end;
	 end;
	 %if &pl ne 3 %then %do;
 	  outprof[(rby-1)*&nprof*&profn+(&j-1)*&profn+(1:3),8]=minl-(penlike-minl)/10;
      outprof[(rby-1)*&nprof*&profn+(&j-1)*&profn+(1:3),10]=minl-1.5*(penlike-minl)/10;
	 %end;
	 %else %do;
 	  outprof[(rby-1)*&nprof*&profn+(&j-1)*&profn+(1:3),8]=-minl/10;
      outprof[(rby-1)*&nprof*&profn+(&j-1)*&profn+(1:3),10]=-minl/10*1.5;
	 %end;
    %end;
	ioarray=iosave;
   %end;
  %end;
 %end;
 

%end;

%if &outs0 ne %then %do;
 _shelp=repeat(0,nrow(xall),4);
 indexsh=0;
 do rby=1 to maxrby;
  do i=1 to nby[rby];
   indexsh=indexsh+1;
   _shelp[indexsh,(1:4)]=rby||timeall[indexsh,2]||censall[indexsh]||0;
   do j=1 to &nvar;
    _shelp[indexsh,4]=_shelp[indexsh,4]+xall[indexsh,j]*outest[(rby-1)*(&nvar+1)+1,j+1];
   end;
  end;
 end;
%end;

create &outest from outest[colname={"_rby_" %do j=1 %to %eval(&nvar+&ntde); "&&vg&j" %end; "_penlik_" "_var_"}];
append from outest;
close &outest;

create &outtab from outtab[colname={"_rby_" "_var_" "_beta_" "_stderr_" "_lo_" "_up_" "_p_" "_bstd_"}];
append from outtab;
close &outtab;

create &outmod from outmod[colname={"_rby_" "_penlik_" "_penli0_" "_modchi_" "_p_" "_it_" "_events_" "_cens_" "_scorechi_" "_scop_" "_waldchi_" "_waldp_"}];
append from outmod;
close &outmod;

%if &ntest ne 0 %then %do;
 create &outtest from outtest[colname={"_rby_" "_penlik_" "_penli0_" "_chisqu_" "_df_" "_p_"}];
 append from outtest;
 close &outtest;
%end;

%if &nprof ne 0 %then %do;
 create &outprof from outprof[colname={"_rby_" "_var_" "_b_" "_profli_" "_normal_" "_refer_" "_bpl_" "_cipl_" "_bwald_" "_ciwald_"}];
 append from outprof;
 close &outprof;
%end;

%if &outs0 ne %then %do;
 create _shelp from _shelp[colname={"_rby_" "_time_" "_cens_" "_lp_"}];
 append from _shelp;
 close _shelp;
%end;

quit;



%if &nprof ne 0 %then %do;
 data &outprof;
 merge &outprof _fby;
 by _rby_;
 length _name_ $ 8;
 %do j=1 %to &nvar+&ntde;
  if _var_=&j then do;
   _name_="&&vg&j";
  end;
 %end;
 if _ciwald_=0 and _bwald_=0 and _bpl_=0 and _cipl_=0 then do;
  _ciwald_=.;
  _bwald_=.;
  _bpl_=.;
  _cipl_=.;
 end;
 %let clevp=%sysevalf(100*(1-&alpha));
 label _name_="Variable" _b_="beta" 
  _profli_="&claname" 
  _cipl_="&clevp.% &claname c.i." 
  _normal_="Wald" _refer_="&clevp.% reference line" _bpl_="beta" _bwald_="beta"
  _ciwald_="&clevp.% Wald c.i.";
 run;

 proc sort data=&outprof;
 by _rby_ _name_;
 run;
%end;

data &outmod;
merge &outmod _fby;
by _rby_;
run;


data &outest;
set &outest;
length _name_ $ 8 _type_ $ 8;
%do k=1 %to &nvar;
 if _var_=&k then do; 
  _name_="%upcase(&&var&k)";
  _type_="COV";
 end;
%end;
if _var_=0 then do;
 _name_="%upcase(&time)";
 _type_="PARMS";
end;
run;


data &outest;
merge &outest _fby;
by _rby_;
run;

%if &ntest ne 0 %then %do;
 data &outtest;
 merge &outtest _fby;
 by _rby_;
 run;
%end;

data &outtab;
set &outtab;
length _NAME_ $ 8 _weight_ $ 8;
%do k=1 %to &nvar;
 if &k=_VAR_ then do;
  _NAME_="%upcase(&&var&k)           ";
  if &&weit&k=1 then do;
   _WEIGHT_="%upcase(&gentype)     ";
  end;
  else do;
   _weight_="NONE    ";
  end;
 end;
%end;
%do k=1 %to &ntde;
 %let numm=&&ftmap&k;
 if _VAR_=&k+&nvar then do;
  _NAME_="%upcase(&&tdename&k)";
  if &&weit&numm=1 then do;
   _WEIGHT_="%upcase(&gentype)     ";
  end;
  else do;
   _weight_="NONE    ";
  end;
 end;
%end;
%if &risk=1 %then %do;
 risk=exp(_beta_);
 rr_lo=exp(_lo_);
 rr_up=exp(_up_);
%end;
%let conflev=%sysevalf(100*(1-&alpha));
label _var_="Variable No." _beta_="Parameter estimate" _lo_="Lower &conflev.% c.l." _up_="Upper &conflev.% c.l."
 _p_="Pr > Chi-Square" _stderr_="Standard error" _name_="Variable" _bstd_="Standardized estimate" _weight_="Weighting"
%if &risk=1 %then %do;
 risk="Risk ratio" rr_lo="Lower &conflev.% c.l." rr_up="Upper &conflev.% c.l."
%end;
;
run;

data &outtab;
merge &outtab _fby;
by _rby_;
run;

%if &outs0 ne %then %do;
 data _shc;
 do _rby_=1 to &maxby;
  _lp_=0;
  output;
 end;
 run;

 proc phreg data=_shelp noprint;
 by _rby_;
 model _time_*_cens_(0)= /offset=_lp_ maxiter=0;
 baseline out=&outs0 covariates=_shc survival=_surv_ cumhaz=_cumhaz_ loglogs=_loglogs_/ method=&s0method nomean;
 run;

 data &outs0;
 merge &outs0 _fby;
 by _rby_;
 label _time_="Survival time";
 run;
%end;

%if &print=1 %then %do;
 data &outmod;
 set &outmod;
 file print;
 if _n_=1 then do;
  put "   ";
  put "  ";
  %picture;
  put "  ";
  put "  ";
  put "  ";
  put "  Author:                     Georg Heinze";
  put "  Version:                    &version ";
  put "   ";
  put "  Documentation:              Heinze, G. (2006). "; 
  put "                              Technical Report 1/2006:";
  put "                              &mynameis.: A program for estimating Cox models";
  put "                              using Firth's penalization";
  put "                              Section of Clinical Biometrics,";
  put "                              Core Unit for Medical Statistics and Informatics,";
  put "                              Medical University of Vienna.";
  put "  ";
  put "  Data set:                   %upcase(&data)";
  %if &time ne %then %do;
   put "  Dependent variable:         %upcase(&time)";
  %end;
  %else %do;
   put "  Dependent variables:        (%upcase(&time1), %upcase(&time2)]";
  %end;
  put "  Censoring indicator:        %upcase(&cens)";
  put "  Censoring value:            &censval";
  put "  Ties handling:              Breslow";
  %if %eval(&nvar+&ntde) > 1 %then %do;
   put "  Time-invariant effects:     &varlist";
  %end;
  %else %do;
   put "  Time-invariant effect:       &varlist &tdenames";
  %end;
  %if &ntde > 0 %then %do;
   put "  Time-dependent effects:";
   %do j=1 %to &ntde;
    %let numm=&&ftmap&j;
    put "                              &&tdename&j = &&var&numm * &&ft&j ";
   %end;
  %end;
  put " ";
  put "  Table with parameter estimates saved in %upcase(&outtab).";
  %if &outest ne ' ' %then %do;
   put "  Estimates and covariance matrix saved in %upcase(&outest).";
   put "  ";
  %end;
*  put "  Covariance matrix is based on inverse Fisher information.";

 end;
 run;

 %do rby=1 %to &maxby;
  data &outmod;
   set &outmod;
   by _rby_;
   _pev_=floor(_events_/(_cens_+_events_)*1000+0.5)/10;
   _pce_=100-_pev_;
   _nobs_=_events_+_cens_;
   label _nobs_="Number of observations" _events_="Number of events" _cens_="Censored"
    _pev_="% events" _pce_="% censored"
    _penlik_="&penalcap likelihood"
    _penli0_="Null &penalize log likelihood" 
    _modchi_="Likelihood ratio Chi-square" _p_="Prob>Chi" _it_=
    "Iterations" _scorechi_="Score Chi-Square" _scop_="Prob>Chi"
    _waldchi_="Wald Chi-sqare" _waldp_="Prob>Chi";
   run;

   proc print noobs label;
   title4 "Model fitting information";
   where _rby_=&rby;
    %if &byproc=1 %then %do;
     by &by;
    %end;
    var _it_ _penlik_ _penli0_ _events_ _cens_ _nobs_ _pev_ _pce_;
	format _p_ pvalue6. _scop_ pvalue6.;
  run;

   data &global;
   set &outmod;
   df=&nvar+&ntde;
   Test="Likelihood Ratio";
   ChiSquared=_modchi_;
   ProbChi=_p_;
   output;
   Test="Scores";
   ChiSquared=_scorechi_;
   ProbChi=_scop_;
   output;
   Test="Wald";
   ChiSquared=_waldchi_;
   ProbChi=_waldp_;
   output;
   label ChiSquared="Chi-Square" ProbChi="Pr > Chi-Square" df="Degrees of freedom";
   keep &by _rby_ Test ChiSquared df ProbChi;
   run;

   proc print noobs label;
   title4 "Testing global null hypothesis: beta=0";
   where _rby_=&rby;
    %if &byproc=1 %then %do;
     by &by;
    %end;
   var Test ChiSquared df ProbChi;
   format ProbChi pvalue6.;
   run; 

  proc print data=&outtab noobs label;
   title4 "&estname estimates, &clname confidence limits";
   title5 "and &testname tests";
  %if &byproc=1 %then %do;
   where _rby_=&rby;
   by &by;
  %end;

  var _name_ 
   %if &ngv>0 %then %do;
    _weight_ 
   %end;
   _beta_ _stderr_ _bstd_ _lo_ _up_ _p_;
   format _p_ pvalue6.;
  run;

  %if &risk=1 %then %do;
   proc print data=&outtab noobs label;
    title4 "&estname risk ratio estimates, &clname confidence limits";
	title5 "and &testname tests";
   %if &byproc=1 %then %do;
    where _rby_=&rby;
    by &by;
   %end;

   var _name_ 
   %if &ngv>0 %then %do;
    _weight_ 
   %end;
    risk RR_lo RR_up _p_;
   format _p_ pvalue6.;
   run;
  %end;



  %if &ntest ne 0 %then %do;
   data &outtest;
   set &outtest;
   _tested_="&test";
   label 
     _penlik_="&penlik" _penli0_="&penli0" 
   _chisqu_="Chi-Square" _df_="Degrees of freedom" _p_="Pr > Chi-Square" _tested_="Tested parameters";
   run;

   proc print data=&outtest noobs label;
   title4 "&pentname test for parameters";
    where _rby_=&rby;
   %if &byproc=1 %then %do;
    by &by;
   %end;
   var _tested_ &restrict _chisqu_ _df_ _p_;
   format _p_ pvalue6.;
   run;
  %end;
  title4;

  %if &nprof ne 0 %then %do;
   symbol1 i=join v=none c=black line=1;
   symbol2 i=join v=none c=black line=2;
   symbol3 i=join v=plus c=black line=1;
   symbol4 i=join v=plus c=black line=3;
   legend1 down=4 label=none
   position=(bottom center outside)
   mode=reserve
   shape=symbol(10,1) 
    value=("&statval" "&clevp.% reference line" 
    "point estimate and &clevp.% &clname confidence interval" 
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

%if %upcase(&call)=EXE %then %do;
 data _ttt;
 file "&mynameis.out.bat";
 put "del delfcfiles.bat";
 put "del &mynameis..bat";
 put "del &mynameis..dll";
 put "del &mynameis.exe.exe";
 put "del copyfc.bat";
 put "exit";
 run;
 x &mynameis.out.bat;
%end;
%ende:;
%mend;
