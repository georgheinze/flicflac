# CFL.SAS: Conditional logistic regression using Firth's correction

## New option to allow computation with continuous covariates

A new option 'condition' was added which allows to include continuous covariates. If you have experienced problems with previous versions when you included continuous variables in the analysis,
you can specify the condition=NETWORKMC to solve these. This will employ a Monte Carlo algorithm to sample from the joint distribution of sufficient statistics that is needed
to compute the CFL estimates and confidence intervals. (If NETWORKMC does not work either, you can try MCMC.) The NETWORKMC option is just directly passed to PROC LOGISTIC's EXACTOPTIONS
statement which controls the computation of the joint distribution.

See the example below:

```
data simhard;
do strata =1 to 20;
	do i=1 to 4;
		if i=1 then y=1;
		else y=0;
		x=ranbin(58301,1,0.2);
		z=rannor(583818);
		output;
	end;
end;
run;

data simhard;
set simhard;
z1=round(z,0.1);
run;

* will not run: condition=DIRECT has always been the default for CFL;
%cfl(data=simhard, y=y, strata=strata, varlist=x z, condition=DIRECT);


* new option condition=NETWORKMC to compute the joint distribution (condition in strata) by the NETWORK-Monte-Carlo method;
%cfl(data=simhard, y=y, strata=strata, varlist=x z, condition=NETWORKMC);
```

Another way to treat the problem is by rounding of the continuous covariate. However, this may result in loss of information.

```
* default works but with rounded continuous variable z1;
%cfl(data=simhard, y=y, strata=strata, varlist=x z1, condition=DIRECT);
```

Further new features include an automatic erasure of previous results tables such that if an error occurs, the macro does not show you 'old' results.

