
*
*	This file contains code to work with synthetic data
*

capture: program drop makeJointData
program define makeJointData, rclass

	syntax [, Groups(integer 5000) Size(integer 7) Varnoise(real -1) Rho1pars(numlist min=2 max=2 >=0.05 <=1e5) CLEAR INDEPendent]
	
	if (`varnoise'<0) {
		local varnoise 1/18
	}
	tempname alpha beta
	if ("`rho1pars'"=="") {
		scalar `alpha' = 1e5
		scalar `beta' = 0.15
	}
	else {
		scalar `alpha' = max(0.05,min(1e5,`: word 1 of `rho1pars''))
		scalar `beta' = max(0.15,min(1e5,`: word 2 of `rho1pars''))
	}
		
	tempname halfwidth
	scalar `halfwidth' = sqrt(3*max(0,`varnoise'))
	
	`clear'
	set obs `=`groups'*`size''
	
	*	Arrival is by definition uniform over (0,1) but it is not actually random
	tempname sajk
	gen ARRIVAL = (_n-1)/(_N-1)
	sum ARRIVAL
	scalar `sajk' = r(sd)
	
	*	We want to pseudo-sort on ARRIVAL
	tempvar sort_noise sort_var
	gen `sort_noise' = runiform(-`halfwidth', `halfwidth')
	gen `sort_var' = ARRIVAL + `sort_noise'
	sort `sort_var'
	gen sheet = _n
	gen precinct = floor((sheet-1)/`size')
	
	*	Some diagnostic information
	tempname max_range
	tempvar amin amax arange
	bys precinct: egen `amin' = min(ARRIVAL)
	bys precinct: egen `amax' = max(ARRIVAL)
	gen `arange' = `amax'-`amin'
	sum `arange'
	scalar `max_range' = r(max)
	
	*	Group means and deviations
	tempname sak srjk
	tempvar ak rjk
	bys precinct: egen `ak' = mean(ARRIVAL)
	gen `rjk' = ARRIVAL-`ak'
	sum `ak'
	scalar `sak' = r(sd)
	sum `rjk'
	scalar `srjk' = r(sd)
		
	*	Maximum available rho
	tempname rho_max
	scalar `rho_max' = `sak'/`sajk'
		
	*	Pick a rho
	tempname rho1 rho
	if ("`independent'"=="independent") {
		scalar `rho1' = 0
	}
	else {
		scalar `rho1' = 2*rbeta(`alpha',`beta')-1
	}
	di `rho1'
	scalar `rho' = `rho1'*`rho_max'
		
	*	Initial vector
	*		must be constant within precinct
	*		currently uniform, mean zero, unit variance
	gen mu = runiform(-sqrt(3),sqrt(3))
	bys precinct: replace mu = mu[1]
	
	*	Take part of vector orthogonal to group means of ARRIVAL
	tempname sek sekrjk
	tempvar ek ekrjk
	reg mu `ak'
	predict `ek', resid
	gen `ekrjk' = `ek'*`rjk'
	sum `ek'
	scalar `sek' = r(sd)
	sum `ekrjk'
	scalar `sekrjk' = r(sd)
	replace mu = `rho1'*`ak'+sqrt(1-`rho1'^2)*`sak'*(`ek'/`sek')
	
	*	Set the mean to 1/2
	sum mu
	replace mu = mu+(1/2-r(mean))
	
	correl `ak' `ek' `rjk'
	
	reg mu ARR
	di `rho'*`rho'/`rho1'
	
	return scalar rho_max = `rho_max'
	return scalar rho = `rho'
	return scalar noise = `varnoise'
	
	di `alpha'
	di `beta'
	
end

capture: program drop addResults
program define addResults, rclass

	syntax [, B5(real 0) noNoise noSlope SHUTdown(numlist min=1 max=2) Means(numlist min=1 max=4 >0 <100) CDMax]

	if (`b5'<-100 | `b5'>100) {
		local b5 0
	}
	if ("`means'"=="") {
		local means 45 40
	}
	tempname s0 s1
	if ("`shutdown'"=="") {
		scalar `s0' = 0
		scalar `s1' = 0
	}
	else {
		scalar `s0' = `: word 1 of `shutdown''
		if (`:list sizeof shutdown'>1) {
			scalar `s1' = `: word 2 of `shutdown''
		}
		else {
			scalar `s1' = 0
		}
	}
	di `s0'
	di `s1'
	
	gen SHUTDOWN = 6*ARRIVAL>5

	*	Simulate turnout
	tempvar turnout_base
	gen `turnout_base' = rbeta(16,4)
	bys precinct: gen turnout = rbinomial(250,`turnout_base'[1])
	
	*	Simulate SI/MAS share
	tempname beta1 A B C D b0 b1 b2 b3 bp bm br p1 m1 cmin cmax dmin dmax
	
	reg mu ARRIVAL
	lincom ARRIVAL
	scalar `beta1' = r(estimate)
	
	if ("`slope'"=="noslope" | abs(`beta1')<1e-5) {
		scalar `bp' = 1
		scalar `bm' = 0
		scalar `br' = 1
	}
	else {
		scalar `bp' = (1+abs(`beta1'))/(2*abs(`beta1'))
		scalar `bm' = (1-abs(`beta1'))/(2*abs(`beta1'))
		scalar `br' = (1+abs(`beta1'))/(1-abs(`beta1'))
	}

	scalar `A' = `: word 1 of `means''
	scalar `B' = `: word 2 of `means''
		
	scalar `A' = max(0,min(100,`A'))
	scalar `B' = max(max(0,-`b5'),min(100+min(0,-`b5'),`B'))
		
	scalar `cmin' = max(`A'/`br'                    , `A'*`br'                    -100/`bm')
	scalar `cmax' = min(`A'*`br'                    , `A'/`br'                    +100/`bp')
	scalar `dmin' = max(`B'/`br'+max(`b5',`b5'/`br'), `B'*`br'+max(`b5',`b5'*`br')-100/`bm')
	scalar `dmax' = min(`B'*`br'+min(`b5',`b5'*`br'), `B'/`br'+min(`b5',`b5'/`br')+100/`bp')
	
	if ("`slope'"=="noslope") {
		scalar `C' = `A'
		scalar `D' = `B'+`b5'
	}
	else {
		scalar `C' = 55
		scalar `D' = 60+`b5'/`beta1'
		if ("`cdmax'"=="cdmax") {
			scalar `C' = `dmax'-5
			scalar `D' = `dmax'
		}
		else {
			if (`: list sizeof means'>2) {
				scalar `C' = `:word 3 of `means''
				if (`: list sizeof means'>3) {
					scalar `D' = `:word 4 of `means''
				}
			}
		}
	}
	scalar `C' = max(`cmin',min(`cmax',`C'))
	scalar `D' = max(`dmin',min(`dmax',`D'))

	scalar `b0' = `bp'*`A'-`bm'*`C'
	scalar `b1' = `bp'*(`B'-`A')-`bm'*(`D'-`C'-`b5')
	scalar `b2' =  (`C'-`A')/`beta1'
	scalar `b3' = ((`D'-`B')-(`C'-`A')-`b5')/`beta1'
	
	if ("`slope'"=="noslope") {
		scalar `b3' = 0
		scalar `b2' = 15
		scalar `b0' = `A'-`b2'/2
		scalar `b1' = `B'-`A'
	}
			
	*	mean is y_ijk-e_ijk
	gen mean2016 =  `b0'      + `b2'      *mu			  +(`s0'+`s1'*ARRIVAL)*SHUTDOWN/2
	gen mean2019 = (`b0'+`b1')+(`b2'+`b3')*mu+`b5'*ARRIVAL+(`s0'+`s1'*ARRIVAL)*SHUTDOWN
	reshape long mean, i(sheet) j(Year)
	*	share is y_ijk
	if ("`noise'"=="nonoise") {
		gen Share = mean
	}
	else {
		gen Share = 100*rbinomial(turnout,mean/100)/turnout
	}
	gen e = Share-mean

	tempname BA Bhat B0 B1 B5 B5e mum em e1m e0m am sa samu sae1 sae0 rho
	sum mu
	scalar `mum' = r(mean)
	sum e, meanonly
	scalar `em' = r(mean)
	sum e if Year==2019, meanonly
	scalar `e1m' = r(mean)
	sum e if Year==2016, meanonly
	scalar `e0m' = r(mean)
	sum ARRIVAL
	scalar `am' = r(mean)
	scalar `sa' = r(sd)
	correl ARRIVAL mu, cov
	scalar `samu' = r(cov_12)
	correl ARRIVAL e if Year==2019, cov
	scalar `sae1' = r(cov_12)
	correl ARRIVAL e if Year==2016, cov
	scalar `sae0' = r(cov_12)
	
	scalar `B0' = `b0'+`b2'*`mum'
	scalar `B5' = `b3'*`samu'/(`sa'*`sa')
	scalar `B1' = `b1'+`b3'*`mum'-`B5'*`am'
	mat `Bhat' = (`B1' \ `B5'+`b5' \ `B0')

	scalar `B5e' = (`sae1'-`sae0')/(`sa'*`sa')
	scalar `B0' = `B0'+`em'-(`e1m'-`e0m')/2
	scalar `B1' = `B1'-`B5e'*`am'+(`e1m'-`e0m')
	scalar `B5' = `B5'+`B5e'
	mat `BA' = (`B1' \ `B5'+`b5' \ `B0')
	di
	di as text " A: " as result `A'
	di as text " B: " as result `B'
	di as text " C: " as result `C'
	di as text " D: " as result `D'
	di
	di as text "b0: " as result `b0'
	di as text "b1: " as result `b1'
	di as text "b2: " as result `b2'
	di as text "b3: " as result `b3'
	di as text "b5: " as result `b5'
	di
	di as text "b0            : " as result `b0'
	di as text "b0   +b2      : " as result `b0'     +`b2'
	di as text "b0+b1+b2+b3+b5: " as result `b0'+`b1'+`b2'+`b3'+`b5'
	di as text "b0+b1      +b5: " as result `b0'+`b1'          +`b5'
	di as text "b0+b1         : " as result `b0'+`b1'
	di
	di as text "Expected coefficient vector:"
	mat li `Bhat'
	di
	di as text "Estimated coefficient vector:"
	mat li `BA'
	
	return scalar D0 = 60+`b5'/`beta1'
	
	return scalar b5 = `b5'
	return scalar b3 = `b3'
	return scalar b2 = `b2'
	return scalar b1 = `b1'
	return scalar b0 = `b0'
	
	return scalar D = `D'
	return scalar C = `C'
	return scalar B = `B'
	return scalar A = `A'

	return scalar beta1 = `beta1'

end

local newgraphs

set seed 20210615
local rngstate = c(rngstate)
local c1 ceprgreen%25
local c2 ceprorange%25

local c1l ceprgreen%10
local c2l ceprorange%10

tempfile geo nongeo

*
*	Data Generation
*
*		note we will need basically two joint data draws
*			one will have arrival be unrelated to geograpghy
*
*	Figure 11
*
*	Connected to geography
clear
makeJointData
save `geo'
*	Generic response to geography
addResults
keep if Year==2019
bys precinct: egen amax = max(ARR)
bys precinct: egen amin = min(ARR)
gen ara = amax-amin
sum ara
levelsof prec if ara==r(max), local(mp) clean
capture: confirm file ${FIGS}/Fig11.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if prec==`mp', mc(ceprorange)) ///
		, xti(" " "ARRIVAL") yti("Vote Share") legend(off) aspect(1)
	graph export ${FIGS}/Fig11.png, replace
}
*
*	Figure 12
*
capture: confirm file ${FIGS}/Fig12.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter Share ARR, mc(`c2') msize(vtiny)) ///
		(scatter Share ARR if prec==`mp', mc(ceprorange)) ///
		, xti(" " "ARRIVAL") yti("Vote Share") legend(off)  aspect(1)
	graph export ${FIGS}/Fig12.png, replace
}
*
*	Figure 13
*
use `geo', clear
*	With added trend
qui addResults, b5(20) means(45 40 55 85)
ret li
scalar rD = r(D)
bys sheet: egen m16 = total(cond(Year==2016,mean,.))
keep if Year==2019
egen mcut = cut(m16), g(11)
bys mcut: egen mm = mean(mean)
bys precinct: egen amax = max(ARR)
bys precinct: egen amin = min(ARR)
gen ara = amax-amin
sum ara
levelsof prec if ara==r(max), local(mp) clean
sum m16 if precinct==`mp'
local mlo = floor(r(min)/2.5)*2.5
local mhi = ceil(r(max)/2.5)*2.5
macro li _mlo _mhi
capture: confirm file ${FIGS}/Fig13.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if precinct==`mp', mc(ceprorange)) ///
		, name(yesb2, replace) xti("") yti("") legend(off)
}
capture: confirm file ${FIGS}/FigC.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if precinct==`mp', mc(ceprorange)) ///
		if inrange(m16,`mlo',`mhi'), name(ryesb2, replace) xti("") yti("") legend(off) ylab(20(20)100)
}

use `geo', clear
*	Same response with no trend
qui addResults, means(45 40 55 `=rD')
ret li
if (abs(r(D)-rD)>0.0001) {
	di as error "Oh noes!!!"
	mata: VVV
}
bys sheet: egen m16 = total(cond(Year==2016,mean,.))
keep if Year==2019
egen mcut = cut(m16), g(4)
bys precinct: egen amax = max(ARR)
bys precinct: egen amin = min(ARR)
gen ara = amax-amin
reg Share i1.SHUT, vce(robust)
reg Share i1.SHUT#ib(0).mcut ib(0).mcut, vce(robust)
predict shat
/*
forvalues ii=0/3 {
	sum ara if mcut==`ii'
	levelsof prec if ara==r(max), local(mp) clean
	local mp: word 1 of `mp'
	twoway ///
		(scatter Share ARR if ~SHUT, mc(`c2') msize(vtiny)) ///
		(scatter Share ARR if SHUT, mc(`c1') msize(vtiny)) ///
		(line shat ARR if ~SHUT, lc(ceprorange)) ///
		(line shat ARR if SHUT, lc(ceprgreen)) ///
		(scatter Share ARR if ~SHUT & precinct==`mp', mc(ceprorange)) ///
		(scatter Share ARR if SHUT & precinct==`mp', mc(ceprgreen)) ///
		if mcut==`ii', name(x`ii') xti("") yti("") legend(off)
}
graph combine x0 x1 x2 x3, ycomm xcomm
mata: CVVV
*/
sum mean if inrange(m16,`mlo',`mhi')
capture: confirm file ${FIGS}/Fig13.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if precinct==`mp', mc(ceprorange)) ///
		, name(nob2, replace) xti("") yti("") legend(off)
	graph combine nob2 yesb2, ycommon b1title("ARRIVAL") l1title("Vote Share")
	graph export ${FIGS}/Fig13.png, replace
}

capture: confirm file ${FIGS}/FigC.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if precinct==`mp', mc(ceprorange)) ///
		if inrange(m16,`mlo',`mhi'), name(rnob2, replace) xti("") yti("") legend(off) ylab(20(20)100)
	graph combine rnob2 ryesb2, ycommon b1title("ARRIVAL") l1title("Vote Share")
	graph export ${FIGS}/FigC.png, replace
}

*
*	Table 2: Difference Models
*	Column 2
*
*	Figure 14
*
*	Not connected to geography
makeJointData, indep clear
save `nongeo'
*	No voting response either
qui addResults, noslope
ret li
keep if Year==2019
bys precinct: egen amax = max(ARR)
bys precinct: egen amin = min(ARR)
gen ara = amax-amin
sum ara
levelsof prec if ara==r(max), local(mp2) clean
capture: confirm file ${FIGS}/Fig14.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR if ~SHUT, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if  SHUT, mc(`c1') msize(vtiny)) ///
		(scatter mean ARR if ~SHUT & precinct==`mp2', mc(ceprorange)) ///
		(scatter mean ARR if  SHUT & precinct==`mp2', mc(ceprgreen)) ///
		, legend(off) xti(" " "ARRIVAL") yti("Vote Share") ylab(20(10)60) aspect(1)
	graph export ${FIGS}/Fig14.png, replace
}
*
*	Figure 15
*
use `nongeo', clear
*	Add shutdown break
qui addResults, nos shut(5)
ret li
keep if Year==2019
capture: confirm file ${FIGS}/Fig15.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR if ~SHUT, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if  SHUT, mc(`c1') msize(vtiny)) ///
		(scatter mean ARR if ~SHUT & precinct==`mp2', mc(ceprorange)) ///
		(scatter mean ARR if  SHUT & precinct==`mp2', mc(ceprgreen)) ///
		, legend(off) yti("Vote Share") ylab(20(10)60) aspect(1)
	graph export ${FIGS}/Fig15.png, replace 
}
*
*	Figure 16
*
reg mean i1.SHUT
predict model
sort ARR
capture: confirm file ${FIGS}/Fig16.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR if ~SHUT, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if  SHUT, mc(`c1') msize(vtiny)) ///
		(scatter mean ARR if ~SHUT & precinct==`mp2', mc(ceprorange)) ///
		(scatter mean ARR if  SHUT & precinct==`mp2', mc(ceprgreen)) ///
		(line model ARR if ~SHUT, pstyle(p2)) ///
		(line model ARR if  SHUT, pstyle(p4)) ///
		, legend(off) yti("Vote Share") ylab(20(10)60) aspect(1)
	graph export ${FIGS}/Fig16.png, replace
}

capture: program drop addEstimates
program define addEstimates

	syntax, Est(numlist) SE(numlist) DF(numlist) Colors(string asis) ci(cilevel) LABels(string asis) XTItle(passthru)
	
	clear
	
	local ne: list sizeof est
	
	tempname lb ub
	scalar `lb' = `:word 1 of `est''
	scalar `ub' = `lb'
	forvalues ii=1/`ne' {
		scalar `lb' = min(`lb',`:word `ii' of `est''+invt(`:word `ii' of `df'',(100-`ci')/200)*`:word `ii' of `se'')
		scalar `ub' = max(`ub',`:word `ii' of `est''+invt(`:word `ii' of `df'',1-(100-`ci')/200)*`:word `ii' of `se'')
	}
	macro li _ci
	di `lb'
	di `ub'
	
	set obs 1000
	gen yx = .
	gen xx = 0
	gen x0 = `lb'+(`ub'-`lb')*(_n-1)/_N
	gen x1 = `lb'+(`ub'-`lb')*(_n-0)/_N
	gen xp = (x0+x1)/2
	gen xpr = (xp)*(827823+153890)/6137671
	tempname dmax
	scalar `dmax' = 0
	forvalues ii=1/`ne' {
		gen d`ii' = t(`:word `ii' of `df'',(x1-`:word `ii' of `est'')/`:word `ii' of `se'') ///
			-t(`:word `ii' of `df'',(x0-`:word `ii' of `est'')/`:word `ii' of `se'')
		sum d`ii', meanonly
		scalar `dmax' = max(`dmax',r(max))
		tempvar hi`ii' lo`ii'
		gen `hi`ii'' = -`ii'+1.2
		gen `lo`ii'' = -`ii'+0.8
		replace `hi`ii'' = 0.5-`ii'+1*tden(`:word `ii' of `df'',(xp-`:word `ii' of `est'')/`:word `ii' of `se'')/`:word `ii' of `se''
		replace `lo`ii'' = 0.5-`ii'
	}

	local lcoms (scatter yx xx)
	local ylabs ylab(none, noticks)
	tempvar xh xl nj
	gen `xh'=.
	gen `xl'=.
	gen `nj' = 1-_n
	forvalues jj=1/`ne' {
		forvalues ii=1/100 {
			local lcoms `lcoms' (rarea `hi`jj'' `lo`jj'' xp if d`jj'>(`ii'/100)*`dmax', lw(none) color(`:word `jj' of `colors''%2))
		}
		/*
		replace `xl'= `:word `jj' of `est''+invt(`:word `jj' of `df'',0.025)*`:word `jj' of `se'' in `jj'
		replace `xh'= `:word `jj' of `est''+invt(`:word `jj' of `df'',0.975)*`:word `jj' of `se'' in `jj'
		local lcoms `lcoms' (rcap `xl' `xh' `nj' in `jj', horiz lw(thin) color(`:word `jj' of `colors''))
		*/
		local ylabs `ylabs' ylab(`=1-`jj'' "`:word `jj' of `labels''", add custom labcolor(`:word `jj' of `colors''))
	}
	twoway `lcoms', legend(off) ysc(r(`=-(`ne'-1)-0.25' 0.25)) `ylabs' `xtitle' name(efs, replace) aspect() ///
		ti("Escobari and Hoover Table 2: “Fraud” Coefficients") note("Note: Estimate (1) sign reversed")

end


*
*	Figure 18
*
use `geo', clear
qui addResults, b5(0) means(45 40 55 60)
ret li
scalar Delta = (r(beta1)*(r(b2)+r(b3))+r(b5))/2
di as text "Delta: " as result Delta
keep if Year==2019
reg mean i1.SHUT
lincom 1.SHUT
scalar est = r(estimate)
scalar se = r(se)
scalar df = r(df)
predict model
reg mean c.ARR##i1.SHUT
predict trend
sort ARR
capture: confirm file ${FIGS}/Fig18.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR if ~SHUT, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if  SHUT, mc(`c1') msize(vtiny)) ///
		(scatter mean ARR if ~SHUT & precinct==`mp', mc(ceprorange)) ///
		(scatter mean ARR if  SHUT & precinct==`mp', mc(ceprgreen)) ///
		(line trend ARR if ~SHUT, pstyle(p2) lp(dash) lw(thin)) ///
		(line trend ARR if SHUT, pstyle(p4) lp(dash) lw(thin)) ///
		(line model ARR if ~SHUT, pstyle(p2)) ///
		(line model ARR if  SHUT, pstyle(p4)) ///
		, legend(off) xti(" " "ARRIVAL") yti("Vote Share") aspect(1)
	/*
	addEstimates, e(8.286 7.975 16.26 7.243 6.672 0.365 `=est') ///
		se(0.324 0.343 0.653 0.437 0.464 0.194 `=se') ///
		df(34527 34527 34527 34527 34527 34527 `=df') ///
		colors(ceprorange ceprgreen ceprblue ceprblue ceprblue ceprblue ceprorange) ///
		labels("CC (1)" "MAS (2)" "MAS-CC (3)" "MAS-CC (4)" "MAS-CC (5)" "MAS-CC (6)" "Synthetic") ///
		xti(" " "SHUTDOWN") ci(99.99)
	graph export ~/Documents/ptest.png, replace wid(2000)
	mata: VVV
	*/
	graph export ${FIGS}/Fig18.png, replace
}
use `geo', clear
qui addResults, b5(0) means(45 40 55 60)
ret li
scalar Delta = (r(beta1)*(r(b2)+r(b3))+r(b5))/2
di as text "Delta: " as result Delta
bys sheet: egen m16 = total(cond(Year==2016,mean,.))
keep if Year==2019
egen mcut = cut(m16), g(4)
reg Share i1.SHUT, vce(robust)
reg Share i1.SHUT#ib(0).mcut ib(0).mcut, vce(robust)
predict model
reg Share c.ARR##i1.SHUT
predict trend
sort ARR
*
*	Figure 19
*
capture: confirm file ${FIGS}/Fig19.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	gen shut2 = .
	forvalues aa=1(2)7 {
		capture: drop model trend
		replace shut2 = ARR>`aa'/10
		reg mean i1.shut2
		predict model
		reg mean c.ARR##i1.shut2
		predict trend
		sort ARR
		twoway ///
			(scatter mean ARR if ~shut2, mc(`c2l') msize(vtiny)) ///
			(scatter mean ARR if  shut2, mc(`c1l') msize(vtiny)) ///
			(scatter mean ARR if ~shut2 & precinct==`mp', mc(ceprorange)) ///
			(scatter mean ARR if  shut2 & precinct==`mp', mc(ceprgreen)) ///
			(line trend ARR if ~shut2, pstyle(p2) lp(dash) lw(thin)) ///
			(line trend ARR if shut2, pstyle(p4) lp(dash) lw(thin)) ///
			(line model ARR if ~shut2, pstyle(p2)) ///
			(line model ARR if  shut2, pstyle(p4)) ///
			, legend(off) yti("") xti("") name(fig14_`aa') aspect(1)
	}
	graph combine fig14_1 fig14_3 fig14_5 fig14_7, xcommon ycommon b1t("ARRIVAL") l1t("Vote Share")
	graph export ${FIGS}/Fig19.png, replace
}
*
*	Table 2
*	Column 6
*
*	Figure 20
*
capture: drop model
areg mean i1.SHUT, abs(precinct)
predict fe, d
gen ashare = Share-fe
gen amean = mean-fe
reg amean c.ARR##i1.SHUT
predict model
sort ARR
capture: confirm file ${FIGS}/Fig20.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter ashare ARR if ~SHUT, mc(`c2') msize(vtiny)) ///
		(scatter ashare ARR if  SHUT, mc(`c1') msize(vtiny)) ///
		(line model ARR if ~SHUT, pstyle(p2)) ///
		(line model ARR if  SHUT, pstyle(p4)) ///
		, legend(off) xti(" " "ARRIVAL") yti("(Precinct-Adjusted) Vote Share") aspect(1)
	graph export ${FIGS}/Fig20.png, replace
}
*
*	Table 3: Difference-In-Difference Models
*	Column 5
*
*	Figure 21
*
use `nongeo', clear
qui addResults, nos shut(4)
table SHUT Year, c(mean mean)
reg mean i1.SHUT##i.Year
predict model
sort ARR
capture: confirm file ${FIGS}/Fig21.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR if Year==2019, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if Year==2016, mc(`c1') msize(vtiny)) ///
		(scatter mean ARR if Year==2019 & precinct==`mp2', mc(ceprorange)) ///
		(scatter mean ARR if Year==2016 & precinct==`mp2', mc(ceprgreen)) ///
		(line model ARR if Year==2019 & SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & SHUT, pstyle(p4)) ///
		(line model ARR if Year==2019 & ~SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & ~SHUT, pstyle(p4)) ///
		, legend(order(6 5) label(5 "2019") label(6 "2016")) xti(" " "ARRIVAL") yti("Vote Share") ylab(20(10)60) aspect(1)
	graph export ${FIGS}/Fig21.png, replace
}
*
*	Figure 22
*
use `geo', clear
qui addResults
ret li
reg mean i1.SHUT##i.Year
predict model
reg mean c.ARR##i.Year
predict trend
sort ARR
capture: confirm file ${FIGS}/Fig22.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter mean ARR if Year==2019, mc(`c2') msize(vtiny)) ///
		(scatter mean ARR if Year==2016, mc(`c1') msize(vtiny)) ///
		(scatter mean ARR if Year==2019 & precinct==`mp', mc(ceprorange)) ///
		(scatter mean ARR if Year==2016 & precinct==`mp', mc(ceprgreen)) ///
		(line trend ARR if Year==2019, pstyle(p2) lp(dash) lw(thin)) ///
		(line trend ARR if Year==2016, pstyle(p4) lp(dash) lw(thin)) ///
		(line model ARR if Year==2019 & SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & SHUT, pstyle(p4)) ///
		(line model ARR if Year==2019 & ~SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & ~SHUT, pstyle(p4)) ///
		, legend(order(8 7) label(7 "2019") label(8 "2016")) xti(" " "ARRIVAL") yti("Vote Share") aspect(1)
	graph export ${FIGS}/Fig22.png, replace
}
*
*	Table 3
*	Column 6
*
*	Figure 22
*
capture: drop model trend
areg mean i1.SHUT##i.Year, abs(precinct)
predict model
predict fe, d
gen amean = mean-fe
reg amean c.ARR##i.Year
predict trend
sort ARR
capture: confirm file ${FIGS}/Fig23.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter amean ARR if Year==2019, mc(`c2') msize(vtiny)) ///
		(scatter amean ARR if Year==2016, mc(`c1') msize(vtiny)) ///
		(scatter amean ARR if Year==2019 & precinct==`mp', mc(ceprorange)) ///
		(scatter amean ARR if Year==2016 & precinct==`mp', mc(ceprgreen)) ///
		(line trend ARR if Year==2019, pstyle(p2) lp(dash) lw(thin)) ///
		(line trend ARR if Year==2016, pstyle(p4) lp(dash) lw(thin)) ///
		(line model ARR if Year==2019 & SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & SHUT, pstyle(p4)) ///
		(line model ARR if Year==2019 & ~SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & ~SHUT, pstyle(p4)) ///
		, legend(order(8 7) label(7 "2019") label(8 "2016")) xti(" " "ARRIVAL")  yti("(Precinct-Adjusted) Vote Share") ylab(45(5)55) aspect(1)
	graph export ${FIGS}/Fig23.png, replace
}
*
*	Table 3
*	Column 4
*
*	Figure 24
*
use `geo', clear
qui addResults, b5(0) means(45 40 55 60)
ret li
areg mean i1.SHUT##i.Year, abs(sheet)
predict model
predict fe, d
gen amean = mean-fe
reg amean c.ARR##i1.SHUT##i.Year
predict trend
sort ARR
capture: confirm file ${FIGS}/Fig24.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter amean ARR if Year==2019, mc(`c2') msize(vtiny)) ///
		(scatter amean ARR if Year==2016, mc(`c1') msize(vtiny)) ///
		(scatter amean ARR if Year==2019 & precinct==`mp', mc(ceprorange)) ///
		(scatter amean ARR if Year==2016 & precinct==`mp', mc(ceprgreen)) ///
		(line trend ARR if Year==2019, pstyle(p2) lp(dash) lw(thin)) ///
		(line trend ARR if Year==2016, pstyle(p4) lp(dash) lw(thin)) ///
		(line model ARR if Year==2019 & SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & SHUT, pstyle(p4)) ///
		(line model ARR if Year==2019 & ~SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & ~SHUT, pstyle(p4)) ///
		, legend(order(8 7) label(7 "2019") label(8 "2016")) xti(" " "ARRIVAL") yti("(Station-Adjusted) Vote Share") ylab(45(5)55) aspect(1)
	graph export ${FIGS}/Fig24.png, replace
}
*
*	Figure 25
*
reg fe c.ARR##i1.SHUT if Year==2019
predict fehat
capture: confirm file ${FIGS}/Fig25.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter fe ARR if ~SHUT, mc(`c2') msize(vtiny)) ///
		(scatter fe ARR if SHUT, mc(`c1') msize(vtiny)) ///
		(line fehat ARR if ~SHUT, pstyle(p2)) ///
		(line fehat ARR if SHUT, pstyle(p4)) ///
		if Year==2019, legend(off) xti(" " "ARRIVAL") yti("Polling-Station Effect") aspect(1)
	graph export ${FIGS}/Fig25.png, replace
}
*
*	Table 4
*	Column 3
*
*	Figure 26
*
use `geo', clear
qui addResults, b5(5) means(45 40 55 60)
ret li
areg mean i1.SHUT##i.Year, abs(sheet)
predict model
predict fe, d
gen amean = mean-fe
reg amean c.ARR##i1.SHUT##i.Year
predict trend
sort ARR
capture: confirm file ${FIGS}/Fig26.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter amean ARR if Year==2019, mc(`c2') msize(vtiny)) ///
		(scatter amean ARR if Year==2016, mc(`c1') msize(vtiny)) ///
		(scatter amean ARR if Year==2019 & precinct==`mp', mc(ceprorange)) ///
		(scatter amean ARR if Year==2016 & precinct==`mp', mc(ceprgreen)) ///
		(line trend ARR if Year==2019, pstyle(p2) lp(dash) lw(thin)) ///
		(line trend ARR if Year==2016, pstyle(p4) lp(dash) lw(thin)) ///
		(line model ARR if Year==2019 & SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & SHUT, pstyle(p4)) ///
		(line model ARR if Year==2019 & ~SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & ~SHUT, pstyle(p4)) ///
		, legend(order(8 7) label(7 "2019") label(8 "2016")) xti(" " "ARRIVAL") yti("(Station-Adjusted) Vote Share") ylab(45(5)55) aspect(1)
	graph export ${FIGS}/Fig26.png, replace
}
*
*	Column 4
*
*	Figure 27
*
use `geo', clear
qui addResults, b5(0) means(45 40 55 60)
ret li
areg mean 1.SHUT#2019.Year 2019.Year c.ARR#1.SHUT#2019.Year, abs(sheet)
predict model
predict fe, d
gen amean = mean-fe
reg amean c.ARR##i1.SHUT##i.Year
predict trend
sort ARR
capture: confirm file ${FIGS}/Fig27.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter amean ARR if Year==2019, mc(`c2') msize(vtiny)) ///
		(scatter amean ARR if Year==2016, mc(`c1') msize(vtiny)) ///
		(scatter amean ARR if Year==2019 & precinct==`mp', mc(ceprorange)) ///
		(scatter amean ARR if Year==2016 & precinct==`mp', mc(ceprgreen)) ///
		(line trend ARR if Year==2019, pstyle(p2) lp(dash) lw(thin)) ///
		(line trend ARR if Year==2016, pstyle(p4) lp(dash) lw(thin)) ///
		(line model ARR if Year==2019 & SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & SHUT, pstyle(p4)) ///
		(line model ARR if Year==2019 & ~SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & ~SHUT, pstyle(p4)) ///
		, legend(order(8 7) label(7 "2019") label(8 "2016")) xti(" " "ARRIVAL") yti("(Station-Adjusted) Vote Share") ylab(45(5)55) aspect(1)
	graph export ${FIGS}/Fig27.png, replace
}
*
*	Table 5
*	Column 1
*
*	Figure 28
*
capture: drop model fe amean trend
areg mean 1.SHUT#2019.Year c.ARR#2019.Year 2019.Year, abs(sheet)
predict model
predict fe, d
gen amean = mean-fe
reg amean c.ARR##i1.SHUT##i.Year
predict trend
sort ARR
capture: confirm file ${FIGS}/Fig28.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter amean ARR if Year==2019, mc(`c2') msize(vtiny)) ///
		(scatter amean ARR if Year==2016, mc(`c1') msize(vtiny)) ///
		(scatter amean ARR if Year==2019 & precinct==`mp', mc(ceprorange)) ///
		(scatter amean ARR if Year==2016 & precinct==`mp', mc(ceprgreen)) ///
		(line trend ARR if Year==2019, pstyle(p2) lp(dash) lw(thin)) ///
		(line trend ARR if Year==2016, pstyle(p4) lp(dash) lw(thin)) ///
		(line model ARR if Year==2019 & SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & SHUT, pstyle(p4)) ///
		(line model ARR if Year==2019 & ~SHUT, pstyle(p2)) ///
		(line model ARR if Year==2016 & ~SHUT, pstyle(p4)) ///
		, legend(order(8 7) label(7 "2019") label(8 "2016")) xti(" " "ARRIVAL") yti("(Station-Adjusted) Vote Share") ylab(45(5)55) aspect(1)
	graph export ${FIGS}/Fig28.png, replace
}
