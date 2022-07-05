set more off
log close _all
clear all
set type double
set maxvar 10000
*	- help scheme entries - 
set scheme ceprrebrand2
*	This program will seek and create files in the current working directory
*	Please adjust as desired and ensure that the contents of the LocalSources directory on GitHub are in the subdirectory indicated as "RAW" below 
cd ~/Documents/GitHub/Bolivia-Escobari-Hoover/

*	For raw source data not simply available on web
global RAW LocalSources
*	For raw source data
global DATA Data
*	For intermediate files
global INPUT Input
*	For final results
global OUTPUT Results
global FIGS ${OUTPUT}/Figures

capture: confirm file ${DATA}
if (_rc~=0) {
	*	Ensure that the data in LocalSources is available in the raw source directory
	*	This is to ensure that less-easily obtainable data is available if the RAW (Data), INPUT, and OUTPUT (Results) directories need to be reconstructed)
	*	 - That is, the full results should be replicable from only the LocalSources directory and its contents, plus the Stata do-files, and internet connectivity
	! cp -R ${RAW} ${DATA}
}
capture: mkdir ${INPUT}
capture: mkdir ${OUTPUT}
capture: mkdir ${FIGS}

include newEscobari_utilities
include newEscobari_matching
include newEscobari_construction
include newEscobari_monteCarlo

capture: program drop buildBaseData
program define buildBaseData

	syntax [, REPLACE]
	
	capture: use ${OUTPUT}/baseData, clear
	if (_rc~=0 | "`replace'"=="replace") {
	
		tempfile gnt
		getResumptionTSE
		keep NúmeroMesa
		save `gnt'
	
		buildFullMatches
		mergePadronFinal
		addMoreGeo
		addPublicTiming
		
		merge 1:1 NúmeroMesa using `gnt'
		gen inResumption = _merge==3
		drop _merge
		
		addTSE
		compress
		save ${OUTPUT}/baseData, replace
	}
	
end

capture: program drop buildFullData
program define buildFullData

	syntax [, REPLACE]

	buildBaseData, `replace'

	*	Unique geographic identifiers
	gen isForeign = Pa~="Bolivia"
	egen forcode = group(isForeign), missing
	drop isForeign
	egen ruecode = group(forc RUE), missing
	egen paicode = group(ruec idPa), missing
	egen depcode = group(paic idDep), missing
	egen procode = group(depc idPro), missing
	egen muncode = group(proc idMun), missing
	egen loccode = group(proc idLoc), missing
	egen discode = group(locc idDis), missing
	egen zoncode = group(disc idZon), missing
	egen reccode = group(zonc idRec), missing
	egen shecode = group(recc NúmeroMesa), missing

	bys recc: egen num_mesas = total(1)
	bys recc: egen num_early = total(reported==1)
	label def before_type -1 "Early" 0 "Split (Early)" 1 "Split (Late)" 2 "Late"
	gen before_type:before_type = cond(num_mesas==num_early, -1, cond(num_early==0, 2, reported<1))
	lab var before_type "Precinct split at announcement"
	drop num_mesas num_early

	egen sumV = rowtotal(CC-PANB)
	lab var sumV "Actual valid votes recorded"
	gen margin = 100*(MAS-CC)/sumV
	lab var margin "MAS lead over CC (% of valid votes)"
	gen r_margin = 100*(mS-mNO)/(mS+mNO)
	lab var r_mar "Sí lead over No (% of valid votes)"

end

capture: program drop doEscobariTables
program define doEscobariTables

	qui {
		gen isForeign = Pa~="Bolivia"
		egen fcode = group(isForeign)
		egen pcode = group(fcode idPai)
		egen dcode = group(pcode idDep)
		egen vcode = group(dcode idPro)
		egen mcode = group(vcode idMun)
		egen lcode = group(mcode idLoc)
		egen rcode = group(lcode idRec)
		egen scode = group(rcode NúmeroMesa)

		/*
		We exploit the fact that we have 34,555 polling stations divided into 455 municipalities, 
		which are then further divided into 3,593 localities and 5,288 precincts. Er.... ?

		Definitely failed to distinguish different localities with same name. Not sure how this compounded into precincts
		Not sure this wasn't corrected in regressions
		*/		
		
		foreach var of varlist CC-PANB {
			gen s_`var'= 100*`var'/VotosV
		}
		gen smargin = s_MAS-s_CC
		gen fmargin = MAS-CC
		gen SHUT = reported~=1
		gen NOTREP = reported==-1

		gen c_16 = 100*mNO/mVAL
		gen c_19 = s_CC
		gen m_16 = 100*mS/mVAL
		gen m_19 = s_MAS
		gen d_16 = 100*(mS-mNO)/mVAL
		gen d_19 = s_MAS-s_CC


		reshape long c_ m_ d_, i(shec) j(Year)

		gen Y2019 = Year==19
		keep if Y2019
	}
	
	qui {
		reg c_ 1.SHUT if Y2019, vce(robust)
		est store t2c1
		reg m_ 1.SHUT if Y2019, vce(robust)
		est store t2c2
		reg d_ 1.SHUT if Y2019, vce(robust)
		est store t2c3
		areg d_ 1.SHUT if Y2019, abs(mcod) vce(robust)
		est store t2c4
		tsset mcod scod
		xtreg  d_ 1.SHUT if Y2019, fe
		est store t2c4_f
		areg d_ 1.SHUT if Y2019, abs(lcod) vce(robust)
		est store t2c5
		tsset lcod scod
		xtreg  d_ 1.SHUT if Y2019, fe
		est store t2c5_f
		areg d_ 1.SHUT if Y2019, abs(rcod) vce(robust)
		est store t2c6
		tsset rcod scod 
		xtreg  d_ 1.SHUT if Y2019, fe
		est store t2c6_f
	}
	
	qui {
		noi di
		noi di as text "Regression results using areg"
		noi est table t2c?, stats(N r2) b(%5.3fc) se(%5.3fc)
		noi di
		noi di as text "F-tests (F_f) using xtreg, fe"
		noi est table t2c?_f, stats(F_f) b(%5.3fc)
	}


end

local ptype mlw(none) msize(vtiny)
local ptypeo `ptype' mcolor(ceprorange%\`mpct')
local ptypeg `ptype' mcolor(ceprgreen%\`mpct')

set seed 202110927
local newgraphs 0

buildFullData

qui {
log using ${OUTPUT}/newEscobari.txt, text replace

noi di as text `"-- "There are 20 coins in a bucket..." --"'
noi di as text "Pr(8+ nickels) = 1 in 1/(1-hypergeometric(20,10,10,7)) = 1 in " as result 1/(1-hypergeometric(20,10,10,7))
noi di

noi di as text `"-- "Preliminary Count Results..." --"'
egen TSE_sumV = rowtotal(TSE_CC-TSE_PANB)
lab var TSE_sumV "Announcement: Actual valid votes recorded"
gen TSE_absolute_margin = TSE_MAS-TSE_CC
lab var TSE_absolute "Announcement: MAS lead over CC"
gen TSE_margin = 100*TSE_absolute_margin/TSE_sumV
lab var TSE_margin "Announcement: MAS lead over CC (% of valid votes)"
count if reported==1
noi di as text "Percent of actas included: " as result 100*`r(N)'/`=_N'
noi table before_type if before_type<1, c(sum TSE_MAS sum TSE_CC sum TSE_absolute_margin sum TSE_sumV) row col format(%20.0fc)
noi table before_type if before_type<1 [aw=TSE_sumV], c(mean TSE_margin) row col format(%5.2fc)
noi di

egen rtag = tag(reccode)
noi table before_type if rtag

*	Margin needed post-announcement
scalar total_sumV = 6137671
sum TSE_sumV if before_type<1
scalar current_sumV = r(sum)
sum TSE_absolute_margin if before_type<1
scalar current_am = r(sum)
noi di as text "Margin needed in rest of count: " as result 100*(total_sumV/10-current_am)/(total_sumV-current_sumV)
noi di

noi di as text `"-- "Such a swing had precedent..." --"'
scalar si = 2028351
scalar no = 2363548
scalar sif = 2546135
scalar nof = 2682517
noi di as text "Percent of actas included in Opinión: " as result 100*24290/29224
noi di as text "2016 Margin included in Opinión: " as result  100*(si-no)/(si+no)
noi di as text "2016 Margin (final): " as result  100*(sif-nof)/(sif+nof)
noi di as text "2016 Margin (swing): " as result  100*((sif-nof)/(sif+nof)-(si-no)/(si+no))
noi di as text `"    [footnote]"'
noi di as text "2016 Margin (late): " as result 100*(sif-si-nof+no)/(sif-si+nof-no)
noi di as text "2016 Margin (shift): " as result 100*((sif-si-nof+no)/(sif-si+nof-no)-(si-no)/(si+no))
noi di as text "2019 Margin (swing) required: " as result 100*((total_sumV/10-current_am)/(total_sumV-current_sumV)-current_am/current_sumV)
noi di as text "2019 Margin (late, with half swing): " as result 100*((sif-si-nof+no)/(sif-si+nof-no)/2-(si-no)/(si+no)/2+current_am/current_sumV)
noi di as text "2019 Margin (final, with half swing): " as result 100*(((sif-si-nof+no)/(sif-si+nof-no)/2-(si-no)/(si+no)/2+current_am/current_sumV)*(total_sumV-current_sumV)+current_am)/total_sumV
noi di

noi di as text `"-- "Combining data from 2016..." --"'
noi tab isBig4 [w=Ins]
noi li reported Dist Zona idRec Rec MAS CC match if isBig4 & idDis==33
capture: confirm file ${FIGS}/Fig1.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	doBig4
	graph export ${FIGS}/Fig1.png, replace
}
noi di
noi di as text "	[Polling stations opposing referendum]"
noi tab reported if Loc=="Cochabamba" & r_margin<0
noi di
noi di as text "	[Polling stations favoring referendum]"
noi tab reported if Loc=="Cochabamba" & r_margin>0
capture: confirm file ${FIGS}/Fig2.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	doSplit
	graph export  ${FIGS}/Fig2.png, replace
}
noi di
noi di as text `"-- "... Were Not Representative..." --"'
noi table before_type, c(freq sum Ins) row col format(%12,0fc)
noi table before_type [aw=sumV], c(mean margin mean r_margin) row col format(%5.2fc)
noi table before_type if reported==1 [aw=sumV], c(mean margin mean r_margin) row col format(%5.2fc)
noi table before_type if reported~=1 [aw=sumV], c(mean margin mean r_margin) row col format(%5.2fc)
sum r_margin [aw=sumV] if before_type==1
scalar bt1 = r(mean)
sum r_margin [aw=sumV] if before_type==2
scalar dbt = r(mean)-bt1
noi di

noi di as text `"-- "Under pressure from the EOM..." --"'
noi table reported if inResum [w=sumV], c(mean r_margin) row col
preserve
getResumptionTSE
egen sumV = rowtotal(CC-PANB)
gen margin = 100*(MAS-CC)/sumV
noi mean margin [w=sumV]
restore
noi di
		
noi di as text `"-- "...Precinct Colegio Sebastían Pagador" --"'
noi table before_type if idPa==32 & idDep==3 & idLoc==821 & idRec==174 [aw=sumV], c(freq rawsum MAS rawsum CC rawsum sumV mean mar) row col cellw(12)
noi di
		
noi di as text `"-- "... Precincts Leaned Predictably..." --"'
noi di as text "	[Note that frequencies in tables below do not include four stations with zero valid votes]"
bys recc: egen before_sumV = total(cond(reported==1,sumV,.))
bys recc: egen before_margin = total(cond(reported==1,MAS-CC,.))
replace before_margin  = 100*before_margin/before_sumV
noi table before_type if reported==1 [aw=sumV], c(freq rawsum sumV mean before_margin mean margin) row col 
sum sumV
scalar tsumV = r(sum)
di tsumV
gen mard = cond(reported==1, MAS-CC, sumV*before_margin/100)
sum mard if before_type~=2
scalar cmar = r(sum)
sum sumV if before_type==2
scalar isumV = r(sum)
scalar marx = 100*(tsumV/10-cmar)/isumV
noi di as text "Margin needed in rest of count: " as result marx
sum before_margin [aw=sumV] if before_type==1
scalar xfill = r(mean)+dbt
gen xmargin = cond(reported==1, margin, cond(before_type==2, xfill, before_margin))
noi table before_type [aw=sumV], c(freq rawsum sumV mean xmargin mean margin) row col 
noi di

noi di as text `"-- "Escpbari and Hoover Require Massive Fraud..." --"'
sum MAS
scalar dvot = r(sum)
sum CC
scalar dvot = dvot-r(sum)
scalar ehd = dvot-5580*(13.77+15.57)
sum MAS if before_type~=2
scalar dvot2 = r(sum)
sum CC if before_type~=2
*	Note, only subtracting 4460, not 4462, because two relevant polling stations were annulled
scalar dvot2 = dvot2-r(sum)-4460*2.309
scalar marx3 = 100*(ehd-dvot2)/isumV
noi di as text "Margin needed in rest of count: " as result marx3
gen xmargin3 = cond(reported==1, margin, cond(before_type==2, marx3, 100*(MAS-CC-2.309)/sumV))
gen dxmar = xmargin3-r_margin
noi table before_type if reported==1 [aw=sumV], c(rawsum sumV mean r_margin mean xmargin3 mean margin mean dxmar) row col 
noi table before_type [aw=sumV], c(rawsum sumV mean r_margin mean xmargin3 mean margin mean dxmar) row col 
noi di

noi di as text `"-- "... transmitted every five minutes." --"'
gen ru = runiform()
sort UltTrans ru
gen ARRIVAL = (_n-1)/(_N-1)
gen t1000 = 1000*5*60*1000/(UltTrans-UltTrans[_n-1000])
capture: confirm file ${FIGS}/Fig3.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	line t1000 ARR, lw(thin) xti(" " "ARRIVAL of most recent {it:acta}") ///
		yti("Rate of last 1,000 arrivals per five minutes") ylab(0(500)1500)
	graph export  ${FIGS}/Fig3.png, replace
}
noi di

noi di as text `"-- "Hypothetical shift in support..." --"'
capture: confirm file ${FIGS}/Fig4.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway (function 10+50*x+(100-10-50*x)/5, col(ceprgreen) base(100)) (function 10+50*x, col(ceprorange)) ///
		(function 10+50*x+(100-10-50*x)/5, recast(area) lw(none) col(ceprgreen%50) base(100)) (function 10+50*x, recast(area) lw(none) col(ceprorange%50)) ///
		, yti("Share of Valid Votes (%)") xti(" " "ARRIVAL") legend(off) ylab(0(25)100) ///
		text(20 0.7 "Hard Support (MAS)") text(30 0.22 "Soft Support") text(70 0.4 "Opposition ({it:NO})") aspect(1)
	graph export  ${FIGS}/Fig4.png, replace
}
noi di

noi di as text `"-- "... the Direction of Fraud is Ambiguous..." --"'
capture: confirm file ${FIGS}/Fig5.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(function 10+50*x+(100-10-50*x)/5, recast(area) color(white) lc(ceprgreen) base(100)) ///
		(function 10+50*x, recast(area) lw(none) col(ceprorange%50)) ///
		(function 10+(1-1/5)*50*x, recast(area) col(white) lc(ceprorange)) ///
		, yti("Share of Valid Votes (%)") xti(" " "ARRIVAL") legend(off) ylab(0(25)100) ///
		text(15 0.65 "“Actual” Support (MAS)") text(70 0.4 "Opposition ({it:NO})") aspect(1) name(fraud1, replace)
	twoway ///
		(function 10+50*x+(100-10-50*x)/5, recast(area) color(white) lc(ceprgreen) base(100)) ///
		(function 10+(1-4/5)*50+(1-1/5)*50*x, recast(area) lw(none) col(ceprorange%50)) ///
		(function 10+50*x, recast(area) col(white) lc(ceprorange)) ///
		, yti("Share of Valid Votes (%)") xti(" " "ARRIVAL") legend(off) ylab(0(25)100) ///
		text(15 0.65 "Official Support (MAS)") text(70 0.4 "Opposition ({it:NO})") aspect(1) name(fraud2, replace)
	graph combine fraud1 fraud2, ycommon
	graph export  ${FIGS}/Fig5.png, replace
}
noi di

noi di as text `"-- "Actual Trends" --"'
gen x = .
gen mas = 100*MAS/sumV
gen sis = (r_margin+100)/2
gen isEarly = reported==1
capture: confirm file ${FIGS}/Fig6.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	plotlp mas sis ARR [aw=sumV], range(0 1) ///
		pre((rarea x x ARR, lw(none) color(ceprgreen%50)) (rarea x x ARR, lw(none) color(ceprorange%50)) (scatter mas ARR, msize(vtiny) mcolor(ceprorange%10)) ///
			(scatter sis ARR, msize(vtiny) mcolor(ceprgreen%10))) ///
		post(, legend(order(1 2) label(1 "2016 ({it:Sí})") label(2 "2019 (MAS)") col(1)) aspect(1) name(p1) yti("Vote Share") xti(" " "ARRIVAL"))
	graph export  ${FIGS}/Fig6.png, replace
}
capture: confirm file ${FIGS}/Fig7.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	plotlp mas sis [aw=sumV], by(isEarly) range(0 100) ///
		pre((rarea x x sis, lw(none) color(ceprgreen%50)) (rarea x x sis, lw(none) color(ceprorange%50)) (scatter mas sis if isEarly, msize(vtiny) mcolor(ceprorange%10)) ///
			(scatter mas sis if ~isEarly, msize(vtiny) mcolor(ceprgreen%10))) ///
		post(, legend(order(1 2) label(2 "Pre-SHUTDOWN") label(1 "SHUTDOWN") col(1)) aspect(1) name(p2) xti(" " "2016 {it:Sí} Vote Share") yti("2019 MAS Vote Share"))

	graph combine p1 p2, ycommon
	graph export ${FIGS}/Fig7.png, replace
}
noi di 
noi di as text `"-- "Fewest Ballots Report First" --"'
gen Emitidos = sumV+Bla+Nul
bys recc: egen nacta = total(1)
bys recc: egen amin = min(PriTrans)
gen isfirst = PriTrans==amin
bys recc (Emit): gen eorder = _n-1
bys recc (ARR): gen aorder = sum(ARR>ARR[_n-1])
bys recc: egen e0 = total(cond(eorder==0,Emit,.))
bys recc: egen e1 = total(cond(eorder==1,Emit,.))
gen de = e1-e0
replace de = e0 if nacta==1

egen ctmp = cut(nacta) if rtag & RUE~=2, g(10)
table ctmp if eorder==0, c(freq min nacta max nacta)
bys recc: egen cnac = mode(ctmp)
label var cnac "Number of Polling Stations in Precinct"
lab def cnac 3 "1 Station" 4 "2 Stations" 5 "3 Stations" 6 "4-5 Stations" 7 "6-10 Stations" 8 "12-16 Stations" 9 "17+ Stations", modify
lab val cnac cnac
table cnac if eorder==0, c(freq min nacta max nacta)
bys cnac: egen pmin = max(nacta)
bys cnac: egen pmax = min(nacta)
replace pmin = 100/pmin-0.0001
replace pmax = 100/pmax

probit isfirst c.de if eorder==0 & nacta>1 & RUE~=2
predictnl ihat = 100*predict(), ci(il ih)
noi ameans nacta if e(sample)
gen ph = 100/r(lb_h)
gen pl = 100/r(ub_h)
gen p = 100/nacta
sort de
capture: confirm file ${FIGS}/Fig8.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
			(rarea ph pl de, lw(none) color(ceprorange%50)) ///
			(rarea il ih de, pstyle(ci)) (line ihat de, pstyle(p1)) ///
		if e(sample), legend(order(3 1) label(3 "Observed") label(1 "At Random")) ///
			ylab(0(20)100) xti(" " "Next-Most Ballots Cast - Fewest Ballots Cast") yti("Percent") aspect(1)
	graph export ${FIGS}/Fig8.png, replace		
}
noi di
noi di as text `"-- and "By Number of Stations" --"'
probit isfirst c.de##ib(3).cnac if eorder==0 & nacta>1 & RUE~=2
predictnl ihat0 = 100*predict(), ci(il0 ih0)
sort de
capture: confirm file ${FIGS}/Fig9.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
			(rarea il0 ih0 de, pstyle(ci)) (line ihat0 de, pstyle(p1)) ///
			(line pmax de if cnac<6, color(ceprorange%20)) ///
			(rarea pmin pmax de, lw(none) color(ceprorange%20)) ///
		if e(sample), by(cnac, note("")) legend(order(1 3) label(1 "Observed") label(3 "At Random")) ///
			xti(" " "Next-Most Ballots Cast - Fewest Ballots Cast") yti("Percent")
	graph export ${FIGS}/Fig9.png, replace
}
gen isTrans = ~mi(UltTrans)
bys recc: egen hasmTrans = max(eorder==0 & mi(UltTrans))
bys recc: egen numTrans = total(~mi(UltTrans))
noi table isTrans if eorder==0 & RUE~=2, row col
noi table RUE if eorder==0 & RUE~=2 & mi(UltTrans), row col
noi table nacta if eorder==0 & RUE==0 & mi(UltTrans), row col
noi table Reci if RUE==0 & nacta==2 & hasmTrans, c(sum isTrans) row col
noi table Reci if RUE==0 & nacta>2 & hasmTrans, c(sum isTrans) row col
noi di

noi di as text `"-- Arrival Order of Polling Stations --"'
local mpct 50
capture: confirm file ${FIGS}/Fig10.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter x ARR, mcolor(ceprorange)) ///
		(scatter x ARR, mcolor(ceprgreen)) ///
		(scatter sis ARR if reported==1, `ptypeo') ///
		(scatter sis ARR if reported~=1, `ptypeg') ///
		, legend(order(1 2) label(1 "Pre-SHUTDOWN") label(2 "SHUTDOWN")) ///
		xti(" " "ARRIVAL") yti("2016 {it:Sí} Vote Shate") aspect(1)
	graph export ${FIGS}/Fig10.png, replace
}
noi di 
noi di as text `"-- Difference Model on 2016 Results --"'
gen SHUTDOWN = reported~=1
reg sis i1.SHUT
predict shat
capture: confirm file ${FIGS}/Fig17.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	twoway ///
		(scatter sis ARR if ~SHUT, `ptypeo') ///
		(scatter sis ARR if SHUT, `ptypeg') ///
		(line shat ARR if ~SHUT, lc(ceprorange)) ///
		(line shat ARR if SHUT, lc(ceprgreen)) ///
		, legend(col(1) order(3 4) label(3 "Pre-SHUTDOWN") label(4 "SHUTDOWN")) ///
		xti(" " "ARRIVAL")  yti("2016 {it:Sí} Vote Shate") aspect(1)
	graph export  ${FIGS}/Fig17.png, replace
}
noi di

noi di as text `"-- Analysis of Synthetic Data --"'
qui if ("`newgraphs'"=="newgraphs" | 1) {
	preserve
	include newEscobari_simulations
	restore
}

noi di
noi di as text `"-- Difference Model on 2016 Results --"'

egen mtag = tag(reported recc)
gen mts = 100*MTS/sumV
gen ccs = 100*CC/sumV
gen p2s = 100*p21/sumV
gen nos = 100-sis
noi correl MAS MTS CC p21F if reported==1

egen ncut = cut(nos) if reported==1 & mtag, g(5)
noi table ncut, c(freq min nos max nos)
egen ncut2 = cut(nos) if reported==1, at(0 19 32 46 63 100)
noi table ncut2 if mtag, c(freq min nos max nos)
gen pvalm = .
gen svalm = .
gen pvalc = .
gen svalc = .
levelsof ncut2, local(nl) clean
local lastlab 0
local labcom
foreach n of local nl {
	correl mas mts if reported==1 & mtag & ncut2==`n'
	replace pvalm = r(rho) if reported==1 & mtag & ncut2==`n'
	correl ccs p2s if reported==1 & mtag & ncut2==`n'
	replace pvalc = r(rho) if reported==1 & mtag & ncut2==`n'
	spearman mas mts if reported==1 & mtag & ncut2==`n'
	replace svalm = r(rho) if reported==1 & mtag & ncut2==`n'
	spearman ccs p2s if reported==1 & mtag & ncut2==`n'
	replace svalc = r(rho) if reported==1 & mtag & ncut2==`n'
	if (`n') {
		local labcom `labcom' `lastlab' "`lastlab'-`n'"
	}
	local lastlab `n'
}
lab def ncut2 `labcom' `lastlab' "`lastlab'-100"
lab val ncut2 ncut2
capture: confirm file ${FIGS}/Fig29.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	graph bar pvalm pvalc if reported==1 & mtag, bar(1, color(ceprblue%50)) bar(2, color(ceprorange%50)) over(ncut2) b1t("NO Vote Share in 2016") ///
		legend(label(1 "MAS vs MTS") label(2 "CC vs 21F")) l1t("Correlation in Vote Shares")
	graph export ${FIGS}/Fig29.png, replace
}
noi di

noi di `"-- "Minor Parties..." --"'
gen mps = mts+p2s
gen mpm = mts-p2s
gen diff = margin-2*mps
noi table before_type [w=sumV], c(mean margin mean mps mean mpm mean diff) row col

noi di
noi di "*** Begin Appendices ***"
noi di

noi di as text "Appendix A"
noi di as text `"-- "Regression discontinuity" --"'

bys reccode: egen vsum = total(VotosV), missing
bys reccode: egen rtmar = total(100*(MAS-CC)/vsum), missing

gen tmar = 100*(MAS-CC)/VotosV
gen smar = 100*(mS-mNO)/mVAL
rdrobust tmar ARR if ~mi(UltTrans), p(2) c(0.95)
capture: confirm file ${FIGS}/FigA1.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	rdplot tmar ARR if ~mi(UltTrans) & ARR>=0.75, p(2) c(0.95) graph_options(ysc(r(-20 80)) legend(off) ti("") xti(" " "Arrival of the Polling Stations at the Electoral Court") yti("MAS-CC"))
	graph export ${FIGS}/FigA1.png, replace
}
capture: confirm file ${FIGS}/FigA2.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	rdplot smar ARR if ~mi(UltTrans) & ARR>=0.75, p(2) c(0.95) graph_options(ysc(r(-20 80)) legend(off) ti("") xti(" " "Arrival of the Polling Stations at the Electoral Court") yti("{it:SÍ-NO}"))
	graph export ${FIGS}/FigA2.png, replace
}

noi di
noi di as text "Appendix D"
noi di as text `"-- "Replication of Escobari and Hoover" --"'
noi doEscobariTables

noi di
noi di as text "Appendix F"
noi di as text `"-- "Monte Carlo Imputation of Late Results" --"'
tempfile tf
getAnnouncementTSE
egen other = rowtotal(CC-PANB)
gen nonv = Ins-other
replace other = other-(MAS+CC)
keep NúmeroMesa MAS CC other nonv
compress
save `tf'
buildFullData
keep NúmeroMesa forc-shec Ins
compress
merge 1:1 NúmeroMesa using `tf', nogen
doMC MAS CC other nonv [w=Ins], ids(forcode-shecode) dir(yTSE) i(1000) keep(recc)
compressIters, dir(yTSE)
capture: confirm file ${FIGS}/FigF.png
if (_rc~=0 | "`newgraphs'"=="newgraphs") {
	graph export ${FIGS}/FigF.png, replace
}
noi codebook MCm

log close
