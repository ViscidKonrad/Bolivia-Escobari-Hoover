
capture: program drop multiDraw
program define multiDraw

	syntax varlist [if] [in], Prefix(string asis)
	
	marksample touse
	
	tempvar nn rr
	gen `nn' = `: word 1 of `varlist'' if `touse'
	gen `rr' = 1 if `touse'
	local nv: list sizeof varlist
	forvalues ii=2/`=`nv'-1' {
		local var: word `ii' of `varlist'
		macro li _var
		gen `prefix'`var' = cond(`nn'==0,0,cond(`var'/`rr'<1e-8,0,cond(`var'/`rr'>1-1e-8,`nn',rbinomial(`nn',`var'/`rr')))) if `touse' & `rr'>0
		replace `nn' = `nn'-`prefix'`var' if `touse' & `rr'>0
		replace `rr' = `rr'-`var' if `touse' & `rr'>0
	}
	gen `prefix'`: word `nv' of `varlist'' = `nn'

end

capture: program drop dirichletDraw
program define dirichletDraw

	syntax varlist [if] [in], Prefix(string asis)
	
	marksample touse
	tempvar tgamma

	*	Dirichlet draw
	local dlist
	foreach var of local varlist {
		gen `prefix'`var' = cond(`var'<1e-4,0,rgamma(`var',1)) if `touse'
		local dlist `dlist' `prefix'`var'
	}
	egen `tgamma' = rowtotal(`dlist')
	foreach var of local dlist {
		replace `var' = `var'/`tgamma'
	}

end

capture: program drop makeLists
program define makeLists

	syntax [anything(name=numlevel)] [if] [, IDs(varlist) Code(integer 1)]
	
	marksample touse
	
	local lmax: list sizeof ids

	qui {
		if ("`numlevel'"=="") {
			local numlevel 1
		}
		if (`numlevel'==1) {
			foreach l of local ids {
				gen `l'_list = ""
				gen `l'_len = 0
			}
		}
		if (`numlevel'<6) {
			local pref noi
		}

		local id: word `numlevel' of `ids'	
		`pref' levelsof `id' if `touse', local(lcl)
		replace `id'_lis = "`lcl'" in `code'
		replace `id'_len = `: list sizeof lcl' in `code'		

		if (`numlevel'<`lmax') {
			foreach lc of local lcl {
				noi makeLists `=`numlevel'+1' if `touse' & `id'==`lc', id(`ids') code(`lc')
			}
		}
	}

end

capture: program drop doDataPrep
program define doDataPrep

	syntax varlist [if] [in]
	
	marksample touse
		
	tokenize `varlist'
	local plevel `1'
	macro shift
	local rlevels `*'
	foreach level of local rlevels {
		gen MC`plevel' = `plevel' if `level'_len[`plevel']~=0
		local plevel `level'
	}
	gen MC`plevel' = `plevel' if `touse'

end

capture: program drop doImpute
program define doImpute

	*
	*	If strict is specified, imputations will be made strict within the smallest available geography
	*	That is, the default is to have some probability of sourcing imputations from outside that geography
	*		i.e., if we have one sheet in a precinct of 12, we may impute the other 11 from a different precinct
	*
	
	syntax varlist, STRICT(varlist)
		
	tokenize `varlist'
	local level `1'
	macro shift
	local rlevels `*'

	replace MC`level' = . if `level'_tag & ~`: list level in strict' & ln(runiform())<`level'_nsheets*ln1m(1/`level'_len[1])
	replace MC`level' = real(word(`level'_list[1],runiformint(1,`level'_len[1]))) if `level'_tag & mi(MC`level')
	replace MC`level' = MC`level'[`level'_rep] // if mi(MC`level')
	local plevel `level'
	foreach level of local rlevels {
		replace MC`level' = . if `level'_tag & (MC`plevel'~=`plevel' | (~`: list level in strict' & ln(runiform())<`level'_nsheets*ln1m(1/`level'_len[MC`plevel'])))
		replace MC`level' = real(word(`level'_list[MC`plevel'], runiformint(1,`level'_len[MC`plevel']))) if `level'_tag & mi(MC`level')
		replace MC`level' = MC`level'[`level'_rep] // if mi(MC`level')
		local plevel `level'
	}
	replace MC`plevel' = `plevel' if hasData
	
end

capture: program drop doMC
program define doMC

	syntax varlist [if/] [in/] [aweight/] [, STRICT(passthru) REPREP REPLACE DETAIL IDs(varlist) Directory(string asis) Keep(varlist) Iters(integer 50)]

	marksample touse
	
	*	number of levels
	local nls: list sizeof ids
	*	observation id code
	local obsid: word `nls' of `ids'
	*	precinct-level id code
	local groupid: word `=`nls'-1' of `ids'
	*	level ids other than observation
	local idsml: list ids - obsid
	
	if ("`directory'"=="") {
		local directory test
	}
	if ("`reprep'"=="reprep") {
		local replace replace
	}
	if ("`strict'"=="") {
		local strict strict(`groupid')
	}
	if ("`exp'"=="") {
		gen weight = 1
	}
	else {
		gen weight = `exp'
	}
		
	capture: confirm file ${MC}/`directory'/prepped.dta
	if (_rc~=0 | "`reprep'"=="reprep") {

		capture: mkdir ${MC}
		capture: mkdir ${MC}/`directory'

		gen hasData = `touse'
		keep `ids' hasData weight `varlist' `keep'
		compress
		save ${MC}/`directory'/basedata, replace

		collapse (sum) hasData `varlist', by(`idsml')
		foreach var of local varlist {
			gen alpha`var' = cond(hasData,1+`var',.)
		}
	
		*	Fix data order to allow table lookup
		sort `groupid'
		gen impute_order = _n
		sort impute_order

		*	Tag one of each geography as a representative
		foreach id of local idsml {
			*	Tag representative
			egen `id'_tag = tag(`id')
			*	Create pointer to representative
			bys `id': egen `id'_rep = total(`id'_tag*impute_order)
			*	Count number of observations in geography
			bys `id': egen `id'_nsheets = total(hasData)
		}
		
		makeLists if hasData, id(`idsml')
		doDataPrep `idsml' if hasData
		compress
		save ${MC}/`directory'/prepped, replace
		
	}
	else {
		use ${MC}/`directory'/prepped.dta, clear
	}
	
	preserve
	local mcvars
	local pvars
	foreach var of local varlist {
		local mcvars `mcvars' MC`var'
		local pvars `pvars' p_MC`var'
	}
	
	qui forvalues iter = 1/`iters' {
	
		capture: confirm file ${MC}/`directory'/iterations/`iter'.dta
		if (_rc~=0 | "`replace'"=="replace") {
	
			capture: mkdir ${MC}/`directory'/iterations
			restore, preserve
			noi di as text "Iteration: " as result `iter'
			doImpute `idsml', `strict'

			foreach var of local varlist {
				gen MC`var' = alpha`var'[MC`groupid']
			}
			*	Draw probability vectors
			dirichletDraw `mcvars', prefix(p_)
			keep `groupid' `pvars'
			merge 1:n `groupid' using ${MC}/`directory'/basedata, nogen
			
			*	Draw vote vectors
			multiDraw weight p_* if ~hasData, prefix(n_)
			ren n_p_* *
			foreach var of local varlist {
				replace MC`var' = `var' if hasData
			}
			keep hasData `obsid' MC* `keep'
			gen iter = `iter'
			order iter `obsid' hasData `keep'
			compress
			save ${MC}/`directory'/iterations/`iter', replace
		}
	}
	restore, not

end

capture: program drop compressIters
program define compressIters

	syntax [, Directory(string asis)]

	if ("`directory'"=="") {
		local directory test
	}

	tempfile tf
	clear
	save `tf', emptyok
	local ilist: dir "${MC}/`directory'/iterations/" files "*"

	local tmesas
	local nmesas
	qui foreach i of local ilist {
		use ${MC}/`directory'/iterations/`i', clear
		if ("`tmesas'"=="") {
			count 
			local tmesas = r(N)
			count if hasData
			local nmesas = r(N)
		}
		collapse (sum) MCCC MCMAS MCother, by(iter)
		append using `tf'
		save `tf', replace
	}
	gen MCmargin = 100*(MCMAS-MCCC)/(MCM+MCC+MCo)
	twoway kdensity MCmargin, recast(area) lw(none) color(ceprblue%25) ///
		xti(" " "Final MAS margin (simulated, percent of valid vote)") ysc(off) ///
		caption("`:di %5.0fc _N' iterations based on results from `:di %6.0fc `nmesas'' of `:di %6.0fc `tmesas'' mesas")


end

