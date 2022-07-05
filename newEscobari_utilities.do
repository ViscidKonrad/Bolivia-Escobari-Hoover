

*
*	This document contains various utilities, particularly regarding:
*	 - the retrieval of raw data
*	 - the most basic processing of files into usable format
*

*
*	PROBLEM WITH -IMPORT EXCEL-
*		any timestamp data imported as strings comes in with a 12 hour clock and no AM/PM indicator
*		any timestamp data imported as numeric are fine, but not necessarily formatted usefully
*	WHY THIS IS A PROBLEM
*		information may be lost if -IMPORT EXCEL- imports as string, say if misisng values are listed as "null"
*	SOLUTION
*		reach into the underlying XML files and pull data directly as string
*
*	CONSEQUENCE
*		this fix is probably NOT something Nooruddin did for the OAS, so some timestamps may differ
*

*
*	GENERIC UTILITIES
*
*	Generic Data Utilities
*
*	Simply unzips an XLSX file to the specified location
capture: program drop unpackXLSX
program define unpackXLSX

	syntax anything(name=xlsxfile), [To(string asis)] [REPLACE]
	
	*
	*	Note that XLSX files are just other stuff packaged in a zip
	*
		
	capture: confirm file `to'
	if (_rc~=0 | "`replace'"=="replace") {
		! unzip `xlsxfile' -d `to'
	}

end

*	Read in an Excel worksheet (from XML)
capture: program drop getWorksheet
program define getWorksheet

	syntax [anything(name=sheetname)], From(string asis)

	*
	*
	
	if ("`sheetname'"=="") {
		local sheetname sheet1
	}
	
	*
	*	First, get any shared strings
	*
	tempfile ss
	*	Add newlines after each cell for easier parsing
	filefilter `from'/xl/sharedStrings.xml `ss', from("</si>") to("</si>\U") replace
	*	Read in the data
	import delimited using `ss', delim("</si>\n", asstring) enc("utf-8") clear stringcols(_all)
	gen sharedString = regexs(1) if regexm(v1,"<si><t[^>]*>([^<]*)")
	drop v1
	drop if mi(sharedString)
	gen stringCode = _n-1
	compress
	save `ss', replace

	*
	*	Read in the sheet XML
	*
	tempfile fft xml
	*	Add newlines after each cell for easier parsing
	filefilter `from'/xl/worksheets/`sheetname'.xml `fft', from("\Q/><c r=\Q") to("\Q/>\U<c r=\Q") replace
	filefilter `fft' `xml', from("</c>") to("</c>\U") replace
	*	Read in the data
	import delimited using `xml', delim("</c>\n", asstring) enc("utf-8") clear stringcols(_all)

	*	Parse out the basics
	gen cell = regexs(1) if regexm(v1,`"c r="([^"]*)""')
	gen row = regexs(1) if regexm(cell,"[A-Z]+([0-9]+)")
	gen col = regexs(1) if regexm(cell,"([A-Z]+)[0-9]+")
	keep if ~mi(row) & ~mi(col)
	gen value = regexs(1) if regexm(v1,"<v>([^<]*)</v>")
	gen isSString = regexm(v1,`"t="s""')
	gen isIString = regexm(v1,`"t="inlineStr""')

	*	Handle shared strings
	gen stringCode = value if isSString
	destring stringCode, force replace
	compress
	merge n:1 stringCode using `ss', keep(match master)
	li if ~mi(stringCode) & _merge==1
	replace sharedString = "p21F" if sharedString=="21F"
	replace sharedString = "PANBOL" if sharedString=="PAN-BOL"
	replace value = sharedString if isSString

	*	Handle inline strings
	gen inlineString = regexs(1) if regexm(v1,"<is><t>([^<]*)</t></is?")
	replace inlineString = "p21F" if inlineString=="21F"
	replace inlineString = "PANBOL" if inlineString=="PAN-BOL"
	replace value = inlineString if isIString

	*	Reshape into rows
	drop v1 cell _merge
	drop isSString-inlineString
	destring row, force replace
	replace col = strlower(col) if length(col)>1
	reshape wide value, i(row) j(col) string
	drop row

	*
	*	Cleanup
	*
	*	Set variable names from first row
	foreach var of varlist value* {
		capture: ren `var' `=substr(subinstr(`var'[1]," ","",.),1,32)'
	}
	drop in 1
	compress

end

*	End Generic Data Utilities

*
*	Generic Analytic Utilities
*

*	Convert lat/long to x/y in meters
capture: program drop LambertConformalEllipsoid
program define LambertConformalEllipsoid

	syntax varlist(min=2 max=2) [if] [in], [A(real 6378137) Finv(real 298.257223563)] ///
		Standards(numlist min=2 max=2 sort) Origin(numlist min=2 max=2) ///
		[GENerate(namelist min=2 max=2)] [REPLACE] [FALSEEasting(real 0)]  [FALSENorthing(real 0)]

	marksample touse

	*
	*	Ellipsoidal formulas for Lambert Conformal Conic (WGS84)
	*		Lat/Long in degrees -> x/y in meters
	*
	
	*	https://pubs.usgs.gov/pp/1395/report.pdf
	*	starts bottom page 107
	
	tempname p4 pc ecc2 ecc phi1 phi2 phi0 lam0 m1 m2 t1 t2 t0 n Fish rho0
	tempvar phi lam tp rho theta
		
	*	Conversion from degrees to radians
	scalar `p4' = _pi/4
	scalar `pc' = _pi/180
	
	*	Conversion from inverse flattening to eccentricity
	scalar `ecc2' = (2-1/`finv')/`finv'
	scalar `ecc'= sqrt(`ecc2')
		
	*	Standard Parallels
	tokenize `standards'
	scalar `phi1' = (`1')*`pc'
	scalar `phi2' = (`2')*`pc'
	
	*	Origin
	tokenize `origin'
	scalar `phi0' = (`1')*`pc'
	scalar `lam0' = (`2')*`pc'
		
	*	Constants
	forvalues i=1/2 {
		scalar `m`i'' = cos(`phi`i'')/sqrt(1-`ecc2'*(sin(`phi`i''))^2)
	}
	forvalues i=0/2 {
		scalar `t`i'' = tan(`p4'-`phi`i''/2)/exp(`ecc'*ln((1-`ecc'*sin(`phi`i''))/(1+`ecc'*sin(`phi`i'')))/2)
	}
	scalar `n' = ln(`m1'/`m2')/ln(`t1'/`t2')
	scalar `Fish' = `a'*`m1'/`n'
	scalar `rho0' = `Fish'*exp(`n'*ln(`t0'/`t1'))
	
	*	Conversion to radians
	tokenize `varlist'
	gen `phi' = (`1')*`pc' if `touse'
	gen `lam' = (`2')*`pc' if `touse'
	
	*	Point-specific parameters
	gen `tp' = tan(`p4'-`phi'/2)/exp(`ecc'*ln((1-`ecc'*sin(`phi'))/(1+`ecc'*sin(`phi')))/2) if `touse'
	gen `rho' = `Fish'*exp(`n'*ln(`tp'/`t1')) if `touse'
	gen `theta' = `n'*(`lam'-`lam0') if `touse'
	
	*	Coordinates
	tokenize `generate'
	gen `1' = `rho'*sin(`theta')+`falseeasting' if `touse'
	gen `2' = `rho0'-`rho'*cos(`theta')+`falsenorthing' if `touse'
		
end

*	End Generic Analytic Utilities

*
*	DATA SPECIFIC UTILITIES
*
*	Specific Retrieval Utilities
*
*	Fetch TSE data and un-rar <-- NOTE MOSTLY FOR REFERENCE
capture: program drop fetchTSE
program define fetchTSE

	syntax [, REPLACE]
	
	local directory ${DATA}/TSE
	local filename `directory'/datos_trep_computo_eg2019
	capture: confirm file `filename'
	if (_rc~=0 | "`replace'"=="replace") {
		capture: mkdir `directory'
		local URL https://github.com/ViscidKonrad/Bolivia-Elections-2019/raw/master/datos_trep_computo_eg2019.rar
		copy `"`URL'"' `filename'.rar, `replace'
	}
	/*
	 *	The resulting archive must be unpacked and re-unpacked to produce the following files in particular:
	 *	${DATA}/TSE/datos_trep_computo_eg2019/computo/2.RecepcionSobres_final.xlsx
	 *	${DATA}/TSE/datos_trep_computo_eg2019/trep/OEA_TREP_LogCompleto_19G_2019.11.03.06.53.xlsx 
	 *
		! brew install rar
		! /usr/local/bin/unrar e `filename'.rar `filename'/
		! /usr/local/bin/unrar e `filename'/computo.rar `filename'/computo/
		! /usr/local/bin/unrar e `filename'/trep.rar `filename'/trep/
	 */

end	

*	Fetch OEP public reports <- NOTE MOSTLY FOR REFERENCE
capture: program drop getPublicReports
program define getPublicReports

	syntax, [TREP] [COMPUTO] [REPLACE]
	
	local count `trep'`computo'
	if ("`count'"=="trepcomputo") {
		local count
	}
	if ("`count'"=="") {
		di as error "Must choose one of TREP or COMPUTO"
		local count computo
	}
	capture: confirm file ${DATA}/OEP/reports_`count'.zip
	if (_rc~=0 | "`replace'"=="replace") {
		capture: mkdir ${DATA}/OEP
		local URL https://archive.org/download/election_bolivia_2019/reports_`count'.zip
		copy `"`URL'"' ${DATA}/OEP/reports_`count'.zip, `replace'
	}
	 *	unzip to ${DATA}/OEP/reports_`count'
	
end

*	Fetch Padrón Electoral for 2016 or 2019, internal or exterior
capture: program drop getPadron
program define getPadron

	syntax anything(name=year) [, REPLACE EXterior]
	
	if (`year'==2016) {
		if ("`exterior'"=="exterior") {
			local dtafile ${INPUT}/2016/OEP/padron16x
			capture: use `dtafile', clear
			if (_rc~=0 | "`replace'"=="replace") {
				capture: mkdir ${INPUT}/2016
				capture: mkdir ${INPUT}/2016/OEP
				local xlsxfile ${DATA}/2016/OEP/padron_ref_const_2016.xlsx
				capture: confirm file `xlsxfile'
				if (_rc~=0) {
					capture: mkdir ${DATA}/2016
					capture: mkdir ${DATA}/2016/OEP
					local URL https://www.oep.org.bo/wp-content/uploads/2017/02/padron_ref_const_2016.xlsx
					copy `"`URL'"' `xlsxfile', replace
				}
				import excel using `xlsxfile', sheet(Ciu_hab_por mesa INTER) cellrange(a5) firstrow case(preserve) clear				
				ds *, has(type string)
				foreach var of varlist `r(varlist)' {
					replace `var' = ustrtrim(`var')
				}
				compress
				save `dtafile', replace
			}
		}
		else {
			local dtafile ${INPUT}/2016/OEP/padron16
			capture: use `dtafile', clear
			if (_rc~=0 | "`replace'"=="replace") {
				capture: mkdir ${INPUT}/2016
				capture: mkdir ${INPUT}/2016/OEP
				local xlsxfile ${DATA}/2016/OEP/padron_ref_const_2016.xlsx
				capture: confirm file `xlsxfile'
				if (_rc~=0) {
					local URL https://www.oep.org.bo/wp-content/uploads/2017/02/padron_ref_const_2016.xlsx
					copy `"`URL'"' `xlsxfile', replace
				}
				import excel using `xlsxfile', sheet(Ciu_hab_por mesa NACIONAL) cellrange(a5) firstrow case(preserve) clear
				ds *, has(type string)
				foreach var of varlist `r(varlist)' {
					replace `var' = trim(`var')
				}
				replace NombreRec = `"Unidad Educativa "Simón Bolivar La Palca""' if NombreRec==`"Unidad Educativa "Simón Bolivar La Palca"'
				compress
				save `dtafile', replace
			}
		}
	}
	if (`year'==2019) {
		if ("`exterior'"=="exterior") {
			local dtafile ${INPUT}/OEP/padron19x
			capture: use `dtafile', clear
			if (_rc~=0 | "`replace'"=="replace") {
				capture: mkdir ${INPUT}/OEP
				local xlsxfile ${DATA}/OEP/Estadisticas_Exterior_EG_2019.xlsx
				capture: confirm file `xlsxfile'
				if (_rc~=0) {
					capture: mkdir ${DATA}/OEP
					local URL https://www.oep.org.bo/wp-content/uploads/2019/09/Estadisticas_Exterior_EG_2019.xlsx
					copy `"`URL'"' `xlsxfile', replace
				}
				import excel using `xlsxfile', sheet(Hab_inhab_depu_X_Recinto) cellrange(a7) firstrow case(preserve) clear
				keep if ~mi(pais)
				destring dep, force replace
				ds *, has(type string)
				foreach var of varlist `r(varlist)' {
					replace `var' = trim(`var')
				}
				*	Error: Granada (Localidad) should be 4393, not 4392
				replace idloc=4393 if pais==67 & dep==30 & PROV==8 & sec==1 & idloc==4392
				compress
				save `dtafile', replace	
			}
		}
		else {
			local dtafile ${INPUT}/OEP/padron19
			capture: use `dtafile', clear
			if (_rc~=0 | "`replace'"=="replace") {
				capture: mkdir ${INPUT}/OEP
				local xlsxfile ${DATA}/OEP/Estadisticas_Nacional_EG_2019.xlsx
				capture: confirm file `xlsxfile'
				if (_rc~=0) {
					capture: mkdir ${DATA}/OEP
					local URL https://www.oep.org.bo/wp-content/uploads/2019/09/Estadisticas_Nacional_EG_2019.xlsx
					copy `"`URL'"' `xlsxfile', replace
				}
				import excel using `xlsxfile', sheet(Hab_Inhab_depu_X_Mesa) cellrange(a5) firstrow case(preserve) clear
				keep if ~mi(dep)
				ds *, has(type string)
				foreach var of varlist `r(varlist)' {
					replace `var' = trim(`var')
				}
				compress
				save `dtafile', replace
			}
		}
	}

end

*	Fetch historical election data
capture: program drop getDatosAbiertos
program define getDatosAbiertos

	syntax anything(name=year) [, EXterior FETCH]
	
	local suffix
	local slist 9
	if ("`exterior'"=="exterior") {
		local suffix _exterior
		local slist 3 5 6
	}

	local csvfile ${DATA}/`year'/OEP/votos_totales`suffix'.csv
	capture: confirm file `csvfile'
	if (_rc~=0 | "`fetch'"=="fetch") {
		
		capture: mkdir ${DATA}/`year'
		capture: mkdir ${DATA}/`year'/OEP

		local index
		if (`year'==2014) {
			local index = cond("`exterior'"=="exterior", 85, 17)
		}
		if (`year'==2016) {
			local index = cond("`exterior'"=="exterior", 83, 38)
		}

		tempfile tf tf2
		copy `"https://atlaselectoral.oep.org.bo/descarga/`index'/votos_totales.csv"' `tf', text
		filefilter `tf' `tf2', replace from("LA PALCA") to("LA PALCA\Q")
		copy `tf2' `csvfile', replace
	
	}
	import delimited using `csvfile', delim("|") clear varnames(1) case(preserve) enc("utf8") stringcols(`slist')
	
end

*	Fetch a single data file from the OEP Geoelectoral site
*		noting if it is the last file in the data set or not
capture: program drop getSingleGeo
program define getSingleGeo

	syntax anything(name=server) [, Num(integer 0) REPLACE FETCH]
	
	if ("`fetch'"=="fetch") {
		local replace replace
	}
	
	local dtafile ${INPUT}/Geo/`server'/`num'
	capture: use `dtafile', clear
	if (_rc~=0 | "`replace'"=="replace") {
		capture: mkdir ${INPUT}
		capture: mkdir ${INPUT}/Geo
		capture: mkdir ${INPUT}/Geo/`server'
		local jfile ${DATA}/Geo/`server'/`num'.json
		capture: confirm file `jfile'
		if (_rc~=0 | "`fetch'"=="fetch") {
			local URLHEAD https://geoelectoral.oep.org.bo/oep/rest/services/
			local URLTAIL &outFields=*&returnGeometry=true&f=json
			if ("`server'"=="RE16") {
				local URLSERV Recintos_Electorales/Recintos_Electorales_18_04_2017/MapServer/0/query?where=(ObjectID>=`num')
				*GeografiaElectoral/Asientos_Recintos_19_09_17/MapServer/1
			}
			if ("`server'"=="UE16") {
				local URLSERV UnidadesEducativas/UnidadesEducativas_2016/MapServer/0/query?where=(ObjectID>=`num')
			}
			if ("`server'"=="EN19") {
				*	Resultados 20/10/19
				local URLSERV EleccionesGrales2019/Elecciones_Nacionales_2019/MapServer/2/query?where=(ObjectID>=`num')
			}
			if ("`server'"=="ENX19") {
				*	Resultados 20/10/19
				local URLSERV EleccionesGrales2019/Elecciones_Nacionales_2019_Exterior/MapServer/1/query?where=(FID>=`num')
			}
			if ("`server'"=="UE") {
				local URLSERV UnidadesEducativas/UnidadesEducativas_2016/MapServer/0/query?where=(OBJECTID>=`num')
			}
			capture: mkdir ${DATA}/Geo/
			capture: mkdir ${DATA}/Geo/`server'
			tempfile rf
			
			copy `"`URLHEAD'`URLSERV'`URLTAIL'"' `rf'
			filefilter `rf' `jfile', replace from("\LQ") to("'")
			
		}
		clear
		
		insheetjson using `jfile', showresponse flatten

		*	Get list of fields
		tempname jf
		file open `jf' using `jfile', read text
		file read `jf' endline
		local alist
		if regexm(`"`endline'"',`""fieldAliases":\{([^\}]*)\},"') {
			local fa = regexs(1)
			macro li _fa
			tokenize `"`fa'"', parse(":,")
			while (`"`1'"'~="") {
				local alist `alist' `1'
				macro shift 4
			}
			macro li _alist
		}
		*	Find out if there is more data after this
		file seek `jf' eof
		file seek `jf' `=r(loc)-29'
		file read `jf' endline
		file close `jf'
		*	Produce empty variables
		foreach a of local alist {
			local cl `cl' "attributes:`a'"
			capture: gen str256 `a'=""
			if (_rc~=0) {
				gen str256 x`a'=""
			}
		}
		gen str256 x=""
		gen str256 y=""
		local cl `cl' "geometry:x" "geometry:y"
		*	Read in JSON
		insheetjson * using `jfile', tableselector(features) col(`cl')
		*	Add note indicating if there is more data to be read
		if (regexm(`"`endline'"',`""exceededTransferLimit":true"')) {
			notes: eft_true
		}
		else {
			notes: eft_false
		}
		save `dtafile', replace
	}

end

*	Fetch a full data set from the OEP Geoelectoral site
capture: program drop getAllGeo
program define getAllGeo

	syntax anything(name=server) [, REPLACE FETCH]

	local num 0
	if (inlist("`server'","EN19","UE","RE16")) {
		local num 1
	}
	if ("`fetch'"=="fetch") {
		local replace replace
	}
	
	local dtafile ${INPUT}/Geo/`server'
	capture: use `dtafile', clear
	if (_rc~=0 | "`replace'"=="replace") {
	
		tempfile tf
		clear
		save `tf', emptyok
		local eft_true 1
		while (`eft_true') {
			getSingleGeo `server', num(`num') `replace' `fetch'
			if ("`:char _dta[note1]'"=="eft_false") {
				local eft_true 0
			}
			append using `tf'
			save `tf', replace
			local num = `num'+1000
		}
		foreach var of varlist * {
			replace `var' = subinstr(`var',"\u00A0","",.)
			replace `var' = subinstr(`var',"\","",.)
		}
		capture: replace As = `"Rodeo "A""' if As=="Rodeo A"
		capture: replace NombreRec = `"U.E. Ramon Dario Gutierrez "A""' if NombreRec=="U.E. Ramon Dario Gutierrez A"
		capture: replace NombreRec = `"Public School 69 Jackson Heights"' if NombreRec=="Public School 69"
		save `dtafile', replace
		
	}

end

*	End Specific Retrieval Utilities

*
*	Data Construction Utilities
*

*	Read in Chumacero's assignments
capture: program drop readChumacero
program define readChumacero

	*
	*	NOTE: Data is not available online
	*

	syntax [, REPLACE]

	capture: use ${INPUT}/DTA/urb_rural, clear
	if (_rc~=0 | "`replace'"=="replace") {
		capture: import excel using ${DATA}/urb_rural.xlsx, firstrow case(preserve) clear
		if (_rc~=0) {
			di as error "Must place " as result "urb_rural.xlsx" as error " in " as result "${DATA}/urb_rural/"
			exit(_rc)
		}
		else {
			compress
			ren Rur RUE
			lab var RUE "Chumacero's assignments of domestic urban/rural localities"
			lab def RUE 0 "Rural" 1 "Urban" 2 "External"
			lab val RUE RUE
			save ${INPUT}/DTA/urb_rural, replace
		}
	}
	
end

*	Read in election update
capture: program drop getUpdate
program define getUpdate

	syntax anything(name=timestamp) [, TREP COMPUTO REPLACE]
	
	local counts `trep'`computo'
	if ("`counts'"=="trepcomputo") {
		di as error "May choose at most one of TREP or COMPUTO"
		exit(99)
	}
	if ("`counts'"=="") {
		local counts trep computo
	}
	
	local found 0
	foreach count of local counts {
		if (`found'==0) {
			local dtafile OEP/reports_`count'/`timestamp'
			capture: use ${INPUT}/`dtafile', clear
			if (_rc~=0) {
				local xlsxfile = subinstr("`dtafile'","`timestamp'","acta."+subinstr("`timestamp'","_",".",.),.)
				capture: confirm file ${DATA}/`xlsxfile'.xlsx
				if (_rc==0) {
					capture: import excel using ${DATA}/`xlsxfile'.xlsx, firstrow case(preserve) clear
					if (_rc==0 & `=_N') {
						local found 1
						capture: mkdir ${INPUT}/OEP
						capture: mkdir ${INPUT}/OEP/reports_`count'
						keep if ustrregexm(Ele,"Pres")
						compress
						save ${INPUT}/`dtafile', replace
					}
				}
			}
		}
	}
	
end

*	Wrappers for particular updates
*	Wrapper for results in TSE Announcement
capture: program drop getAnnouncementTSE
program getAnnouncementTSE

	getUpdate 2019_10_20_19_40_57, trep

end

*	Wrapper for results immediately upon resumption of TREP
capture: program drop getResumptionTSE
program getResumptionTSE

	getUpdate 2019_10_21_18_29_53, trep

end

*	Wrapper for final TREP
capture: program drop getFinalTREP
program getFinalTREP

	getUpdate 2019_10_25_21_07_40, trep

end

*	Wrapper for final computo
capture: program drop getFinalComputo
program define getFinalComputo

	getUpdate 2019_10_25_21_09_30, computo

end

*	Read in TSE data for 2019
capture: program drop getTSE
program define getTSE

	syntax, [TREP] [COMPUTO] [REPLACE]
	
	local count `trep'`computo'
	if ("`count'"=="trepcomputo") {
		local count
	}
	if ("`count'"=="") {
		di as error "Must choose one of TREP or COMPUTO"
		local count computo
	}
	if ("`count'"=="trep") {
		local upfile ${DATA}/TSE/datos_trep_computo_eg2019/trep/OEA_TREP_LogCompleto_19G_2019.11.03.06.53.xlsx
		local sheetname sheet1
	}
	else {
		local upfile ${DATA}/TSE/datos_trep_computo_eg2019/computo/2.RecepcionSobres_final.xlsx
		local sheetname sheet1
	}
	
	capture: use ${INPUT}/DTA/`count', clear
	if (_rc~=0 | "`replace'"=="replace") {
	
		capture: mkdir ${INPUT}/DTA
	
		*	Unpack XLSX to a working directory
		capture: mkdir ${INPUT}/TSE
		unpackXLSX `upfile', to(${INPUT}/TSE/`count')

		*	Read in the data
		getWorksheet, from(${INPUT}/TSE/`count')
		
		*
		*	Specific cleaning
		*

		*	Detect and decode any timestamps
		foreach var of varlist * {
			if (regexm("`var'","Date|Fecha") | regexm(`var'[1],"^[0-9][0-9][0-9][0-9][0-9]\.[0-9]+$")) {
				destring `var', force replace
				replace `var' = round((`var'+td(30dec1899))*86400)*1000
				format %tcDDmonCCYY_HH:MM:SS.sss `var'
			}
		}
		
		*	Destring election numbers
		destring NumMesa CodVer CC-Ins, force replace
		
		*
		*	SKIP THIS NEXT STEP. Nooruddin does NOT assign zero to missing vote data.
		*	That's obviously a bad thing, but he did what he did
		/*
		foreach var of varlist CC-PANBOL Bla Nul {
			replace `var' = 0 if mi(`var')
		}
		*/
		
		*	Remove variable labels (formerly adjusted cell columns)
		foreach var of varlist * {
			capture: label var `var' ""
		}
		
		compress
		save ${INPUT}/DTA/`count', replace
		
	}

end

*	End Specific Data Utilities

*
*	Miscellaneous Utilities
*

*	Do a fancy plot thingy
capture: program drop plotlp
program define plotlp, rclass

	syntax varlist [if/] [in/] [aweight] [, BY(varname) Range(numlist ascending min=1 max=2) PREplot(string asis) POSTplot(string asis)]
	
	marksample touse

	local vl: list sizeof varlist
	local xvar: word `vl' of `varlist'
	
	local rl: list sizeof range
	local r0 .
	local r1 .
	if (`rl') { 
		local r0: word 1 of `range'
		local r1 `r0'
		if (`rl'>1) {
			local r1: word 2 of `range'
		}
	}
	
	sum `xvar' if `touse'
	local gmin = min(`r0',r(min))
	local gmax = max(`r1',r(max))
	local gran = `gmax'-`gmin'
	
	tempvar gv
	gen `gv' = `gmin'+`gran'*(_n-1)/99 in 1/100
	
	
	local wcom `weight'
	if ("`wcom'"~="") {
		local wcom [`weight'`exp']
	}
	local lcom twoway `preplot'
	local rcolor ceprorange
	forvalues vv=1/`=`vl'-1' {
		tempvar xs0 ss0 xs1 ss1 lb0 ub0 lb1 ub1 acum dcum
		local yvar: word `vv' of `varlist'
		if ("`by'"~="") {
		
			lpoly `yvar' `xvar' `wcom' if ~`by', deg(1) gen(`xs0') se(`ss0') at(`gv') nograph
			gen `lb0' = `xs0'+invnormal(0.025)*`ss0'
			gen `ub0' = `xs0'+invnormal(0.975)*`ss0'
			
			gen `acum' = .
			forvalues ii=1/100 {
				count if ~`by' & `xvar'<=`gv'[`ii']
				replace `acum' = `r(N)' if _n==`ii'
			}
			gen `dcum' = `acum'-`acum'[_n-1]

			local lpcom
			forvalues ii=0(1)99 {
				sum `dcum' if ~mi(`dcum')
				local dcmax = r(max)
				sum `gv' if `dcum'>=`dcmax'*`ii'/100 & ~mi(`dcum')
				local a1 = r(max)
				sum `gv' if `gv'<r(min)
				local a0 = r(max)
				local lpcom `lpcom' (rarea `lb0' `ub0' `gv' if inrange(`gv',`a0',`a1'), lw(none) color(ceprgreen%2))
			}
			local lcom `lcom' `lpcom'
			
			lpoly `yvar' `xvar' `wcom' if `by', deg(1) gen(`xs1') se(`ss1') at(`gv') nograph
			gen `lb1' = `xs1'+invnormal(0.025)*`ss1'
			gen `ub1' = `xs1'+invnormal(0.975)*`ss1'
			
			forvalues ii=1/100 {
				count if `by' & `xvar'<=`gv'[`ii']
				replace `acum' = `r(N)' if _n==`ii'
			}
			replace `dcum' = `acum'-`acum'[_n-1]

			local lpcom
			forvalues ii=0(1)99 {
				sum `dcum' if ~mi(`dcum')
				local dcmax = r(max)
				sum `gv' if `dcum'>=`dcmax'*`ii'/100 & ~mi(`dcum')
				local a1 = r(max)
				sum `gv' if `gv'<r(min)
				local a0 = r(max)
				local lpcom `lpcom' (rarea `lb1' `ub1' `gv' if inrange(`gv',`a0',`a1'), lw(none) color(ceprorange%2))
			}
			local lcom `lcom' `lpcom'
		}
		else {

			lpoly `yvar' `xvar' `wcom', deg(1)  gen(`xs0') se(`ss0') at(`gv') nograph
			gen `lb0' = `xs0'+invnormal(0.025)*`ss0'
			gen `ub0' = `xs0'+invnormal(0.975)*`ss0'
			
			gen `acum' = .
			forvalues ii=1/100 {
				count if `xvar'<=`gv'[`ii']
				replace `acum' = `r(N)' if _n==`ii'
			}
			gen `dcum' = `acum'-`acum'[_n-1]

			local lpcom
			forvalues ii=0(1)99 {
				sum `dcum' if ~mi(`dcum')
				local dcmax = r(max)
				sum `gv' if `dcum'>=`dcmax'*`ii'/100 & ~mi(`dcum')
				local a1 = r(max)
				sum `gv' if `gv'<r(min)
				local a0 = r(max)
				local lpcom `lpcom' (rarea `lb0' `ub0' `gv' if inrange(`gv',`a0',`a1'), lw(none) color(`rcolor'%1))
			}
			local lcom `lcom' `lpcom'
			local rcolor ceprgreen
		}
	}
	local lcom `lcom' `postplot'
	macro li _lcom
	`lcom', `name'

end


*	Find all data sets and indicates whether or not each data set contains an update to the previous
capture: program drop findUpdates
program define findUpdates

	args count
	local ccount = ustrtitle("`count'")

	capture: use ${INPUT}/all`ccount'Updates, clear
	if (_rc~=0) {
		*	get list of files, hashes
		capture: use ${INPUT}/all`ccount'Hashes, clear
		if (_rc~=0) {
			local tfiles: dir "${DATA}/OEP/reports_`count'" files "acta.*.xlsx"
			clear
			set obs `: list sizeof tfiles'
			gen xlsfile = ""
			local rn 1
			foreach file of local tfiles {
				replace xlsfile = `"${DATA}/OEP/reports_`count'/`file'"' in `rn'
				local ++rn
			}
			gen timestamp = clock(substr(xlsfile,24+strlen("`count'"),19),"YMDhms")
			format %tc timestamp
			sort timestamp
			tostring timestamp, generate(dtafile) force format(%tcCCYY!_NN!_DD!_HH!_MM!_SS)
			replace dta = "${INPUT}/OEP/reports_`count'/"+dtafile+".dta"
			hash xlsfile, gen(hash) kind(md5) file
			order times xl dt hash
			save ${INPUT}/all`ccount'Hashes, replace
		}
		li in 1/10
		duplicates li hash

		*	produce any missing dta files
		capture: mkdir ${INPUT}/OEP
		capture: mkdir ${INPUT}/OEP/reports_`count'
		qui forvalues i=1/`=_N' {
			noi di as text "File: " as result `i' as text " of " as result `=_N'
			local tss: di %tcCCYY!_NN!_DD!_HH!_MM!_SS timestamp[`i']
			noi macro li _tss
			capture: confirm file `=dtafile[`i']'
			if (_rc~=0) {
				preserve
				getUpdate `tss' , `count'
				clear
				restore
			}
		}

		sort timestamp
		*	list files with new info
		preserve
		local cl0 = "${INPUT}/OEP/reports_computo/2019_10_20_00_23_02.dta"
		local cl
		forvalues i=1/`=_N' {
			local cl `cl' `=dtafile[`i']'
		}
		local cl: list cl-cl0
		local rclist 0
		local ufiles `:word 1 of `cl''
		use `:word 1 of `cl'', clear
		forvalues ii=2/`:list sizeof cl' {
			local cn: word `ii' of `cl'
			di as text "File " as result `ii' as text " : " as result "`cn'"
			capture: cf _all using `cn', verbose
			local rcc = _rc
			local rclist: list rclist | rcc
			if (_rc==9) {
				local ufiles `ufiles' `cn'
				use `cn', clear
			}
		}
		macro li _rclist
		di `: list sizeof ufiles'
		restore
		gen update = 0
		foreach cn of local ufiles {
			replace update = 1 if dtafile=="`cn'"
		}
		save ${INPUT}/all`ccount'Updates, replace
	}
	
end

*	End Miscellaneous Utilities

*	End document
