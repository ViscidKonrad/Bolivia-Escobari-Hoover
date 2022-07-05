
*
*	This file contains programs pertaining to the assembly of data
*

*	merge padron and  final results
capture: program drop mergePadronFinal
program define mergePadronFinal

	tempfile bfm
	buildFullMatches
	save `bfm'

	tempfile px
	getPadron 2019, ex
	ren pa idPai
	ren NombrePa País
	ren dep idDep
	ren nomdep Departamento
	ren PROV idPro
	ren nomprov Provincia
	ren sec idMun
	ren nombremun Municipio
	ren idloc idLoc
	ren Asi Localidad
	ren dist idDis
	ren nomdist Distrito
	ren zona idZon
	ren nomzona Zona
	drop CIR
	ren RECI idRec
	ren nombrerec Recinto
	ren cant mesas
	drop Inh Depur
	compress
	li in 1
	save `px'

	qui {
	tempfile gfc
	getFinalComputo
	ren Númerodep idDep
	ren Númeromun idMun
	drop Cir Ele Estado
	gen idPai = substr(CódigoMesa,3,3), before(Pa)
	gen cidDep = substr(CódigoMesa,6,2)
	gen idLoc = substr(CódigoMesa,8,4), before(Loc)
	gen idRec = substr(CódigoMesa,12,5), before(Rec)
	gen mesa = substr(CódigoMesa,17,2), after(Rec)
	destring idPai cidDep idLoc idRec mesa, force replace
	count if cidDep~=idDep
	if (r(N)==0) {
		drop cidDep
	}
	compress
	preserve
	keep if Pa=="Bolivia"
	save `gfc'
	restore
	keep if Pa~="Bolivia"
	}
	li in 1

	tempfile fullData
	merge n:1 idPai idDep idLoc idRec using `px', nogen
	save `px', replace
	bys idPa idDep idLoc idRec: egen tmes = total(1)
	bys idPa idDep idLoc idRec: egen tIns = total(Ins)
	count if tIns~=Hab | tmes~=mesas
	if (r(N)==0) {
		drop tmes tIns mesas Hab
	}
	compress
	li in 1
	save `fullData'

	getPadron 2019
	gen País = "Bolivia", before(dep)
	lab var Pa "Name: País"
	gen idPai = 32, before(Pa)
	lab var idPa "Code: País"
	ren dep idDep
	lab var idDep "Code: Departamento"
	ren nomdep Departamento
	lab var Depa "Name: Departmento"
	ren PROV idPro
	lab var idPro "Code: Provincia"
	ren nomprov Provincia
	lab var Prov "Name: Provincia"
	ren sec idMun
	lab var idMun "Code: Municipio"
	ren nombremun Municipio
	lab var Muni "Name: Municipio"
	ren idloc idLoc
	lab var idLoc "Code: Localidad"
	ren Asi Localidad
	lab var Local "Name: Localidad"
	ren dist idDis
	lab var idDis "Code: Distrito"
	ren nomdist Distrito
	lab var Dist "Name: Distrito"
	ren zona idZon
	lab var idZon "Code: Zona"
	ren nomzon Zona
	lab var Zona "Name: Zona"
	drop CIR
	ren REC idRec
	lab var idRec "Code: Recinto"
	ren nombrerec Recinto
	lab var Reci "Name: Recinto"
	ren nummesa CódigoMesa
	lab var CódigoMesa "Long mesa code"
	ren MESA mesa
	lab var mesa "Mesa number within precinct"
	li in 1
	
	merge 1:1 Pa idDep idLoc idRec mesa using `gfc', keep(match) nogen
	li Pa idDep idLoc idRec mesa Ins Hab-Depu if Hab~=Ins
	drop Hab-Depu
	lab var Ins "Eligible voters on acta"
	foreach var of varlist CC-PANB Bla Nul {
		lab var `var' "Official votes: `:var label `var''"
	}
	lab var VotosV "(unofficial) valid vote total reported on acta"
	li in 1

	append using `fullData'
	merge n:1 idPa idDep idLoc idRec using `bfm', nogen
	compress
		
end

*	Adds a few details: Chumacero's assignments, capital city localities
capture: program drop addMoreGeo
program define addMoreGeo

	tempfile fulldata
	li in 1
	save `fulldata'

	readChumacero
	li in 1

	merge 1:1 NúmeroMesa using `fulldata', nogen
	
	local rcond & RUE==1
	local dlist Beni Chuquisaca Cochabamba Paz Oruro Pando Potosí Cruz Tarija
	local mlist Trinidad Sucre Cochabamba Señora Oruro Cobija Potosí Sierra Tarija
	lab def isCapital 0 "Not a capital" 1 "Capital"
	gen isCapital:isCapital = 0, after(Loc)
	forvalues ii=1/`:list sizeof dlist' {
		replace isCap = regexm(Dep,"`:word `ii' of `dlist''") & regexm(Loc,"`:word `ii' of `mlist''") `rcond' if isCap==0
	}
	lab var isCap "1 = Capital Locality"
	
	lab def isBig4 0 "Not in biggest" 1 "In 4 biggest"
	gen isBig4:isBig4 = (ustrregexm(Loc,"Cocha|Nuestra|Sierra") & isCap) | Loc=="El Alto", after(isCap)
	lab var isBig4 "1 = Biggest 4 localities"
	
	compress
	
end

*
capture: program drop addPublicTiming
program define addPublicTiming

	preserve
	
	tempfile gft
	getFinalTREP
	keep NúmeroMesa
	save `gft'
	li in 1
	
	getAnnouncementTSE
	keep NúmeroMesa CC-Nulos
	ren F p21F
	ren * TSE_*
	ren TSE_NúmeroMesa NúmeroMesa
	foreach var of varlist TSE_* {
		lab var `var' `"Announced: `:var label `var''"'
	}
	lab var TSE_VotosV `"Announced: (unofficial) valid vote total reported on acta"'
	merge 1:1 NúmeroMesa using `gft', gen(timing)
	save `gft', replace
	li in 1
	
	restore
	merge 1:1 NúmeroMesa using `gft', nogen
	lab def reported 1 "Announced" 0 "After announcement" -1 "Not in public TREP"
	gen reported:reported = cond(timing==3,1,cond(timing==2,0,-1)), after(NúmeroMesa)
	lab var reported "Public status in TREP"
	drop timing
 
	compress
	li in 1

end

*
capture: program drop addTSE
program define addTSE

	preserve
	
	tempfile gt
	getTSE, computo
	keep if ustrregexm(Ele,"Pres")
	keep NumMesa CC-Ins Fecharecepción Fechaapertura Fechacustodia RegComputoDate ComputoDate AprobComputoDate PrimeraTransmisiónImagen UltimaTransmisiónImagen
	*	Computo timestamps
	lab var Fecharecepción				"Computo: Envelope received at TED"
	lab var Fechaapertura				"Computo: Envelope opened at TED"
	lab var PrimeraTransmisiónImagen	"Computo: First time image scanned"
	lab var UltimaTransmisiónImagen		"Computo: Last time image scanned"
	lab var Fechacustodia				"Computo: Archived"
	lab var RegComputoDate				"Computo: First transcription"
	lab var ComputoDate					"Computo: Second transcription"
	lab var AprobComputoDate			"Computo: Review if transcriptions do not match"
	
	compress
	save `gt'
	sum NumMesa in 10000
	local na = r(mean)
	li in 10000

	getTSE, trep
	keep if ustrregexm(Ele,"Pres")
	keep NumMesa Estado FechaReg PriRegistroDate UltRegistroDate PriTransmisionDate UltTransmisionDate VerificadorDate AprobadorDate
	*	TREP timestamps
	lab var EstadoActa			"TREP: Last action taken on acta"
	lab var FechaRegistroenLog	"TREP: (EstadoActa) timestamp"	
	lab var PriRegistroDate		"TREP: First time on phone when taken/sent"
	lab var UltRegistroDate		"TREP: Last time on phone when taken/sent"
	lab var PriTransmisionDate	"TREP: First time image and transcription received"
	lab var UltTransmisionDate	"TREP: Last time image and transcription received"
	lab var VerificadorDate		"TREP: Transcribed by SERECI"
	notes VerificadorDate: Two actas listed in the 2019.10.20.19.40.57 report were Verified after that time
	notes VerificadorDate: 2433  was (originally) verified at 13:21:05
	notes VerificadorDate: 61373 was verified at 19:40:59 but was included in the 19:40:57 report
	lab var AprobadorDate		"TREP: Reviewed by SERCEI"
	notes: This data only includes the LAST entry for each acta
	notes: The 2019.10.20.19.40.57 report includes all most recently verified no later than 20oct2019 19:40:57, but votes differ on 2433
	notes: 2433 was originally APPROVED at 18:32:26 but with zero valid votes (though registered at 12:05:42 with votes, restored on the 22nd at 09:16:16)

	li if NumMesa==`na'
	bys NumMesa (FechaReg): keep if _n==_N
	compress

	merge 1:1 NumMesa using `gt', nogen
	ren NumMesa NúmeroMesa
	ren MAS MASIPSP
	foreach var of varlist CC-Ins {
		replace `var' = 0 if mi(`var')
		ren `var' tse_`var'
	}
	save `gt', replace
	
	restore
	ren F p21F
	merge 1:1 NúmeroMesa using `gt', nogen
	
	count if VotosV~=tse_Validosen
	local terr = r(N)
	foreach var of varlist Ins-PANBOL Bla Nul {
		count if `var'~=tse_`var'
		local terr = `terr'+r(N)
	}
	if (`terr'==0) {
		drop tse_*
	}
	sum

end

*	add Districto and Locality-level 2016 results for 2016 and show some results
capture: program drop doBig4
program define doBig4

	preserve
	levelsof idLoc if isBig4, local(llist) clean

	tempfile pd
	prepDomestic2016
	keep if inlist(IdLoc,`=usubinstr("`llist'"," ",",",.)')
	collapse (sum) NO SÍ, by(idDep IdLoc Dist)
	bys idDep IdLoc: egen loc_NO = total(NO)
	bys idDep IdLoc: egen loc_SÍ = total(SÍ)
	ren NO dis_NO
	ren SÍ dis_SÍ
	gen idPai = 32
	ren IdLoc idLoc
	ren Dist idDis
	compress
	save `pd'
	
	restore, preserve
	merge n:1 idPa idDep idLoc idDis using `pd', nogen
	
	collapse (sum) MAS CC sumV (firstnm) dis* loc* if isBig4 & reported==1, by(idLoc Loc idDis Dist)
	bys idLoc: egen loc_MAS = total(MAS)
	bys idLoc: egen loc_CC = total(CC)
	bys idLoc: egen loc_sumV = total(sumV)
	gen dmargin = 100*(MAS-CC)/sumV
	gen lmargin = 100*(loc_MAS-loc_CC)/loc_sumV
	gen drmargin = 100*(dis_S-dis_N)/(dis_S+dis_N)
	gen lrmargin = 100*(loc_S-loc_N)/(loc_S+loc_N)
	egen ltag = tag(idLoc)
	gen sn = Loc
	replace sn = "La Paz" if ustrregexm(Loc,"Nuestra")
	replace sn = "Santa Cruz" if ustrregexm(Loc,"Cruz")
	lab var sn "Locality"
	twoway (scatter lmargin lrmargin if ltag, mlw(none) mc(ceprgreen%50) msize(huge)) ///
		(scatter dmargin drmargin, mlw(none) mc(ceprorange%50) msize(small)) ///
			, aspect(1) legend(order(1 2) label(1 "Entire Locality") label(2 "Districts within Locality")) ///
			xti(" " "{it:Sí} Margin in 2016") yti("MAS Margin in 2019") by(sn, caption("") note(""))
	restore
			

end

capture: program drop doSplit
program define doSplit

	local dvar before_type
	tab `dvar' [w=sumV]

	preserve
	collapse (sum) Ins sumV (firstnm) r_margin, by(recc `dvar')
	gen x = (_n-101) in 1/201
	forvalues i=0/3 {
		kdensity r_margin [aw=sumV] if `dvar'==`i'-1, generate(nx`i' nd`i') at(x) nograph
		*sum Ins if `dvar'==`i'-1
		*replace nd`i' = nd`i'*r(sum)
	}

	forvalues i=0/3 {
		sum r_margin [fw=sumV] if `dvar'==`i'-1
		local ir`i' = r(mean)
		macro li _ir`i'
	}
	twoway (area nd0 nx0, color(ceprorange%25) lw(none)) ///
		(area nd3 nx3, color(ceprgreen%25) lw(none)) ///
		, legend(label(1 "Early Stations") label(2 "Late Stations")) ysc(off) xti("{it:Sí} margin in 2016") ///
		xline(`ir0', lc(ceprorange)) xline(`ir3', lc(ceprgreen)) name(alls) ti("Early/Late Precincts")
	twoway (area nd1 nx1, color(ceprorange%25) lw(none)) ///
		(area nd2 nx2, color(ceprgreen%25) lw(none)) ///
		, legend(label(1 "Early Stations") label(2 "Late Stations")) ysc(off) ///
		xline(`ir1', lc(ceprorange)) xline(`ir2', lc(ceprgreen)) name(splits) ti("Split Precincts")
	grc1leg2 alls splits, xtob xti(alls) ring(1) ycommon
	restore

end
