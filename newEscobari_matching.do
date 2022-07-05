
*
*	This file contains the programs to produce a crosswalk between 2016 and 2019 precincts
*

*	Merge domestic 2016 data with padron at the RECINTO level
capture: program drop prepDomestic2016
program define prepDomestic2016

	*	Recode padron for matching
	tempfile gp gpm
	getPadron 2016
	gen idDep = Dep
	gen idPro = Pro
	gen idMun = Sec
	gen mergeName = ustrupper(NombreRec)
	bys idDep idPro idMun mergeName (IdLoc Reci): gen g = sum(IdLoc~=IdLoc[_n-1] | Reci~=Reci[_n-1])
	bys idDep idPro idMun mergeName g (nummesa): gen Mesa = _n
	compress
	li if idDep==4 & IdLoc==3021 & Reci==2
	save `gp'
	
	*	Collapse padron into precincts
	collapse (count) Mesas=Ciu (sum) TotalC=Ciu, by(idDep idPro idMun mergeName g)
	reshape wide Mesas TotalC, i(idDep idPro idMun mergeName) j(g)
	save `gpm'
	
	*	Recode results
	getDatosAbiertos 2016
	gen Ciudadanos = INS
	replace Ciu = Ciu+1 if inlist(CodigoMESA,27639,27864,32483,33651,73681,73685,72740,71277,74561,80672)
	gen idDep = CodigoDEP
	gen idPro = CodigoPRO-100*CodigoDEP
	gen idMun = CodigoMUN-100*CodigoPRO
	gen mergeName = ustrupper(NombreREC)
	
	*	Match precincts
	merge n:1 idDep idPro idMun mergeName using `gpm', nogen
	gen g = .
	forvalues ii=1/3 {
		bys idDep idPro idMun mergeName g (CodigoMESA): gen msum = sum(1)
		bys idDep idPro idMun mergeName g (CodigoMESA): gen isum = sum(Ciu)
		bys idDep idPro idMun mergeName g (CodigoMESA): egen t1 = total(msum==Mesas1 & isum==TotalC1) if msum<=Mesas1
		bys idDep idPro idMun mergeName g (CodigoMESA): egen t2 = total(msum==Mesas2 & isum==TotalC2) if msum<=Mesas2
		bys idDep idPro idMun mergeName g (CodigoMESA): egen t3 = total(msum==Mesas3 & isum==TotalC3) if msum<=Mesas3
		replace g = cond(t1==1,1,cond(t2==1,2,cond(t3==1,3,.),.),.)
		capture: drop msum isum t1 t2 t3
	}
	drop Mesas1-TotalC3 Ciu
	bys idDep idPro idMun mergeName g (CodigoMESA): gen Mesa = _n
	compress

	*	Merge back with padron and clean
	merge 1:1 idDep idPro idMun mergeName g Mesa using `gp', nogen
	
	*	Clean up and collapse into precincts
	order CodigoPA NombrePA idDep Dep NomDep CodigoDEP NombreDEP idPro Pro NomPro CodigoPRO NombrePRO idMun Sec NombreMun CodigoMUN ///
		IdLoc Asi Dist NomDis Zon NomZon Reci NombreRec NombreREC ///
		CodigoREC Mesa nummesa CodigoMESA NombreMESA Ciu NO-VAL
	drop mergeName g
	collapse (count) Mesa (max) Mesas=Mesa (sum) Ciu-VAL, by(CodigoPA-NombreREC)
	count if Mesa~=Mesas
	drop Mesa

end

*	Merge foreign 2016 data with padron at the MESA level
capture: program drop prepForeign2016
program define prepForeign2016

	*	Recode padron for matching
	tempfile gp
	getPadron 2016, ex
	gen mergeNum = nummesa
	gen mergeC = Ciu
	save `gp'
	
	*	Recode results
	getDatosAbiertos 2016, ex
	gen mergeNum = CodigoMESA
	gen mergeC = INS
	
	*	Easy matching
	merge 1:1 mergeNum mergeC using `gp', nogen
	
	*	Clean and collapse into precincts
	drop mergeNum mergeC
	order idpa nompa CodigoPA NombrePA dep-nombrerec NombreREC CodigoREC numm CodigoMES NombreMES ///
		Ciu NO-VAL
	bys idpa-NombreREC (CodigoMES): gen Mesa = _n
	order Mesa, before(numm)
	compress
	collapse (count) Mesa (max) Mesas=Mesa (sum) Ciu-VAL, by(idpa-NombreREC)
	count if Mesa~=Mesas
	drop Mesa

end

*	Match foreign 2016 to foreign 2019
capture: program drop matchForeign
program define matchForeign

	tempfile mf
	prepForeign2016
	gen idPai = idpa
	gen idDep = dep
	gen idLoc = idloc
	gen idDis = dis
	gen idZon = zon
	gen idRec = reci
	gen ndep = nomdep
	gen nn = NombreREC
	foreach var of varlist NO-VAL {
		ren `var' x`var'
	}

	*	Argentina > Buenos Aires > Esc Primaria Nro 3 Jose manuel Estrada gets its own Locality
	replace idLoc = 4129 if idPa==11 & dep==10 & idLoc==3611 & idRec==25192

	*	Argentina > Buenos Aires > Distrito Provincia Buenos Aires gets broken up
	replace idLoc = cond(zona==206, cond(inrange(idRec,25060,25178),4130, ///
											cond(idRec==25179,4131,cond(idRec==25180,4132,4132+idRec-25181))), ///
					cond(zona==207, 4135, ///
					cond(zona==208, cond(idRec==25062, 4136, cond(idRec==25184, 4137, 4135)), ///
					cond(zona==209, cond(idRec==25063, 4138, 4139), ///
					cond(zona>=210, 4140+zona-210, ///
					idLoc))))) if idPa==11 & dep==10 & idLoc==3611

	*	Argentina > Mendoza > Distrito Mendoza gets broken up
	replace idLoc = 4143+idRec-25148 if idPa==11 & dep==11 & idLoc==3612 & idRec>25065

	*	Argentina > Cordoba > Zona Villa Maria gets its own Locality
	replace idLoc = 4147 if idPa==11 & dep==21 & idLoc==3762 & zona==103

	*	Brazil > Sao Paulo > Distrito Sao Paulo gets broken up
	replace idLoc = cond(zona==103, 4151, ///
					cond(zona==106, cond(idRec==25083, 4155, idLoc), ///
					idLoc)) if idPa==34 & dep==18 & idLoc==3619

	*	USA > Virginia > Distrito Virginia gets broken up
	replace idLoc = 4156 if idPa==62 & dep==13 & idLoc==3614 & idRec==25068
					
	save `mf'

	getPadron 2019, ex
	gen idPai = pa
	gen idDep = dep
	gen idLoc = idloc
	gen idDis = dis
	gen idZon = zon
	gen idRec = RECI
	gen nn = nombrerec

	merge 1:1 idPa idLoc idRec using `mf', nogen update
	compress
	sort idPa dep idLoc dist zona idRec Ciu
	order idPa idDep idLoc idDis idZon idRec
	li NombrePa idPa-idRec nn Ciu Hab, sepby(idPa dep)

	bys idPai-idZon: egen midCount = total(mi(Hab))
	bys idPai-idZon: egen midRec = mean(cond(mi(Hab),idRec,.)) if midCount==1

	*	unique id for each observation
	gen ID = _n

	*	full precinct matches
	gen mRecID = ID if ~mi(Ciu) & ~mi(Hab)
	*	clear precinct matches
	bys idPai-idZon midRec (Ciu): replace mRecID = ID[1] if mi(mRecID) & ~mi(Hab) & ~mi(Ciu[1])

	gen idFor = 1, before(idPai)
	local lg Rec
	local geos Zon Dis Loc Dep Pai For
	local gcodes zon dis idloc dep idpa idFor
	forvalues n=1/6 {
		local g: word `n' of `geos'
		local d: word `n' of `gcodes'
		*	Match to paired precincts if possible
		sort ID
		gen m`g'ID = m`lg'ID
		gen tmp`g' = `d'[m`lg'ID]
		gen tmpCiu = Ciu[m`lg'ID]
		bys idFor-id`g' tmp`g': egen tmpINS = total(-tmpCiu)
		bys idFor-id`g' (tmpINS): replace m`g'ID = m`lg'ID[1] if mi(m`g'ID) & ~mi(Hab)
		drop tmpCiu tmpINS tmp`g'
		*	Match to whatever is available in the area
		sort ID
		bys idFor-id`g' `d': egen tmpINS = total(-Ciu)
		bys idFor-id`g' (tmpINS): replace m`g'ID = ID[1] if mi(m`g'ID) & ~mi(Hab) & tmpINS<0
		drop tmpINS
		local lg `g'
	}
	
	sort ID
	li idPa dep idLoc dis zon idRec nn Ciu Hab midRec ID mRecID mZonID mDisID mLocID mDepID mPaiID mForID, sepby(idPa)

	count if ~mi(Hab)
	count if ~mi(Hab) & mi(Ciu)
	count if ~mi(Hab) & mi(mRecID)
	count if ~mi(Hab) & mi(mZonID)
	count if ~mi(Hab) & mi(mDisID)
	count if ~mi(Hab) & mi(mLocID)
	count if ~mi(Hab) & mi(mDepID)
	count if ~mi(Hab) & mi(mPaiID)
	count if ~mi(Hab) & mi(mForID)

end

*	Add lat/long to 2019 data
capture: program drop addGeo2019
program define addGeo2019

	tempfile gp
	getPadron 2019
	gen idDep = dep
	gen idLoc = idloc
	gen idRec = RECI
	collapse (count) MESAS=MESA (sum) INSCRITOS=Hab, by(idDep-idRec dep-nombrerec)
	compress
	save `gp'
	getAllGeo EN19
	drop p_CC-GAN
	destring OBJ Dep Prov ProvSec Sec Loc-SecLoc IdLocRec Reci DistZona-Zona ZonaRec Cir ///
		MESAS-EMIT x y, force replace
	gen idDep = Dep
	gen idLoc = IdLoc
	gen idRec = Reci
	keep idDep-idRec x y
	merge 1:1 idDep-idRec using `gp', keep(match) nogen
	order x y, last

end

*	Just fills in some blanks in the 2016 lat/long
capture: program drop addCoordinates
program define addCoordinates

	*
	*	Add some location information from other sources
	*
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/80480017
	replace latit = -19.03800011 if Dep==1 & IdLoc==2794 & Reci==109
	replace longi = -65.25968933 if Dep==1 & IdLoc==2794 & Reci==109
	*	https://www.google.com/maps/place/School+of+Medicine/@-19.0460397,-65.2683243,17z/data=!4m12!1m6!3m5!1s0x93fbcf4bbc671ce7:0xb4133184e494c22a!2sSchool+of+Medicine!8m2!3d-19.0460397!4d-65.2661356!3m4!1s0x93fbcf4bbc671ce7:0xb4133184e494c22a!8m2!3d-19.0460397!4d-65.2661356
	replace latit = -19.045826715393662 if Dep==1 & IdLoc==2794 & Reci==518
	replace longi = -65.26611414421708 if Dep==1 & IdLoc==2794 & Reci==518
	*	https://gitlab.agetic.gob.bo/adsib/geoelectoral-api/-/blob/agetic-v1.3.1/public/elecciones_2016_26.csv
	replace latit = -19.03716 if Dep==1 & IdLoc==2794 & Reci==525
	replace longi = -65.24286 if Dep==1 & IdLoc==2794 & Reci==525
	*	https://www.google.com/maps/place/Sucre+Alcantari+International+Airport+(SRE)/@-19.2463533,-65.1532287,17z/data=!3m1!4b1!4m5!3m4!1s0x93fbd4534a037e1b:0xdab9c59c8a72e36d!8m2!3d-19.2463533!4d-65.15104
	replace latit = -19.246160832053977 if Dep==1 & IdLoc==2794 & Reci==218
	replace longi = -65.15101854421515 if Dep==1 & IdLoc==2794 & Reci==218
	*	(alternatively) https://gitlab.agetic.gob.bo/adsib/geoelectoral-api/-/blob/agetic-v1.3.1/public/elecciones_2016_26.csv
	replace latit = -19.01331 if Dep==1 & IdLoc==2794 & Reci==218
	replace longi = -65.29323 if Dep==1 & IdLoc==2794 & Reci==218
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/80730442
	replace latit = -16.5028991 if Dep==2 & IdLoc==729 & Reci==406
	replace longi = -68.12283325 if Dep==2 & IdLoc==729 & Reci==406
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/80860057
	replace latit = -18.26210022 if Dep==3 & IdLoc==1066 & Reci==1309
	replace longi = -65.38819122 if Dep==3 & IdLoc==1066 & Reci==1309
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/80860056
	replace latit = -18.3022995 if Dep==3 & IdLoc==1068 & Reci==1310
	replace longi = -65.37464142 if Dep==3 & IdLoc==1068 & Reci==1310
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/80860011
	replace latit = -18.19070053 if Dep==3 & IdLoc==1069 & Reci==1312
	replace longi = -65.38265228 if Dep==3 & IdLoc==1069 & Reci==1312
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/80860064
	replace latit = -18.14959908 if Dep==3 & IdLoc==3361 & Reci==9062
	replace longi = -65.33539581 if Dep==3 & IdLoc==3361 & Reci==9062
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/80860030
	replace latit = -18.20359993 if Dep==3 & IdLoc==3866 & Reci==1
	replace longi = -65.32512665 if Dep==3 & IdLoc==3866 & Reci==1
	*	https://gitlab.agetic.gob.bo/adsib/geoelectoral-api/-/blob/agetic-v1.3.1/public/elecciones_2016_26.csv
	replace latit = -17.93559 if Dep==4 & IdLoc==2059 & Reci==25041
	replace longi = -67.08788 if Dep==4 & IdLoc==2059 & Reci==25041
	*	https://gitlab.agetic.gob.bo/adsib/geoelectoral-api/-/blob/agetic-v1.3.1/public/elecciones_2016_26.csv
	replace latit = -21.51948 if Dep==5 & IdLoc==1856 & Reci==462
	replace longi = -67.35812 if Dep==5 & IdLoc==1856 & Reci==462
	*	https://www.google.com/maps/place/Notar%C3%ADa+No.+4/@-21.2641343,-63.4725606,17z/data=!4m12!1m6!3m5!1s0x94090519e9c33201:0xa3eeb671e1980876!2sNotar%C3%ADa+No.+4!8m2!3d-21.2641609!4d-63.4725681!3m4!1s0x94090519e9c33201:0xa3eeb671e1980876!8m2!3d-21.2641609!4d-63.4725681
	replace latit = -21.263702794259267 if Dep==6 & IdLoc==3509 & Reci==2
	replace longi = -63.47253292751586 if Dep==6 & IdLoc==3509 & Reci==2
	*	U.E. SIXTO SUAREZ NAVARRO probably merged in 2019
	sum latit if Dep==7 & IdLoc==2234 & Reci==11065
	replace latit = r(mean) if Dep==7 & IdLoc==2234 & Reci==25424
	sum longi if Dep==7 & IdLoc==2234 & Reci==11065
	replace longi = r(mean) if Dep==7 & IdLoc==2234 & Reci==25424
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/41980031
	replace latit = -17.9993 if Dep==7 & IdLoc==2257 & Reci==11024
	replace longi = -63.38502121 if Dep==7 & IdLoc==2257 & Reci==11024
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/41980095
	replace latit = -17.9878006 if Dep==7 & IdLoc==2257 & Reci==25406
	replace longi = -63.38105011 if Dep==7 & IdLoc==2257 & Reci==25406
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/61880173
	*	https://gitlab.agetic.gob.bo/adsib/geoelectoral-api/-/blob/agetic-v1.3.1/public/elecciones_2016_26.csv
	replace latit = -16.33947 if Dep==7 & IdLoc==4047 & Reci==25378
	replace longi = -63.06822 if Dep==7 & IdLoc==4047 & Reci==25378
	*	https://www.google.com/maps/place/INFOCAL/@-14.8233066,-64.9044994,17z/data=!3m1!4b1!4m5!3m4!1s0x93dd6e28fd389a05:0x64c0a234f31b2810!8m2!3d-14.8233066!4d-64.9023107
	replace latit = -14.823099154086536 if Dep==8 & IdLoc==2526 & Reci==25
	replace longi = -64.90228924425435 if Dep==8 & IdLoc==2526 & Reci==25
	*	http://seie.minedu.gob.bo/reportes/mapas_unidades_educativas/ficha/ver/82190131
	replace latit = -14.99527957 if Dep==8 & IdLoc==2596 & Reci==1
	replace longi = -65.64478475 if Dep==8 & IdLoc==2596 & Reci==1
	*	https://gitlab.agetic.gob.bo/adsib/geoelectoral-api/-/blob/agetic-v1.3.1/public/elecciones_2016_26.csv
	replace latit = -12.65806 if Dep==8 & IdLoc==2615 & Reci==1
	replace longi = -64.4332 if Dep==8 & IdLoc==2615 & Reci==1
	
end

*	Add lat/long to 2016 data
capture: program drop addGeo2016
program define addGeo2016

	tempfile gp
	prepDomestic2016
	gen idLoc = IdLoc, before(IdLoc)
	gen idRec = Reci, before(Reci)
	save `gp'
	
	getAllGeo RE16
	qui destring OBJ-sec idloc reci cir cant_m-y, force replace
	gen idDep = dep
	gen idLoc = idloc
	gen idRec = reci
	keep idDep idLoc idRec lati longi
	merge 1:1 idDep idLoc idRec using `gp', keep(match using) nogen
	addCoordinates
	LambertConformalEllipsoid lati longi, s(-21.5 -11.5) o(-24 -64) falsee(1000000) gen(x y)
	sort Dep Prov Sec IdLoc Dist Zon Reci


end

*	Do a little data/graph dance to check work if desired
capture: program drop showArea
program define showArea

	syntax [if/] [in/] [, ADDvars(varlist)]
	
	marksample touse

	li idLoc idDis idZon idRec nomRec Mesas MESAS Ciu INS `addvars' if `touse', sepby(idDep-idZon) h(30)
	twoway ///
			(scatter y x, m(Oh) mlw(vthin) /*mlab(zst) mlabs(tiny) mlabpos(12)*/) ///
			(scatter y x if mi(Ciu), m(Oh) mlw(vthin) mlab(nomRec) mlabs(tiny)) ///
			(scatter y x if mi(INS), m(Oh) mlw(vthin) mlab(nomRec) mlabpos(8) mlabs(tiny)) ///
			(pcarrow my mx y x if mi(Ciu)) ///
		if `touse', legend(off)

end

*	This goes through the data, checking on stuff that didn't match exactly
capture: program drop showTentativeMatches
capture: program drop showTentativeMatches
program define showTentativeMatches

	syntax [if/] [in/] [, ADDvars(varlist)]
	
	marksample touse

	pause on
	levelsof idDep if `touse', local(dl) clean
	foreach d of local dl {
		levelsof idPro if `touse' & idDep==`d', local(pl) clean
		foreach p of local pl {
			levelsof idMun if `touse' & idDep==`d' & idPro==`p', local(ml) clean
			foreach m of local ml {
				count if `touse' & idDep==`d' & idPro==`p' & idMun==`m' & mi(Ciu)
				if (`r(N)') {
					count if `touse' & idDep==`d' & idPro==`p' & idMun==`m' & mi(INS)
					if (`r(N)') {
						di `"showArea `d' `p' `m'"'
						if inlist(10000*`d'+100*`p'+`m',10101,20101,60101,60301,70101) {
							levelsof idLoc if `touse' & idDep==`d' & idPro==`p' & idMun==`m', local(ll) clean
							foreach l of local ll {
								count if `touse' & idDep==`d' & idPro==`p' & idMun==`m' & idLoc==`l' & mi(Ciu)
								if (`r(N)') {
									count if `touse' & idDep==`d' & idPro==`p' & idMun==`m' & idLoc==`l' & mi(INS)
									if (`r(N)') {
										showArea if `touse' & idDep==`d' & idPro==`p' & idMun==`m' & idLoc==`l', addvar(`addvars')
										pause
									}
								}
							}
						}
						else {
							showArea if `touse' & idDep==`d' & idPro==`p' & idMun==`m', addvar(`addvars')
							pause
						}
					}
				}
			}
		}
	}
	pause off
	
end

*	Build tentative matching database
capture: program drop matchDomestic
program define matchDomestic

	tempfile g16
	addGeo2016
	gen idDis = Dist, before(Dist)
	gen idZon = Zona, before(Zona)
	foreach var of varlist NO-VAL {
		ren `var' x`var'
	}
	
	*	These precincts appears to have merged in 2019
	li if regexm(NombreRec,"U.E. LIC. ZACHARY MACY")
	replace idRec = 25405 if Dep==7 & IdLoc==2257 & Reci==11024
	li if regexm(NombreRec,"SIXTO")
	replace idRec = 11065 if Dep==7 & IdLoc==2234 & Reci==25424
	collapse (firstnm) Reci (mean) lati longi x y (sum) Mesas-xVAL, by(idDep-NomZona NombreRec NombreREC)
	li if regexm(NombreRec,"SIXTO|ZACHARY")	

	*	Account for Raqaypampa breaking off from Mizque
	replace idMun = 4 if idDep==3 & idPro==13 & idMun==1 & inlist(idLoc,1066,1068,1069,3361,3866)
	
	*	Hard matches
	replace idLoc = 4227 if idDep==7 & idLoc==2329 & idRec==25445
	replace idLoc = 2773 if idDep==9 & idLoc==2750
	*

	*	Pretty indisputable
	replace idZon = 1001 if idDep==2 & idLoc==729 & idDis==8
	*
	
	*	Various gentle nudges
	replace idZon = 147 if idDep==3 & idLoc==821 & idRec==110
	replace idZon = 0 if idDep==3 & idLoc==3725 & idRec==1
	replace idDis = 0 if idDep==3 & idLoc==3725 & idRec==1
	replace idLoc = 4242 if idDep==3 & idLoc==3725 & idRec==1
	replace idZon = 257 if idDep==4 & idLoc==2059 & idZon==248		// Arguably could match Zona 261 instead
	replace idZon = 258 if idDep==4 & idLoc==2059 & idRec==25042
	replace idLoc = 4356 if idDep==6 & idLoc==1867
	replace idZon = 99 if idDep==6 & idLoc==2680 & idRec==91
	replace idZon = 91 if idDep==6 & idLoc==2680 & idZon==93
	replace idZon = 1 if idDep==6 & idLoc==1949 & idRec==14
	replace idLoc = 4357 if idDep==6 & idLoc==2005
	replace idZon = 2 if idDep==7 & idLoc==2342 & idRec==22005
	replace idDis = 0 if idDep==8 & idLoc==2573
	replace idZon = 2 if idDep==8 & idLoc==2573
	replace idZon = 0 if idDep==9 & IdLoc==1170
	replace idLoc = 4173 if idDep==9 & IdLoc==1170
	*
	
	*	Really just matching up stragglers
	replace idZon = 0 if idDep==9 & idLoc==1124
	replace idLoc = 4367 if idDep==9 & idLoc==1124
	*
		
	save `g16'

	addGeo2019
	gen idPro = PROV, after(idDep)
	gen idMun = sec, after(idPro)
	gen idDis = dist, after(idLoc)
	gen idZon = zona, after(idDis)
	gen nomRec = nombrerec

	merge 1:1 idDep idLoc idRec using `g16', nogen update
	replace nomRec = NombreRec if mi(nomRec)
	sort idDep-nomRec

	bys idDep-idZon: egen midCount = total(mi(INS))
	bys idDep-idZon: egen midRec = mean(cond(mi(INS),idRec,.)) if midCount==1
	bys idDep-idZon: egen my = mean(cond(idRec==midRec,y,.)) if midCount==1
	bys idDep-idZon: egen mx = mean(cond(idRec==midRec,x,.)) if midCount==1
	compress
	local addvars ID mRecID mZonID mDisID mLocID mMunID nmm ism

	*	unique id for each observation
	replace midRec = idRec if ~mi(Ciu) & ~mi(INS)
	sort idDep-idRec Ciu
	gen ID = _n

	*	full precinct matches
	gen mRecID = ID if ~mi(Ciu) & ~mi(INS)
	*	clear precinct matches
	bys idDep-idZon midRec (Ciu): replace mRecID = ID[1] if mi(mRecID) & ~mi(INS) & ~mi(Ciu[1])

	local lg Rec
	local geos Zon Dis Loc Mun 
	local gcodes Zon Dis IdLoc Sec 
	forvalues n=1/4 {
		local g: word `n' of `geos'
		local d: word `n' of `gcodes'
		sort ID
		gen m`g'ID = m`lg'ID
		gen tmp`g' = `d'[m`lg'ID]
		gen tmpCiu = Ciu[m`lg'ID]
		bys idDep-id`g' tmp`g': egen tmpINS = total(-tmpCiu)
		bys idDep-id`g' (tmpINS): replace m`g'ID = m`lg'ID[1] if mi(m`g'ID) & ~mi(INS)
		drop tmpCiu tmpINS tmp`g'
		local lg `g'
	}

	count if ~mi(INS)
	count if ~mi(INS) & mi(Ciu)
	count if ~mi(INS) & mi(mRecID)
	count if ~mi(INS) & mi(mZonID)
	count if ~mi(INS) & mi(mDisID)
	count if ~mi(INS) & mi(mLocID)
	count if ~mi(INS) & mi(mMunID)

	gen ism = ~mi(INS) & mi(mLocID)
	bys idDep-idMun: egen nmm = total(ism)
	sort ID
	*	This is something I used to help nudge the matching along
	*showTentativeMatches if nmm>0, addv(`addvars')

end

*	Find best matches to 2016 for each domestic 2019 precinct
capture: program drop addMatchedDomestic
program define addMatchedDomestic

	matchDomestic
	
	*	Generate 2016 aggregates
	local byl Dep Pro
	local addby Sec IdLoc Dis Zon Rec
	local nn Mun Loc Dis Zon Rec
	unab vl: xNO-xVAL
	forvalues i=1/5 {
		local byl `byl' `:word `i' of `addby''
		local g: word `i' of `nn'
		macro li _byl
		macro li _g
		foreach var of local vl {
			bys `byl': egen `var'`g' = total(`var')
		}
	}
	
	*	Match to 2019
	sort ID
	local nvars nombremun Asien nomdist nomzon nombrerec
	local nvars NombreMun AsientoEl NomDist NomZon NombreRec
	gen level2016 = ""
	gen match2016 = ""
	foreach var of local vl {
		gen m`var' = .
	}
	forvalues i=1/5 {
		local g: word `i' of `nn'
		replace level2016 = "`g'" if ~mi(m`g'ID)
		replace match2016 = `:word `i' of `nvars''[m`g'ID] if ~mi(m`g'ID)
		foreach var of local vl {
			replace m`var' = `var'`g'[m`g'ID] if ~mi(m`g'ID)
		}
	}
	tab level2016 if ~mi(INS), m

end

*	Find best matches to 2016 for each domestic 2019 precinct
capture: program drop addMatchedForeign
program define addMatchedForeign

	matchForeign
	gen NombreFor = "Foreign"
	
	*	Generate 2016 aggregates
	local byl
	local addby idFor idpa dep idloc dis zon rec
	local nn For Pai Dep Loc Dis Zon Rec
	unab vl: xNO-xVAL
	forvalues i=1/7 {
		local byl `byl' `:word `i' of `addby''
		local g: word `i' of `nn'
		macro li _byl
		macro li _g
		foreach var of local vl {
			bys `byl': egen `var'`g' = total(`var')
		}
	}
	
	*	Match to 2019
	sort ID
	local nvars NombreFor NombrePa nomprov Asi nomdist nomzon nombrerec
	gen level2016 = ""
	gen match2016 = ""
	foreach var of local vl {
		gen m`var' = .
	}
	qui forvalues i=1/7 {
		local g: word `i' of `nn'
		replace level2016 = "`g'" if ~mi(m`g'ID)
		replace match2016 = `:word `i' of `nvars''[m`g'ID] if ~mi(m`g'ID)
		foreach var of local vl {
			replace m`var' = `var'`g'[m`g'ID] if ~mi(m`g'ID)
		}
	}
	lab var level "2016/19 level of data match"
	lab var match "2016: Geography to which precinct is matched"
	lab var mxNO "2016: No votes in match2016"
	lab var mxS "2016: Sí votes in match2016"
	lab var mxB "2016: Blank votes in match2016"
	lab var mxE "2016: All votes cast in match2016"
	lab var mxI "2016: Eligible voters in match2016"
	lab var mxNU "2016: Null votes in match2016"
	lab var mxV "2016: Valid votes in match2016"
	tab level2016 if ~mi(Hab), m

end

*	This builds the final crosswalk of matches
capture: program drop buildFullMatches
program define buildFullMatches

	syntax [, REPLACE]

	capture: use ${INPUT}/DTA/matches, clear
	if (_rc~=0 | "`replace'"=="replace") {
	
		capture: mkdir ${INPUT}/DTA
	
		tempfile df
		addMatchedDomestic
		keep if ~mi(INS)
		keep idDep-idRec level2016-mxVAL
		gen idPai = 32, before(idDep)
		drop idPro idMun idDis idZon
		ren mx* m*
		order idPa-idRec
		compress
		li in 1
		save `df'

		qui addMatchedForeign
		keep if ~mi(Hab)
		keep idPai-idRec lev mat mx*
		drop idDis idZon
		ren mx* m*
		qui compress
		li in 1
		
		append using `df'
		save ${INPUT}/DTA/matches, replace
		
	}
	
end
