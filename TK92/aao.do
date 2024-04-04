

clear all
cd "C:\Users\Yao Thibaut Kpegli\Desktop\AAO\All-at-Once-A-Comprehensive-and-Tractable-Semi-Parametric-Method"  /* need to change the path of the data here */ 

import excel TK92.xlsx, sheet("T&K1992") firstrow


eststo clear
drop if x==.
sort Nature p

tab p, gen(I_)


gen x_mix_four=.
gen y_mix_four=.
gen z_mix_four=.
gen w_mix_four=.

replace x_mix_four=112 in 1
replace x_mix_four=301 in 2
replace y_mix_four=-50 in 1
replace y_mix_four=-125 in 2
replace z_mix_four=50 in 1
replace z_mix_four=150 in 2
replace w_mix_four=-20 in 1
replace w_mix_four=-50 in 2


gen G=(ce>=0)
gen L=(ce<0)


****************************************************************************************************************************************************************
******************************************************************************Semi-parametric*******************************************************************
****************************************************************************************************************************************************************


            ************* Step 0: checking for heteroscedasticity
			
			
program define all_at_once_hetero
args lnf alpha w01 w05 w10 w25 w50 w75 w90 w95 w99 sigma lna

tempvar ce x y I01 I05 I10 I25 I50 I75 I90 I95 I99 Ig Il ux uy w eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `I01' = $ML_y4
generate double `I05' = $ML_y5
generate double `I10' = $ML_y6
generate double `I25' = $ML_y7
generate double `I50' = $ML_y8
generate double `I75' = $ML_y9
generate double `I90' = $ML_y10
generate double `I95' = $ML_y11
generate double `I99' = $ML_y12
generate double `Ig' = $ML_y13
generate double `Il' = $ML_y14


*****Rmk 1: I constrain utility curvature be positif to improve numerical derivatives
generate double `ux' = ((`Ig'-`Il')*`x')^(exp(`alpha'))
generate double `uy' = ((`Ig'-`Il')*`y')^(exp(`alpha'))

*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w'= normal(`w01'*`I01' + `w05'*`I05' +`w10'*`I10'  +`w25'*`I25' +`w50'*`I50'  +`w75'*`I75'  +`w90'*`I90'  +`w95'*`I95'  +`w99'*`I99')

generate double `eror'= `ce'- ((`Ig'-`Il')*((`ux'-`uy')*`w'+`uy')^(1/exp(`alpha')))

replace `lnf' =  ln(normalden(`eror',0,`sigma'*abs(`x'-`y')^exp(`lna')))

}

end

ml model lf all_at_once_hetero (alpha: ce x y I_1 I_2 I_3 I_4 I_5 I_6 I_7 I_8 I_9 G L= G L, nocons) (w01: G L, nocons) (w05: G L, nocons) (w10: G L, nocons) (w25: G L, nocons) (w50: G L, nocons) (w75: G L, nocons) (w90: G L, nocons) (w95: G L, nocons) (w99: G L, nocons) (sigma: G L, nocons) (lna: ) if Nature!="Mixed"

ml maximize

estat ic

************************ Testing heteroscedasticity coming from outcome range*********************

di  "The value of the heteroscedasticity parameter"
nlcom (a: exp(_b[lna:_cons]))

*It turn out that heteroscedasticity from outcomes is rejected at median level

************************ Testing dependance of the variance to the domain
test _b[sigma:G]=_b[sigma:L]
* we accept that variance are equal for both domain at the median level 



nlcom (cur_gain: 1/(1+exp(_b[alpha:G])) )(cur_loss: - 1/(1+exp(_b[alpha:L])) ) (alpha: exp(_b[alpha:G]))(beta: exp(_b[alpha:L])) ///
(w01: normal(_b[w01:G])) (w05: normal(_b[w05:G])) (w10: normal(_b[w10:G])) (w25: normal(_b[w25:G])) (w75: normal(_b[w50:G])) (w90: normal(_b[w75:G])) (w95: normal(_b[w95:G])) (w99: normal(_b[w99:G])) ///
(wl01: normal(_b[w01:L])) (wl05: normal(_b[w05:L])) (wl10: normal(_b[w10:L])) (wl25: normal(_b[w25:L])) (wl75: normal(_b[w50:L])) (wl90: normal(_b[w75:L])) (wl95: normal(_b[w95:L])) (wl99: normal(_b[w99:L]))(sigmaG: _b[sigma:G])(sigmaL: _b[sigma:L])(lna: _b[lna:_cons]) , post 

test _b[cur_gain]=0.5
test _b[cur_loss]=-0.5
test _b[cur_gain]=-_b[cur_loss]




            **************************************************1a/  Power utility specification (no contrain)*****************************************


***First step (utility and proba weighting function) with homoscedasticity (it holds from previous tests: step 0)


eststo clear 

program define all_at_once_homoscedastic
args lnf alpha w01 w05 w10 w25 w50 w75 w90 w95 w99 sigma

tempvar ce x y I01 I05 I10 I25 I50 I75 I90 I95 I99 Ig Il ux uy w eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `I01' = $ML_y4
generate double `I05' = $ML_y5
generate double `I10' = $ML_y6
generate double `I25' = $ML_y7
generate double `I50' = $ML_y8
generate double `I75' = $ML_y9
generate double `I90' = $ML_y10
generate double `I95' = $ML_y11
generate double `I99' = $ML_y12
generate double `Ig' = $ML_y13
generate double `Il' = $ML_y14


*****Rmk 1: I constrain utility curvature be positif to improve numerical derivatives
generate double `ux' = ((`Ig'-`Il')*`x')^(exp(`alpha'))
generate double `uy' = ((`Ig'-`Il')*`y')^(exp(`alpha'))
*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w'= normal(`w01'*`I01' + `w05'*`I05' +`w10'*`I10'  +`w25'*`I25' +`w50'*`I50'  +`w75'*`I75'  +`w90'*`I90'  +`w95'*`I95'  +`w99'*`I99')

generate double `eror'= `ce'- ((`Ig'-`Il')*((`ux'-`uy')*`w'+`uy')^(1/exp(`alpha')))

replace `lnf' = -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end


ml model lf all_at_once_homoscedastic (alpha: ce x y I_1 I_2 I_3 I_4 I_5 I_6 I_7 I_8 I_9 G L= G L, nocons) (w01: G L, nocons) (w05: G L, nocons) (w10: G L, nocons) (w25: G L, nocons) (w50: G L, nocons) (w75: G L, nocons) (w90: G L, nocons) (w95: G L, nocons) (w99: G L, nocons) (sigma: ) if Nature!="Mixed" 

ml maximize

estat ic


********Result MLE
eststo : nlcom (alphp: exp(_b[alpha:G]))(betap: exp(_b[alpha:L])) ///
(w1: normal(_b[w01:G])) (w5: normal(_b[w05:G])) (w10: normal(_b[w10:G])) (w25: normal(_b[w25:G])) (w50: normal(_b[w50:G])) (w75: normal(_b[w75:G])) (w90: normal(_b[w90:G]))(w95: normal(_b[w95:G])) (w99: normal(_b[w99:G])) ///
(wl1: normal(_b[w01:L])) (wl5: normal(_b[w05:L])) (wl10: normal(_b[w10:L])) (wl25: normal(_b[w25:L])) (wl50: normal(_b[w50:L])) (wl75: normal(_b[w75:L])) (wl90: normal(_b[w90:L])) (wl95: normal(_b[w95:L])) (wl99: normal(_b[w99:L])) (sigma: _b[sigma:_cons]) , post


****Test: Partial reflection

test _b[alphp]=_b[betap]

**** Test: linear utility function
test _b[alphp]=1 
test _b[betap]=1 

**** Test of linearity of PW

*gain
test _b[w1]=0.01
test _b[w5]=0.05
test _b[w10]=0.1
test _b[w25]=0.25
test _b[w50]=0.5
test _b[w75]=0.75
test _b[w90]=0.9
test _b[w95]=0.95
test _b[w99]=0.99

*** Loss
test _b[wl1]=0.01
test _b[wl5]=0.05
test _b[wl10]=0.1
test _b[wl25]=0.25
test _b[wl50]=0.5
test _b[wl75]=0.75
test _b[wl90]=0.9
test _b[wl95]=0.95
test _b[wl99]=0.99


*gain
test (_b[w1]=0.01)(_b[w5]=0.05)(_b[w10]=0.1)(_b[w25]=0.25)(_b[w50]=0.5)(_b[w75]=0.75)(_b[w90]=0.9)(_b[w95]=0.95)(_b[w99]=0.99)

*loss
test (_b[wl1]=0.01)(_b[wl5]=0.05)(_b[wl10]=0.1)(_b[wl25]=0.25)(_b[wl50]=0.5)(_b[wl75]=0.75)(_b[wl90]=0.9)(_b[wl95]=0.95)(_b[wl99]=0.99)
*all
test (_b[w1]=0.01)(_b[w5]=0.05)(_b[w10]=0.1)(_b[w25]=0.25)(_b[w50]=0.5)(_b[w75]=0.75)(_b[w90]=0.9)(_b[w95]=0.95)(_b[w99]=0.99)(_b[wl1]=0.01)(_b[wl5]=0.05)(_b[wl10]=0.1)(_b[wl25]=0.25)(_b[wl50]=0.5)(_b[wl75]=0.75)(_b[wl90]=0.9)(_b[wl95]=0.95)(_b[wl99]=0.99)

**** test of homogeneity in probability weight
test _b[w1]=_b[wl1]
test _b[w5]=_b[wl5]
test _b[w10]=_b[wl10]
test _b[w25]=_b[wl25]
test _b[w50]=_b[wl50]
test _b[w75]=_b[wl75]
test _b[w90]=_b[wl90]
test _b[w95]=_b[wl95]
test _b[w99]=_b[wl99]

test (_b[w1]=_b[wl1])(_b[w5]=_b[wl5])(_b[w10]=_b[wl10])(_b[w25]=_b[wl25])(_b[w50]=_b[wl50])(_b[w75]=_b[wl75])(_b[w90]=_b[wl90])(_b[w95]=_b[wl95])(_b[w99]=_b[wl99])

**** test of duality
test _b[w1]+_b[wl99]=1
test _b[w5]+_b[wl95]=1
test _b[w10]+_b[wl90]=1
test _b[w25]+_b[wl75]=1
test _b[w50]+_b[wl50]=1
test _b[w75]+_b[wl25]=1
test _b[w90]+_b[wl10]=1
test _b[w95]+_b[wl5]=1
test _b[w99]+_b[wl1]=1

test (_b[w1]+_b[wl99]=1)(_b[w5]+_b[wl95]=1)(_b[w10]+_b[wl90]=1)(_b[w25]+_b[wl75]=1)(_b[w50]+_b[wl50]=1)(_b[w75]+_b[wl25]=1)(_b[w90]+_b[wl10]=1)(_b[w95]+_b[wl5]=1)(_b[w99]+_b[wl1]=1)


*perfect reflection test

test (_b[alphp]=_b[betap])(_b[w1]=_b[wl1])(_b[w5]=_b[wl5])(_b[w10]=_b[wl10])(_b[w25]=_b[wl25])(_b[w50]=_b[wl50])(_b[w75]=_b[wl75])(_b[w90]=_b[wl90])(_b[w95]=_b[wl95])(_b[w99]=_b[wl99])


 sort proba_fi
 gen wpower_g=.
 replace wpower_g=_b[w1]    in 1
 replace wpower_g=_b[w5]    in 2
 replace wpower_g=_b[w10]   in 3
 replace wpower_g=_b[w25]   in 4
 replace wpower_g=_b[w50]   in 5
 replace wpower_g=_b[w75]   in 6
 replace wpower_g=_b[w90]   in 7
 replace wpower_g=_b[w95]   in 8
 replace wpower_g=_b[w99]   in 9
 
  
    sca wgainp=_b[w50]
	sca alphp=_b[alphp]

	 gen wpower_l=.
 replace wpower_l=_b[wl1]    in 1
 replace wpower_l=_b[wl5]    in 2
 replace wpower_l=_b[wl10]   in 3
 replace wpower_l=_b[wl25]   in 4
 replace wpower_l=_b[wl50]   in 5
 replace wpower_l=_b[wl75]   in 6
 replace wpower_l=_b[wl90]   in 7
 replace wpower_l=_b[wl95]   in 8
 replace wpower_l=_b[wl99]   in 9
  
    sca wlossp=_b[wl50]
	sca betap=_b[betap]

	

*Second step: Loss aversion estimation
 
 
gen la_power_nc=((x^(scalar(alphp)))*scalar(wgainp))/(scalar(wlossp)*(-y)^(scalar(betap))) if Nature=="Mixed" 
gen la_four_power_nc=(( z_mix_four^(scalar(alphp)) -  x_mix_four^(scalar(alphp)))*scalar(wgainp))/(scalar(wlossp)*( (-w_mix_four)^(scalar(betap)) - (- y_mix_four)^(scalar(betap)) )) 
 

 ***1b/  Power utility specification (constrain alpha=beta)
 
 
ml model lf all_at_once_homoscedastic (alpha: ce x y I_1 I_2 I_3 I_4 I_5 I_6 I_7 I_8 I_9 G L= ) (w01: G L, nocons) (w05: G L, nocons) (w10: G L, nocons) (w25: G L, nocons) (w50: G L, nocons) (w75: G L, nocons) (w90: G L, nocons) (w95: G L, nocons) (w99: G L, nocons) (sigma: ) if Nature!="Mixed" 

ml maximize

estat ic 

eststo : nlcom (alphp: exp(_b[alpha:_cons])) ///
(w1: normal(_b[w01:G])) (w5: normal(_b[w05:G])) (w10: normal(_b[w10:G])) (w25: normal(_b[w25:G])) (w50: normal(_b[w50:G])) (w75: normal(_b[w75:G])) (w90: normal(_b[w90:G]))(w95: normal(_b[w95:G])) (w99: normal(_b[w99:G])) ///
(wl1: normal(_b[w01:L])) (wl5: normal(_b[w05:L])) (wl10: normal(_b[w10:L])) (wl25: normal(_b[w25:L])) (wl50: normal(_b[w50:L])) (wl75: normal(_b[w75:L])) (wl90: normal(_b[w90:L])) (wl95: normal(_b[w95:L])) (wl99: normal(_b[w99:L])) (sigma: _b[sigma:_cons]) , post





*** Test : linear u

test _b[alphp]=1


***** Test of linearity PW
*gain
test _b[w1]=0.01
test _b[w5]=0.05
test _b[w10]=0.1
test _b[w25]=0.25
test _b[w50]=0.5
test _b[w75]=0.75
test _b[w90]=0.9
test _b[w95]=0.95
test _b[w99]=0.99

*** Loss
test _b[wl1]=0.01
test _b[wl5]=0.05
test _b[wl10]=0.1
test _b[wl25]=0.25
test _b[wl50]=0.5
test _b[wl75]=0.75
test _b[wl90]=0.9
test _b[wl95]=0.95
test _b[wl99]=0.99


 sort proba_fi
 gen wpower_gc=.
 replace wpower_gc=_b[w1]    in 1
 replace wpower_gc=_b[w5]    in 2
 replace wpower_gc=_b[w10]   in 3
 replace wpower_gc=_b[w25]   in 4
 replace wpower_gc=_b[w50]   in 5
 replace wpower_gc=_b[w75]   in 6
 replace wpower_gc=_b[w90]   in 7
 replace wpower_gc=_b[w95]   in 8
 replace wpower_gc=_b[w99]   in 9
 
  
    sca wgainpc=_b[w50]
	sca alphpc=_b[alphp]

 gen wpower_lc=.
 replace wpower_lc=_b[wl1]    in 1
 replace wpower_lc=_b[wl5]    in 2
 replace wpower_lc=_b[wl10]   in 3
 replace wpower_lc=_b[wl25]   in 4
 replace wpower_lc=_b[wl50]   in 5
 replace wpower_lc=_b[wl75]   in 6
 replace wpower_lc=_b[wl90]   in 7
 replace wpower_lc=_b[wl95]   in 8
 replace wpower_lc=_b[wl99]   in 9
  
 sca wlosspc=_b[wl50]
 
 * Second step: Loss aversion estimation
 
 gen la_power_c=((x^(scalar(alphpc)))*scalar(wgainpc))/(scalar(wlosspc)*(-y)^(scalar(alphpc))) if Nature=="Mixed" 

  gen la_four_power_c=(( z_mix_four^(scalar(alphpc)) -  x_mix_four^(scalar(alphpc)))*scalar(wgainpc))/(scalar(wlosspc)*( (-w_mix_four)^(scalar(alphpc)) - (- y_mix_four)^(scalar(alphpc)) )) 

  

  
 ******1c/ unconstrained power + identical proba weighting
 
 
ml model lf all_at_once_homoscedastic (alpha: ce x y I_1 I_2 I_3 I_4 I_5 I_6 I_7 I_8 I_9 G L= G L, nocons) (w01: ) (w05: ) (w10: ) (w25: ) (w50: ) (w75: ) (w90: ) (w95: ) (w99: ) (sigma: ) if Nature!="Mixed" 

ml maximize

estat ic 

eststo : nlcom (alphp: exp(_b[alpha:G]))(betap: exp(_b[alpha:L])) ///
(w1: normal(_b[w01:_cons])) (w5: normal(_b[w05:_cons])) (w10: normal(_b[w10:_cons])) (w25: normal(_b[w25:_cons])) (w50: normal(_b[w50:_cons])) (w75: normal(_b[w75:_cons])) (w90: normal(_b[w90:_cons]))(w95: normal(_b[w95:_cons])) (w99: normal(_b[w99:_cons])) (sigma: _b[sigma:_cons]) , post


 ******1d/ unconstrained power +  proba weighting under duality assumption
 
program define all_at_once_dual
args lnf ap bp w01 w05 w10 w25 w50 w75 w90 w95 w99 sigma

tempvar ce x y Ig Il I01 I05 I10 I25 I50 I75 I90 I95 I99 ux uy w eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `Ig' = $ML_y4
generate double `Il' = $ML_y5
generate double `I01' = $ML_y6
generate double `I05' = $ML_y7
generate double `I10' = $ML_y8
generate double `I25' = $ML_y9
generate double `I50' = $ML_y10
generate double `I75' = $ML_y11
generate double `I90' = $ML_y12
generate double `I95' = $ML_y13
generate double `I99' = $ML_y14


*****Rmk 1: I constrain utility curvature be positif to improve numerical derivatives
generate double `ux' = ((`Ig'-`Il')*`x')^(exp(`ap'*`Ig'+ `bp'*`Il'))
generate double `uy' = ((`Ig'-`Il')*`y')^(exp(`ap'*`Ig'+ `bp'*`Il'))
*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w' = `Ig'*normal(`w01'*`I01'+ `w05'*`I05'+ `w10'*`I10'+ `w25'*`I25'+ `w50'*`I50'+ `w75'*`I75'+ `w90'*`I90'+  `w95'*`I95'+ `w99'*`I99') + `Il'*(1-normal(`w99'*`I01' + `w95'*`I05'+ `w90'*`I10'+ `w75'*`I25'+ `w50'*`I50'+ `w25'*`I75' + `w10'*`I90'+ `w05'*`I95' + `w01'*`I99'))

generate double `eror'= `ce'- ((`Ig'-`Il')*((`ux'-`uy')*`w'+`uy')^(1/exp(`ap'*`Ig'+ `bp'*`Il')))

replace `lnf' = -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end

ml model lf all_at_once_dual (ap: ce x y G L I_1 I_2 I_3 I_4 I_5 I_6 I_7 I_8 I_9=  ) (bp:  ) (w01: )(w05: )(w10: )(w25: )(w50: )(w75: )(w90: )(w95: )(w99: ) (sigma: ) if Nature!="Mixed" 

ml maximize

estat ic

eststo : nlcom (alphp: exp(_b[ap:_cons]))(betap: exp(_b[bp:_cons])) ///
(w1: normal(_b[w01:_cons])) (w5: normal(_b[w05:_cons])) (w10: normal(_b[w10:_cons])) (w25: normal(_b[w25:_cons])) (w50: normal(_b[w50:_cons])) (w75: normal(_b[w75:_cons])) (w90: normal(_b[w90:_cons]))(w95: normal(_b[w95:_cons])) (w99: normal(_b[w99:_cons])) (sigma: _b[sigma:_cons]) , post



 ***2/  Exponential utility specification 

*2a/ Expo sans aucune contrainte sur la proba weight funct.


clear programs

program define  all_at_once_exp_homoscedastic
args lnf ae be w01 w05 w10 w25 w50 w75 w90 w95 w99 sigma

tempvar ce x y I01 I05 I10 I25 I50 I75 I90 I95 I99 Ig Il w  cept eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `I01' = $ML_y4
generate double `I05' = $ML_y5
generate double `I10' = $ML_y6
generate double `I25' = $ML_y7
generate double `I50' = $ML_y8
generate double `I75' = $ML_y9
generate double `I90' = $ML_y10
generate double `I95' = $ML_y11
generate double `I99' = $ML_y12
generate double `Ig' = $ML_y13
generate double `Il' = $ML_y14


*****Rmk 1: I constrain utility curvature be positif to improve numerical derivatives

generate double `w'= (`w01'*`I01' + `w05'*`I05' +`w10'*`I10'  +`w25'*`I25' +`w50'*`I50'  +`w75'*`I75'  +`w90'*`I90'  +`w95'*`I95'  +`w99'*`I99')


generate double `cept' =  ln((exp((`be'*`Il' - `ae'*`Ig')*`x')-exp((`be'*`Il' - `ae'*`Ig')*`y'))*`w'+ exp((`be'*`Il' - `ae'*`Ig')*`y'))/(`be'*`Il' - `ae'*`Ig')

generate double `eror'= `ce'-`cept'

replace `lnf' = -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end


ml model lf all_at_once_exp_homoscedastic (ae: ce x y I_1 I_2 I_3 I_4 I_5 I_6 I_7 I_8 I_9 G L= ) (be: ) (w01: G L, nocons) (w05: G L, nocons) (w10: G L, nocons) (w25: G L, nocons) (w50: G L, nocons) (w75: G L, nocons) (w90: G L, nocons) (w95: G L, nocons) (w99: G L, nocons) (sigma: ) if Nature!="Mixed" 

ml maximize

estat ic

eststo : nlcom (betae: _b[be: _cons])(alphe: _b[ae: _cons]) ///
(w1: _b[w01:G]) (w5: _b[w05:G]) (w10: _b[w10:G]) (w25: _b[w25:G]) (w50: _b[w50:G]) (w75: _b[w75:G]) (w90: _b[w90:G])(w95: _b[w95:G]) (w99: _b[w99:G]) ///
(wl1: _b[w01:L]) (wl5: _b[w05:L]) (wl10: _b[w10:L]) (wl25: _b[w25:L]) (wl50: _b[w50:L]) (wl75: _b[w75:L]) (wl90: _b[w90:L]) (wl95: _b[w95:L]) (wl99: _b[w99:L]) (sigma: _b[sigma:_cons]) , post

**** Test : partial reflection

test _b[alphe]=_b[betae]

**** Test: linear u
test _b[alphe]=0
test _b[betae]=0
test _b[alphe]=_b[betae]=0


    sort proba_fi
 gen wexpo_g=.
 replace wexpo_g=_b[w1]    in 1
 replace wexpo_g=_b[w5]    in 2
 replace wexpo_g=_b[w10]   in 3
 replace wexpo_g=_b[w25]   in 4
 replace wexpo_g=_b[w50]   in 5
 replace wexpo_g=_b[w75]   in 6
 replace wexpo_g=_b[w90]   in 7
 replace wexpo_g=_b[w95]   in 8
 replace wexpo_g=_b[w99]   in 9
  
  sca wgaine=_b[w50]
  sca alphe=_b[alphe]
  
  
 gen wexpo_l=.
 replace wexpo_l=_b[wl1]    in 1
 replace wexpo_l=_b[wl5]    in 2
 replace wexpo_l=_b[wl10]   in 3
 replace wexpo_l=_b[wl25]   in 4
 replace wexpo_l=_b[wl50]   in 5
 replace wexpo_l=_b[wl75]   in 6
 replace wexpo_l=_b[wl90]   in 7
 replace wexpo_l=_b[wl95]   in 8
 replace wexpo_l=_b[wl99]   in 9
 
  
  sca wlosse=_b[wl50]
  sca betae=_b[betae]
   


   
***  Test of linearity PW
   
   *gain
test _b[w1]=0.01
test _b[w5]=0.05
test _b[w10]=0.1
test _b[w25]=0.25
test _b[w50]=0.5
test _b[w75]=0.75
test _b[w90]=0.9
test _b[w95]=0.95
test _b[w99]=0.99

*** Loss
test _b[wl1]=0.01
test _b[wl5]=0.05
test _b[wl10]=0.1
test _b[wl25]=0.25
test _b[wl50]=0.5
test _b[wl75]=0.75
test _b[wl90]=0.9
test _b[wl95]=0.95
test _b[wl99]=0.99

test (_b[w1]=0.01)(_b[w5]=0.05)(_b[w10]=0.1)(_b[w25]=0.25)(_b[w50]=0.5)(_b[w75]=0.75)(_b[w90]=0.9)(_b[w90]=0.9)(_b[w95]=0.95)(_b[w99]=0.99)(_b[wl1]=0.01)(_b[wl5]=0.05)(_b[wl10]=0.1)(_b[wl25]=0.25)(_b[wl50]=0.5)(_b[wl75]=0.75)(_b[wl90]=0.9)(_b[wl90]=0.9)(_b[wl95]=0.95)(_b[wl99]=0.99)



**** test of homogeneity in probability weight
test _b[w1]=_b[wl1]
test _b[w5]=_b[wl5]
test _b[w10]=_b[wl10]
test _b[w25]=_b[wl25]
test _b[w50]=_b[wl50]
test _b[w75]=_b[wl75]
test _b[w90]=_b[wl90]
test _b[w95]=_b[wl95]
test _b[w99]=_b[wl99]

test (_b[w1]=_b[wl1])(_b[w5]=_b[wl5])(_b[w10]=_b[wl10])(_b[w25]=_b[wl25])(_b[w50]=_b[wl50])(_b[w75]=_b[wl75])(_b[w90]=_b[wl90])(_b[w95]=_b[wl95])(_b[w99]=_b[wl99])

**** test of duality
test _b[w1]+_b[wl99]=1
test _b[w5]+_b[wl95]=1
test _b[w10]+_b[wl90]=1
test _b[w25]+_b[wl75]=1
test _b[w50]+_b[wl50]=1
test _b[w75]+_b[wl25]=1
test _b[w90]+_b[wl10]=1
test _b[w95]+_b[wl5]=1
test _b[w99]+_b[wl1]=1

test (_b[w1]+_b[wl99]=1)(_b[w5]+_b[wl95]=1)(_b[w10]+_b[wl90]=1)(_b[w25]+_b[wl75]=1)(_b[w50]+_b[wl50]=1)(_b[w75]+_b[wl25]=1)(_b[w90]+_b[wl10]=1)(_b[w95]+_b[wl5]=1)(_b[w99]+_b[wl1]=1)


*perfect reflection test
test (_b[alphe]=_b[betae])(_b[w1]=_b[wl1])(_b[w5]=_b[wl5])(_b[w10]=_b[wl10])(_b[w25]=_b[wl25])(_b[w50]=_b[wl50])(_b[w75]=_b[wl75])(_b[w90]=_b[wl90])(_b[w95]=_b[wl95])(_b[w99]=_b[wl99])


esttab using result_spm.tex, wide compress r2(4) ar2(4) scalars(F df_m df_r) se star(* 0.1 ** 0.05 *** 0.01) replace

 
 * Loss aversion estimation
 
gen la_expo=((scalar(wgaine))/(scalar(wlosse)))*(scalar(betae)/scalar(alphe))*((1-exp(-scalar(alphe)*x))/(1-exp(scalar(betae)*y))) if Nature=="Mixed" 

gen la_four_expo=((scalar(wgaine))/(scalar(wlosse)))*(scalar(betae)/scalar(alphe))*((exp(-scalar(alphe)*x_mix_four)-exp(-scalar(alphe)*z_mix_four)  )/(exp(scalar(betae)*y_mix_four) - exp(scalar(betae)*w_mix_four) )) 

li la_power_nc  la_power_c  la_expo  if Nature=="Mixed" 
li  la_four_power_nc la_four_power_c la_four_expo if la_four_power_c !=.

sum la_power_nc la_power_c la_expo if Nature=="Mixed" 

 
 *2b/ Expo + identical probab weight funct.


ml model lf all_at_once_exp_homoscedastic (ae: ce x y I_1 I_2 I_3 I_4 I_5 I_6 I_7 I_8 I_9 G L= ) (be: ) (w01: ) (w05: ) (w10: ) (w25: ) (w50: ) (w75: ) (w90: ) (w95: ) (w99: ) (sigma: ) if Nature!="Mixed"

ml maximize

estat ic

eststo : nlcom (betae: _b[be: _cons])(alphe: _b[ae: _cons]) ///
(w1: _b[w01:_cons]) (w5: _b[w05:_cons]) (w10: _b[w10:_cons]) (w25: _b[w25:_cons]) (w50: _b[w50:_cons]) (w75: _b[w75:_cons]) (w90: _b[w90:_cons])(w95: _b[w95:_cons]) (w99: _b[w99:_cons]) ///
(sigma: _b[sigma:_cons]) , post

esttab using result_spm_iden.tex, wide compress r2(4) ar2(4)  se star(* 0.1 ** 0.05 *** 0.01) replace

 *2c/ Expo + duality assumption

program define all_at_once_exp_dual
args lnf ae be w01 w05 w10 w25 w50 w75 w90 w95 w99 sigma

tempvar ce x y Ig Il I01 I05 I10 I25 I50 I75 I90 I95 I99 w cept eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `Ig' = $ML_y4
generate double `Il' = $ML_y5
generate double `I01' = $ML_y6
generate double `I05' = $ML_y7
generate double `I10' = $ML_y8
generate double `I25' = $ML_y9
generate double `I50' = $ML_y10
generate double `I75' = $ML_y11
generate double `I90' = $ML_y12
generate double `I95' = $ML_y13
generate double `I99' = $ML_y14

*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w' = `Ig'*(`w01'*`I01'+ `w05'*`I05'+ `w10'*`I10'+ `w25'*`I25'+ `w50'*`I50'+ `w75'*`I75'+ `w90'*`I90'+  `w95'*`I95'+ `w99'*`I99') + `Il'*(1-(`w99'*`I01' + `w95'*`I05'+ `w90'*`I10'+ `w75'*`I25'+ `w50'*`I50'+ `w25'*`I75' + `w10'*`I90'+ `w05'*`I95' + `w01'*`I99'))

generate double `cept' =  ln((exp((`be'*`Il' - `ae'*`Ig')*`x')-exp((`be'*`Il' - `ae'*`Ig')*`y'))*`w'+ exp((`be'*`Il' - `ae'*`Ig')*`y'))/(`be'*`Il' - `ae'*`Ig')

generate double `eror'= `ce'- `cept'

replace `lnf' = -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end

ml model lf all_at_once_exp_dual (ae: ce x y G L I_1 I_2 I_3 I_4 I_5 I_6 I_7 I_8 I_9=  ) (be:  ) (w01: )(w05: )(w10: )(w25: )(w50: )(w75: )(w90: )(w95: )(w99: ) (sigma: ) if Nature!="Mixed" 

ml maximize

estat ic

eststo : nlcom (betae: _b[be: _cons])(alphe: _b[ae: _cons]) ///
(w1: _b[w01:_cons]) (w5: _b[w05:_cons]) (w10: _b[w10:_cons]) (w25: _b[w25:_cons]) (w50: _b[w50:_cons]) (w75: _b[w75:_cons]) (w90: _b[w90:_cons])(w95: _b[w95:_cons]) (w99: _b[w99:_cons]) (sigma: _b[sigma:_cons]) , post


* esttab using result_spm_dual.tex, wide compress r2(4) ar2(4) se star(* 0.1 ** 0.05 *** 0.01) replace


********************************************************************************
************************************Parametric**********************************
********************************************************************************
 
  
*a)********** Power + T&K function (no constrain)
 
* First step:

clear programs

program define all_at_once_para
args lnf alpha a sigma

tempvar ce x y Ig Il p ux uy w eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `Ig' = $ML_y4
generate double `Il' = $ML_y5
generate double `p' = $ML_y6

generate double `ux' = ((`Ig'-`Il')*`x')^`alpha'
generate double `uy' = ((`Ig'-`Il')*`y')^`alpha'
generate double `w'=  (`p'^normal(`a'))/( ( `p'^normal(`a')+(1-`p')^normal(`a') )^(1/normal(`a')) )   

generate double `eror'= `ce'- ((`Ig'-`Il')*((`ux'-`uy')*`w'+`uy')^(1/`alpha'))

replace `lnf' =  -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end

ml model lf all_at_once_para (alpha: ce x y G L p=G L, nocons  ) (a: G L, nocons) (sigma: ) if Nature!="Mixed", maximize technique(nr)
ml di
estat ic

eststo: nlcom (cur_gain: 1/(1+_b[alpha: G]) )(cur_loss: -1/(1+_b[alpha: L]) )(alphp: _b[alpha: G])(betap: _b[alpha: L])(a_g: normal(_b[a:G]))(a_l: normal(_b[a:L])) (sigma: _b[sigma:_cons] ), post

test _b[cur_gain]=0.5
test _b[cur_loss]=-0.5


test _b[alphp]=1
test _b[betap]=1
test _b[alphp]=_b[betap]

test _b[a_g]=1
test _b[a_l]=1
test _b[a_g]=_b[a_l]

  sca alphtk=_b[alphp]
  sca wgaintk= 0.5^(_b[a_g])/((0.5^(_b[a_g])+(1-0.5)^(_b[a_g]))^(1/(_b[a_g])))
  sca betatk=_b[betap]
  sca wlosstk= 0.5^(_b[a_l])/((0.5^(_b[a_l])+(1-0.5)^(_b[a_l]))^(1/(_b[a_l])))
  
  
 
 * Second step: Loss aversion estimation 
 

 gen la_parametric_nc=((x^(scalar(alphtk)))*scalar(wgaintk))/(scalar(wlosstk)*(-y)^(scalar(betatk))) if Nature=="Mixed" 

 gen la_four_para_nc=(( z_mix_four^(scalar(alphtk)) -  x_mix_four^(scalar(alphtk)))*scalar(wgaintk))/(scalar(wlosstk)*( (-w_mix_four)^(scalar(betatk)) - (- y_mix_four)^(scalar(betatk)) ))
 
 
 
 *b)********** Power + T&K function (constrain)
 
 * First step:
 
 
ml model lf all_at_once_para (alpha: ce x y G L p= ) (a: G L, nocons) (sigma: ) if Nature!="Mixed", maximize technique(nr)
ml di
estat ic

eststo : nlcom (cur_gain: 1/(1+_b[alpha:_cons]) )(cur_loss: - 1/(1+_b[alpha:_cons]) )(alphp: _b[alpha:_cons])(a_g: normal(_b[a:G]))(a_l: normal(_b[a:L])) (sigma: _b[sigma: _cons] ), post


test _b[alphp]=1

test _b[a_g]=1
test _b[a_l]=1
test _b[a_g]=_b[a_l]

sca alphtkc=_b[alphp]
sca wgaintkc= 0.5^(_b[a_g])/((0.5^(_b[a_g])+(1-0.5)^(_b[a_g]))^(1/(_b[a_g])))
sca wlosstkc= 0.5^(_b[a_l])/((0.5^(_b[a_l])+(1-0.5)^(_b[a_l]))^(1/(_b[a_l])))

  
 * Second step: Loss aversion estimation 
 
 gen la_parametric_c=((x^(scalar(alphtkc)))*scalar(wgaintkc))/(scalar(wlosstkc)*(-y)^(scalar(alphtkc))) if Nature=="Mixed" 

 gen la_four_para_c=(( z_mix_four^(scalar(alphtkc)) -  x_mix_four^(scalar(alphtkc)))*scalar(wgaintkc))/(scalar(wlosstkc)*( (-w_mix_four)^(scalar(alphtkc)) - (- y_mix_four)^(scalar(alphtkc)) )) 

 
 *esttab using result_paratevlossaversion.tex, wide compress r2(4) ar2(4) scalars(F df_m df_r) se star(* 0.1 ** 0.05 *** 0.01) replace
 
 
 
**c) ********** expo + T&K function 
 

clear programs

program define all_at_once_para_exp
args lnf ae be a sigma

tempvar ce x y Ig Il p w cept eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `Ig' = $ML_y4
generate double `Il' = $ML_y5
generate double `p' = $ML_y6

*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w'=  (`p'^normal(`a'))/( ( `p'^normal(`a')+(1-`p')^normal(`a') )^(1/normal(`a')) )   

generate double `cept' =  ln((exp((`be'*`Il' - `ae'*`Ig')*`x')-exp((`be'*`Il' - `ae'*`Ig')*`y'))*`w'+ exp((`be'*`Il' - `ae'*`Ig')*`y'))/(`be'*`Il' - `ae'*`Ig')

generate double `eror'= `ce'-`cept'


replace `lnf' =  -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end


*eststo : 
ml model lf all_at_once_para_exp (ae: ce x y G L p=  ) (be:  ) (a: G L, nocons) (sigma: ) if Nature!="Mixed"

ml maximize

estat ic


eststo: nlcom (betae: _b[be: _cons]) (alphe: _b[ae: _cons]) (a_g: normal(_b[a:G])) (a_l: normal(_b[a:L])) (sigma: _b[sigma:_cons]) , post


 esttab using result_parate.tex, wide compress r2(4) ar2(4) scalars(F df_m df_r) se star(* 0.1 ** 0.05 *** 0.01) replace
 
 test _b[alphe]=0
 test _b[betae]=0
   
 test _b[alphe]=_b[betae]
  
 
  sca alphe_para=_b[alphe]
  sca wgainetk= 0.5^(_b[a_g])/((0.5^(_b[a_g])+(1-0.5)^(_b[a_g]))^(1/(_b[a_g])))
  sca betae_para=_b[betae]
  sca wlossetk= 0.5^(_b[a_l])/((0.5^(_b[a_l])+(1-0.5)^(_b[a_l]))^(1/(_b[a_l])))
  
  
 *** loss aversion
  
gen la_para_expo=((scalar(wgainetk))/(scalar(wlossetk)))*(scalar(betae_para)/scalar(alphe_para))*((1-exp(-scalar(alphe_para)*x))/(1-exp(scalar(betae_para)*y))) if Nature=="Mixed" 

gen la_four_para_expo=((scalar(wgainetk))/(scalar(wlossetk)))*(scalar(betae_para)/scalar(alphe_para))*((exp(-scalar(alphe_para)*x_mix_four)-exp(-scalar(alphe_para)*z_mix_four)  )/(exp(scalar(betae_para)*y_mix_four) - exp(scalar(betae_para)*w_mix_four) )) 

li la_parametric_nc  la_parametric_c  la_para_expo  if Nature=="Mixed" 
li  la_four_para_nc la_four_para_c la_four_para_expo if la_four_para_c !=.



   * graph SPM & PM
   
   twoway (scatter wpower_g proba_figure , connect(l)) (function y=x^(0.6351778)/((x^(0.6351778)+(1-x)^(0.6351778))^(1/(0.6351778))))  (function y = x) , name(w_gain_sp_p)  legend(label(1 "Semi-para") label(2 "Parametric") label(3 "45° line")) graphregion(color(white))  ytitle(Weighting function w(p)) xtitle(Probability p) title(Gain domain)

  twoway (scatter wpower_l proba_figure , connect(l))(function y=x^(0.7126445)/((x^(0.7126445)+(1-x)^(0.7126445))^(1/(0.7126445))))(function y = x), name(w_loss_sp_p)   legend(label(1 "Semi-para") label(2 "Parametric") label(3 "45° line")) graphregion(color(white))  ytitle(Weighting function w(p)) xtitle(Probability p) title(Loss domain)
 
 
 graph combine w_gain_sp_p w_loss_sp_p, name(comb_3) graphregion(color(white)) 

 
*/  
  


********************************************************************************    
*************************Parametric: Power + prelecI****************************
******************************************************************************** 


clear programs

program define all_at_once_para_prelecI
args lnf alpha a sigma

tempvar ce x y Ig Il p ux uy w  eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `Ig' = $ML_y4
generate double `Il' = $ML_y5
generate double `p' = $ML_y6

*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w'=  exp(-(-ln(p))^`a')


generate double `ux' = ((`Ig'-`Il')*`x')^exp(`alpha')
generate double `uy' = ((`Ig'-`Il')*`y')^exp(`alpha')


generate double `eror'= `ce'- ((`Ig'-`Il')*((`ux'-`uy')*`w'+`uy')^(1/exp(`alpha')))


replace `lnf' =  -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end

 ml model lf all_at_once_para_prelecI (alpha: ce x y G L p=G L, nocons  ) (a: G L, nocons) (sigma: ) if Nature!="Mixed", maximize technique(nr)
ml di
estat ic

/*
nlcom (alphp: exp(_b[alpha: G]))(betap: exp(_b[alpha: L]))(a_g: _b[a:G])(a_l: _b[a:L]) (sigma: _b[sigma: _cons] ), post

outtex, labels level detail legend  key(stab)

tw (function y = exp(-(-ln(x))^(_b[a_l]) )) || function y=x
*/



********************************************************************************    
************************Parametric: Power + prelecII****************************
******************************************************************************** 

 
 
clear programs

program define all_at_once_para_prelecII
args lnf alpha a b sigma

tempvar ce x y Ig Il p ux uy w  eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `Ig' = $ML_y4
generate double `Il' = $ML_y5
generate double `p' = $ML_y6

*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w'=  exp(-`b'*(-ln(p))^`a')


generate double `ux' = ((`Ig'-`Il')*`x')^exp(`alpha')
generate double `uy' = ((`Ig'-`Il')*`y')^exp(`alpha')


generate double `eror'= `ce'- ((`Ig'-`Il')*((`ux'-`uy')*`w'+`uy')^(1/exp(`alpha')))


replace `lnf' =  -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end

 ml model lf all_at_once_para_prelecII (alpha: ce x y G L p=G L, nocons  ) (a: G L, nocons) (b: G L, nocons) (sigma: ) if Nature!="Mixed", maximize technique(nr)
ml di
estat ic

/*
nlcom (alphp: exp(_b[alpha: G]))(betap: exp(_b[alpha: L]))(a_g: _b[a:G])(b_g: _b[b:G])(a_l: _b[a:L])(b_l: _b[b:L]) (sigma: _b[sigma: _cons] ), post

outtex, labels level detail legend  key(stab)

tw (function y = exp(-_b[b_l]*(-ln(x))^(_b[a_l]) )) || function y=x
*/

********************************************************************************    
*************************Parametric: Power + Latimore***************************
********************************************************************************    
 
clear programs

program define all_at_once_para_lat
args lnf alpha a b sigma

tempvar ce x y Ig Il p ux uy w  eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `Ig' = $ML_y4
generate double `Il' = $ML_y5
generate double `p' = $ML_y6

*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w'=  (exp(`b')*`p'^exp(`a'))/( ( exp(`b')*`p'^exp(`a')+(1-`p')^exp(`a') ) )   


generate double `ux' = ((`Ig'-`Il')*`x')^`alpha'
generate double `uy' = ((`Ig'-`Il')*`y')^`alpha'


generate double `eror'= `ce'- ((`Ig'-`Il')*((`ux'-`uy')*`w'+`uy')^(1/`alpha'))


replace `lnf' =  -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end

 ml model lf all_at_once_para_lat (alpha: ce x y G L p=G L, nocons  ) (a: G L, nocons) (b: G L, nocons) (sigma: ) if Nature!="Mixed", maximize technique(nr)
ml di
estat ic


nlcom (alphp: _b[alpha: G])(betap: _b[alpha: L])(a_g: exp(_b[a:G]))(b_g_: exp(_b[b:G]))(a_l: exp(_b[a:L]))(b_l: exp(_b[b:L])) (sigma: _b[sigma: _cons] ), post

outtex, labels level detail legend  key(stab)



********************************************************************************    
************************Parametric: expo + prelecI******************************
******************************************************************************** 
 
 

clear programs

program define all_at_once_para_exp_PrelecI
args lnf ae be a sigma

tempvar ce x y Ig Il p w cept eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `Ig' = $ML_y4
generate double `Il' = $ML_y5
generate double `p' = $ML_y6

*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w'=  exp(-(-ln(p))^`a')  

generate double `cept' =  ln((exp((`be'*`Il' - `ae'*`Ig')*`x')-exp((`be'*`Il' - `ae'*`Ig')*`y'))*`w'+ exp((`be'*`Il' - `ae'*`Ig')*`y'))/(`be'*`Il' - `ae'*`Ig')

generate double `eror'= `ce'-`cept'


replace `lnf' =  -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end


*eststo : 
ml model lf all_at_once_para_exp_PrelecI (ae: ce x y G L p=  ) (be:  ) (a: G L, nocons) (sigma: ) if Nature!="Mixed"

ml maximize

estat ic

/*
eststo: nlcom (betae: _b[be: _cons])(alphe: _b[ae: _cons]) (a_g: _b[a:G])(a_l: _b[a:L]) (sigma: _b[sigma:_cons]), post
*/



********************************************************************************    
************************Parametric: expo + prelecII****************************
******************************************************************************** 

 
 
clear programs


program define all_at_once_para_exp_PrelecII
args lnf ae be a b sigma

tempvar ce x y Ig Il p w cept eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `Ig' = $ML_y4
generate double `Il' = $ML_y5
generate double `p' = $ML_y6

*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w'=  exp(-`b'*(-ln(p))^`a')  

generate double `cept' =  ln((exp((`be'*`Il' - `ae'*`Ig')*`x')-exp((`be'*`Il' - `ae'*`Ig')*`y'))*`w'+ exp((`be'*`Il' - `ae'*`Ig')*`y'))/(`be'*`Il' - `ae'*`Ig')

generate double `eror'= `ce'-`cept'


replace `lnf' =  -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end


ml model lf all_at_once_para_exp_PrelecII (ae: ce x y G L p=  ) (be:  ) (a: G L, nocons)(b: G L, nocons) (sigma: ) if Nature!="Mixed" , maximize technique(nr)
ml di
estat ic

nlcom (alphp: _b[ae:_cons])(betap: _b[be:_cons])(a_g: _b[a:G])(b_g: _b[b:G])(a_l: _b[a:L])(b_l: _b[b:L]) (sigma: _b[sigma: _cons] ), post

outtex, labels level detail legend  key(stab)





********************************************************************************    
************************Parametric: expo + Latimore*****************************
******************************************************************************** 
 

 
  
clear programs


program define all_at_once_para_exp_lat
args lnf ae be a b sigma

tempvar ce x y Ig Il p w cept eror 
quietly {

*** part 1

generate double `ce' = $ML_y1
generate double `x' = $ML_y2
generate double `y' = $ML_y3
generate double `Ig' = $ML_y4
generate double `Il' = $ML_y5
generate double `p' = $ML_y6

*****Rmk 2: I constrain probability weights be in the range (0,1) to improve numerical derivatives
generate double `w'=  (exp(`b')*`p'^exp(`a'))/( ( exp(`b')*`p'^exp(`a')+(1-`p')^exp(`a') ) )   

generate double `cept' =  ln((exp((`be'*`Il' - `ae'*`Ig')*`x')-exp((`be'*`Il' - `ae'*`Ig')*`y'))*`w'+ exp((`be'*`Il' - `ae'*`Ig')*`y'))/(`be'*`Il' - `ae'*`Ig')

generate double `eror'= `ce'-`cept'


replace `lnf' =  -0.5*(`eror'/`sigma')^2 - ln(`sigma') -0.5*ln(2*c(pi))

}

end


ml model lf all_at_once_para_exp_lat (ae: ce x y G L p=  ) (be:  ) (a: G L, nocons)(b: G L, nocons) (sigma: ) if Nature!="Mixed" , maximize technique(nr)
ml di
estat ic

nlcom (alphp: _b[ae:_cons])(betap: _b[ae:_cons])(a_g: exp(_b[a:G]))(b_g: exp(_b[b:G]))(a_l: exp(_b[a:L]))(b_l: exp(_b[b:L])) (sigma: _b[sigma: _cons] ), post

outtex, labels level detail legend  key(stab)



