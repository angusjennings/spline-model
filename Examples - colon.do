ssc install stpm3
ssc install standsurv
ssc install gensplines

clear

********************************************************************************
*	consSurv Paper Example - STATA
*	Angus Jennings 					30sept24
********************************************************************************


********************************************************************************
*Format data
import excel "C:\Users\rmjwaje\OneDrive - University College London\Angus Jennings\4. U o Leicester\Paper 2_V1.1\colon.xlsx", sheet("Sheet1") firstrow clear

keep if etype==2 & inlist(rx, "Obs", "Lev+5FU")

gen yrs = time/365.25

gen rx2 = rx=="Lev+5FU"
drop rx
rename rx2 rx

stset yrs, failure(status==1)


********************************************************************************
*UNCONSTRAINED STPM3
//for demonstration of code, default stpm3 df/knot locations are used, hence predictions differ sliughtly from those presented in paper 'A spline-based approach to smoothly constrain hazard ratios with a view to apply treatment effect waning'
stpm3 i.rx, scale(lnhazard) df(3) tvc(i.rx) dftvc(3)

*Predictions - unconstrained
predict h1 h2, hazard ci ///
	at1(rx 0) ///
	at2(rx 1) ///
	contrast(ratio) ///
	contrastvar(HR) ///
	timevar(0.01 10, step(0.1)) ///
	frame(uncons_h, replace)
	
frame uncons_h: list
	
predict S1, survival ci ///
	at1(rx 0) ///
	timevar(0.01 10, step(0.1)) ///
	frame(uncons_s, replace)
	
predict S2, survival ci ///
	at1(rx 1) ///
	frame(uncons_s, merge)
	
frame uncons_s: list

*Plots - unconstrained
frame uncons_h: twoway line h1 h2 tt, name(haz_uncons, replace) nodraw

frame uncons_s: twoway line S1 S2 tt, name(surv_uncons, replace) nodraw

frame uncons_h: twoway line HR tt if tt>0.3, name(hr_uncons, replace) nodraw


********************************************************************************
*CONSTRAINED STPM3
*ph term
constraint 1 1.rx
*'forward extrapolation' spline term
constraint 2 1.rx#_ns_tvc3

stpm3 i.rx, scale(lnhazard) df(3) tvc(i.rx) dftvc(3) constraints(1 2)

*Predictions - constrained
predict h1 h2, hazard ci ///
	at1(rx 0) ///
	at2(rx 1) ///
	contrast(ratio) ///
	contrastvar(HR) ///
	timevar(0.01 10, step(0.1)) ///
	frame(cons_h, replace)
	
frame cons_h: list
	
predict S1, survival ci ///
	at1(rx 0) ///
	timevar(0.01 10, step(0.1)) ///
	frame(cons_s, replace)
	
predict S2, survival ci ///
	at1(rx 1) ///
	frame(cons_s, merge)
	
frame cons_s: list

*Plots - constrained
frame cons_h: twoway line h1 h2 tt, name(haz_cons, replace) nodraw

frame cons_s: twoway line S1 S2 tt, name(surv_cons, replace) nodraw

frame cons_h: twoway line HR tt if tt>0.3, name(hr_cons, replace) nodraw


********************************************************************************
*FINAL PLOT
gr combine haz_uncons haz_cons surv_uncons surv_cons hr_uncons  hr_cons, cols(2) imargin(0 0 0 0)


