
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
/*
 * --- ---- --- --- --- --- ---- --- --- ------ ---- --- --- --- --- ---
 * Author: Carlos Andres Zapata Q
 * This .mod file uses Dynare's toolkit to reproduce the results of
 * the Thesis. Cap. 3
 * --- ---- --- --- --- --- ---- --- --- ------ ---- --- --- --- --- ---
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * --- ---- --- --- --- --- ---- --- --- ------ ---- --- --- --- --- ---
 */

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

// 1. Declaration of endogenous variables
var g_hat c_hat cr_hat cnr_hat q_hat u_hat	rk_hat I_hat k_hat w_hat N_hat y_hat  
    R_hat tao_c tao_k tao_w kg_hat Ig_hat trans b_hat emp_hat 
    tao_c_inc tao_k_inc tao_w_inc b_hat_y
    pid_hat pix_hat piim_hat deltaE_hat mc_hat mcx_hat mcmi_hat
    pistar_hat ystar_hat Rstar_hat phitilde_hat //
    data_con data_gov data_r data_inv data_pi data_y //data_pet
    pic_hat r_hat //yoil_hat ytot_hat
    e_a  e_i e_r e_pi e_pix e_pim e_n e_w e_tr e_ig e_tc e_tk 
    e_tw e_g  e_dt //
    x_hat im_hat nx_hat
    data_rstar data_pistar data_ystar   
    ; //  mcc_hat prx_hat pri_hat  // data_ex data_imp
   //data_rev data_expend  data_debt n_ex_hat 
    //data_wages data_hours data_taxc data_taxk data_taxw  data_tr data_inv_g n_hat

// ---------------------------------------------------------------------------------------------------
// Real Economic Variables
// g_hat, c_hat, cr_hat, cnr_hat, q_hat, u_hat, rk_hat, I_hat, k_hat, w_hat, N_hat, y_hat
// These variables represent government spending, consumption, Ricardian and non-Ricardian consumption, 
// Tobin's q, capital utilization, rental rate of capital, investment, capital stock, wage rate, labor, 
// and output respectively.

// Policy and Fiscal Variables:
// R_hat, tao_c, tao_k, tao_w, kg_hat, Ig_hat, trans, b_hat, emp_hat
// These include interest rates, tax rates on consumption, capital, labor, public capital, government
// investment, transfers, debt, and employment.

// Inflation and Prices:
// pid_hat, pix_hat, piim_hat, deltaE_hat, mc_hat, mcx_hat, mcmi_hat, pistar_hat, ystar_hat, Rstar_hat, phitilde_hat
// These represent domestic inflation, export and import inflation, exchange rate depreciation, marginal costs for 
// domestic, export, and import goods, foreign inflation, foreign output, foreign interest rate, and a risk premium
// variable.

// Data Variables:
// data_con, data_gov, data_r, data_inv, data_pi, data_y
// These variables likely represent observable data series used for model estimation or calibration.

// ---------------------------------------------------------------------------------------------------

// 2. Declaration of exogenous variables
varexo  
at wt nt tct tkt pit twt it gt igt trt rt pixt pimt e_pistar e_ystar e_Rstar phitildet dt;// yoilt; %

// These represent the stochastic processes driving the exogenous shocks, such as technological, wage, and tax changes.

// ---------------------------------------------------------------------------------------------------

// 3. Declaration of parameters
parameters  rho_p fc omega_w chi_w e_nu chi_h sigma_c  piT chi_x chi_im
            aa beta  alfa delta omega_d gammak kappa rho_R r_pi r_y sigma_n rho_a alfa_g delta_g phi_g phi_ig
            rho_i rho_w std_at std_it std_pit std_nt rho_n std_wt std_phitilde
            yz_bar R_bar c_bar I_bar g_bar Ig_bar tr_bar wN_bar rky_bar b_bar tao_cbar tao_wbar tao_kbar
            std_tc std_tw std_tk std_g std_tr std_ig omega_im rho_r rho_px rho_pm rk_bar
            psi_g psi_ig phi_tr psi_tr phi_tc psi_tc phi_tk psi_tk phi_tw psi_tw omega_h omega_x
            ela_lw ela_lc ela_lk ela_g ela_ig ela_tr etac etai etax chi_e
            rho_tc rho_tk rho_tw rho_g rho_ig rho_tr //rho_oil omega_oil oil_bar 
            phimc phimi pcd_bar pim_bar phitilde phitildes rho_phitilde b_target 
            rho_pistar rho_ystar rho_Rstar piTstar cm_bar im_bar x_bar // phix
            //a111 a112 a113 a121 a122 a123 a131 a132 a133 a211 a212 a213 a221 a222 a223 a231 a232 a233
            //b11 b12 b13 b21 b22 b23 b31 b32 b33
            ;
            // rho_vt  rho_m std_mt eta_l eta_ll eta_lc  nu_bar bp_bar 
            // eta_lk eta_g eta_ig y_bar 
            
// 4. Parameter estimations

etac         =	   1.1896   ;   // elasticity of substition domestic and  foreign consumption goods
etai	     =	   1.1892	;   // elasticity of substition imported consumption and investments goods 
etax	     =	   1.197	;   // elasticity of substition foreign goods and investments
omega_d      =     0.4841   ;   // share of domestic goods (iqual to omega_c)
omega_w      =     0.5423   ;   // share wages
omega_im     =     0.5148   ;   // share of imports (consump. and investment goods)
omega_x      =     0.4989   ;   // share of expots goods
//omega_oil    =     0.7098 ;   // share of oil product. in gov. revenues
chi_h        =     0.4803   ;   // price indexation
chi_w        =     0.5012   ;   // wage indexation
chi_e        =     0.4997   ;   // emp. indexation
chi_x        =     0.5019   ;   // export indexation
chi_im       =     0.4921   ;   // export indexation
phimc        =	   1.1948   ; 
phimi        =	   1.2003   ;
//phix        =	   1.1      ;
phitilde     =	   0.075    ; // Risk premium
phitildes    =	   0.6717   ; // 0.6 ;   
rho_phitilde  =    0.6969   ; 
fc           =     1.3018   ; 

// Foreign Economy block
rho_pistar   =     0.6736 ;
rho_ystar    =     0.9500 ;
rho_Rstar    =     0.6998 ;

// Persistence parameters
rho_a        =     0.7085  ;
rho_n        =     0.8001  ;
rho_i        =     0.4447  ;
rho_w        =     0.6759  ;
rho_p        =     0.519  ;

rho_r        =     0.6   ;
rho_px       =     0.6987  ;
rho_pm       =     0.6831  ;

// Fiscal policy rules - block

// Persistence of Fiscal rules
rho_tc = 0.5297; % 0.485 ;  //  0;
rho_tk = 0.7393; % 0.7564 ; //
rho_tw = 0.7463; % 0.7572 ; //
rho_g  = 0.6742; % 0.6366 ; //
rho_ig = 0.7634 ; % 0.7611 ; //
rho_tr = 0.7699; % 0.770 ;  //

// Responces to GDP
psi_tc =  0.5276  ;
psi_tk =  0.5076  ;
psi_tw =  0.4942  ;
psi_g  =  0.4350  ;
psi_ig =  0.5075  ;
psi_tr =  0.5050  ;


// Responces to debt
ela_lc       =     0.2115  ;
ela_lk       =     0.2072  ;
ela_lw       =     0.2249  ;
ela_g        =     0.1029  ;
ela_ig       =     0.2041  ;
ela_tr       =     0.2126  ;

// Persistence of Shocks
phi_g        =     0.4488  ;
phi_ig       =     0.7762  ;
phi_tr       =     0.7555  ;
phi_tk       =     0.7597  ;
phi_tw       =     0.7502  ;
phi_tc       =     0.6081  ;  

// --------------------------------------------------
//Calibrated steady state values -29

// General parameters
beta      = 0.9922 ;   // discount factor (6% real annual interest rate - COP TES: 1.06^(-1/4))
delta     = 0.025  ;   // depreciation rate of private capital
delta_g   = 0.015  ;   // depreciation rate of public capital
e_nu      = 1.05   ;    // wage markup (nu-1/nu)
alfa      = 0.3    ;    // Share of capital in production
alfa_g    = 0.01   ;    // Elasticity of output to public capital
gammak    = 6.4    ;    // investment adjustment cost parameter
omega_h   = 0.3    ;    // share non-Ricardian households.
aa        = 0.049  ;    // capital utilisation parameter // 0.049
kappa     = 0.7016 ;    // habit formation
sigma_n   = 1      ;    // Labor supply elasticity (is common in literature).
sigma_c   = 1.2    ;    // coefficient of relative risk aversion
piT       = 1.0074 ;    // domestic Inflation target at 3%
piTstar   = 1.005  ;    // external Inflation target at 2%


// Monetary policy parameters - Based on literature
rho_R        =     0.87   ;  // 0.7444  0.87 ; (fixed following to Christiano et al. (2010)
r_pi         =     1.5685 ;  // 1.1  ;
r_y          =     0.1485;   //0.01 ;   

// Fiscal policy block
// Taxes Based on Rincón-Castro and Delgado-Rojas 2018
tao_cbar  = 0.107;  // steady state labor tax rate % 10.7 %0.07; 
tao_kbar  = 0.186;  // steady state capital tax rate % 18.6 %0.1 
tao_wbar  = 0.154;  // steady state labor tax rate % 15.4 % 0.07

// Long-run relations 
g_bar     = 0.16;    // steady state government consumptio
Ig_bar    = 0.02;     // steady state government investment
b_bar     = 1.8 ;    //  steady state gov. debt (mean 39% annual)
b_target  = 1.8 ;    // steady state government debt target - to Y ratio  (MFMP 23: 55% ancla LP) (2)
yz_bar    = 1.0089;   // Delta GDP: 3.6%
//tr_bar    = 0.025;

// //Steady state values
R_bar    = 1/beta; %1/((1-tao_kbar)*beta); 
rky_bar  = alfa;
wN_bar   = 1-alfa;
rk_bar   = (1/(1-tao_kbar))*(1/beta-(1-delta));
//rk_bar = alfa/(1-alfa)*R_bar/rky_bar*wN_bar;
I_bar    = delta*rky_bar/rk_bar;
c_bar    = 1-g_bar-I_bar-Ig_bar;
%y_bar    = yz_bar^(-alfa)*rk_bar^alfa*wN_bar^(1-alfa);
pcd_bar  = ((1-omega_d) + omega_d *(phimc)^(1-etac))^(1/(1-etac));
pim_bar  = ( omega_im + (1-omega_im)*(1/phimi)^(1-etai))^(1/(1-etac));
cm_bar   = omega_d*pcd_bar^(etac)*c_bar;
im_bar   = omega_im*pim_bar^(etai)*I_bar;
x_bar    = cm_bar+im_bar;
tr_bar   = tao_cbar*c_bar+tao_wbar*wN_bar+tao_kbar*rky_bar-g_bar-Ig_bar-R_bar*b_bar+b_bar;

// Socks - standard deviations
std_at       =     0.0118  ;
std_nt       =     0.0111  ;
std_pit      =     0.0101  ;
std_tc       =     0.0232  ;
std_tk       =     0.0105  ;
std_tw       =     0.0123  ;
std_it       =     0.0735  ;
std_tr       =     0.0130  ;
std_ig       =     0.0113  ;
std_wt       =     0.012   ;
std_g        =     0.0248  ;
std_rt       =     0.01    ;
std_phitilde =     0.0125  ;
std_pixt     =     0.0113  ;
std_pimt     =     0.0269  ;
std_dt       =     0.01     ;
//std_oilt    =     0.1312  ;

//---------------------------------------------------------------------------------------------------------------

// Model equations
model(linear);

//ENDOGENOUS VARIABLES

// 1. consumption of non-Ricardian households
(1+tao_cbar)*c_bar*(cnr_hat+(tao_cbar/(1+tao_cbar))*tao_c)=wN_bar*((1-tao_wbar)*(w_hat+N_hat)-tao_wbar*tao_w)+tr_bar*trans;

// 2. consumption of Ricardian households
cr_hat = ((1/(1+kappa))*cr_hat(+1)) + ((kappa/(1+kappa))*cr_hat(-1)) -((1-kappa)/((1+kappa)*sigma_c))*(R_hat)+((1-kappa)/((1+kappa)*sigma_c))*pic_hat(+1)+((1-kappa)/((1+kappa)*sigma_c))*(e_n-e_n(+1))-((1-kappa)/((1+kappa)*sigma_c))*(tao_cbar/(1+tao_cbar))*(tao_c-tao_c(+1));

// First order conditions of Ricardian households
// 3. Capital 
q_hat=-R_hat+pid_hat(+1)+((1-delta)/(1-delta+(1-tao_kbar)*rk_bar))*q_hat(+1)+((1-tao_kbar)*rk_bar/(1-delta+(1-tao_kbar)*rk_bar))*(rk_hat(+1)-(tao_kbar/(1-tao_kbar))*tao_k(+1));

// 4. Investment
I_hat=((1/(gammak*(1+beta)))*(q_hat)+(1/(1+beta))*I_hat(-1)+(beta/(1+beta))*I_hat(+1)-(1/((1+beta)))*(e_i-beta*e_i(+1)));

// 5. Capital utilisation
u_hat=(1/aa)*(rk_hat-(tao_kbar/(1-tao_kbar))*tao_k);

// 6-8. Wages, hours and employment
w_hat=(beta/(1+beta))*w_hat(+1)+(1/(1+beta))*w_hat(-1)+(beta/(1+beta))*pid_hat(+1)-((1+beta*chi_w)/(1+beta))*pid_hat+(chi_w/(1+beta))*pid_hat(-1)-((((1-omega_w)*(1-beta*omega_w))/((1+((1+e_nu)/(e_nu))*sigma_n)*omega_w))*(1/(1+beta)))*(w_hat-sigma_n*N_hat-(sigma_c/(1-kappa))*(cr_hat-kappa*cr_hat(-1))-(tao_wbar/(1-tao_wbar))*tao_w-(tao_cbar/(1+tao_cbar))*tao_c+e_w);
emp_hat = beta/(1+beta)*emp_hat(+1)+(1/(1+beta))*emp_hat(-1)+((1-chi_e)*(1-beta*chi_e))/((1+beta)*chi_e)*(N_hat-emp_hat); 

// 9. Private capital accumulation 
k_hat=((1-delta)*k_hat(-1)+delta*(I_hat));

// Public capital accumulation 
kg_hat=(1-delta_g)*kg_hat(-1) + delta_g*Ig_hat;

// 11. Production
y_hat=fc*(e_a+alfa*k_hat(-1)+alfa*u_hat+(1-alfa)*N_hat+alfa_g*kg_hat(-1));

// 12-13. combination of First order conditions of Firms
rk_hat+u_hat+k_hat(-1)=w_hat+N_hat;   // + Rf_hat

// 14. Exchange rate’s modified UIP condition
(1-phitildes)*deltaE_hat(+1) - phitildes*deltaE_hat - (R_hat - Rstar_hat) - phitilde*u_hat + phitilde_hat = 0;
//deltaE_hat = ((1-phitildes)*deltaE_hat(+1) - (R_hat - Rstar_hat) - phitilde*u_hat + phitilde_hat)/phitildes;

// 15. New-Keynesian Philips curve for domestic producers
pid_hat = piT+(beta/(1+beta*chi_h))*(pid_hat(+1)-piT)+(chi_h/(1+beta*chi_h))*(pid_hat(-1)-piT)
            + ((((1-omega_d)*(1-beta*omega_d))/(omega_d))*(1/(1+beta*chi_h)))*( mc_hat ) +e_pi;

// 16. New-Keynesian Phillips curve for exporters
pix_hat = piTstar + (beta/(1+chi_x*beta))*(pix_hat(+1)-piTstar) +  (chi_x/(1+chi_x*beta))*(pix_hat(-1)-piTstar) 
            + (((1-omega_x)*(1-beta*omega_x))/(omega_x*(1+chi_x*beta)))*( mcx_hat) + e_pix ;

// 17. New-Keynesian Phillips curve for importers
piim_hat = piT+(beta/(1+chi_im*beta))*(piim_hat(+1)-piT) + (chi_im/(1+chi_im*beta))*(piim_hat(-1)-piT) 
            + (((1-omega_im)*(1-beta*omega_im))/(omega_im*(1+chi_im*beta)))*( mcmi_hat) + e_pim; 

// 18. Resource constraint and net assset position 
y_hat=(c_bar)*c_hat+(I_bar)*I_hat+(1-tao_kbar)*(rky_bar)*u_hat+g_bar*g_hat+Ig_bar*Ig_hat; 

// 19. total consumption
c_hat=(1-omega_h)*cr_hat+omega_h*cnr_hat; 

// 20-21. monetary policy rule:
R_hat = rho_R*R_hat(-1) + (1-rho_R)*(r_pi*(piT+(pid_hat(-1)-piT)) + r_y*(y_hat(-1)))+ rt; //ytot_hat

// Fisher Equation
r_hat= R_hat-pic_hat(+1); 

// 22. Marginal cost (domestic)
mc_hat= (1-alfa)*w_hat+alfa*rk_hat-alfa_g*kg_hat(-1) - e_a ;  

// 23. Inflation (domestic) and relative prices
pic_hat = ((omega_d)*(1/ phimc)^(1-etac)) * pid_hat + ((1-omega_d)*(1/ phimi)^(1-etai)) * piim_hat; //

// 24-25. Imports and Exports
mcx_hat = mcx_hat(-1) + pid_hat - pix_hat - deltaE_hat; //
mcmi_hat = mcmi_hat(-1) + pistar_hat - pid_hat + deltaE_hat;  //-mcx_hat 
x_hat =  -omega_d*(( (1-omega_d)*(1/phimc)^(1-etac) + omega_d )^(1/(1-etac)))*pid_hat - pistar_hat;
im_hat = (c_bar/(c_bar+I_bar))*(c_hat-c_hat(-1))-(c_bar/(c_bar+I_bar))*(etac*(1-omega_d)) + (I_bar/(c_bar+I_bar))*(I_hat-I_hat(-1))-(I_bar/(c_bar+I_bar))*(etai*(1-omega_im));

// Oil block
//phi_oil, omega_oil, oil_bar, oil_hat;
//oil_hat = rho_oil * oil_hat(-1) + e_oil;
//tot_hat = y_hat + oil_hat; 

// 24-27. Government budget constraint and taxes
g_hat*g_bar+Ig_bar*Ig_hat+tr_bar*trans+b_bar*R_bar*(b_hat(-1)+R_hat(-1)-pid_hat) =b_bar*b_hat+tao_cbar*c_bar*(c_hat+tao_c)+tao_wbar*wN_bar*(w_hat+N_hat+tao_w)+tao_kbar*rky_bar*(rk_hat+u_hat+k_hat(-1)+tao_k) ; // omega_oil*yoil_hat*oil_bar
tao_c_inc = (c_hat+tao_c);                          // Consumption tax  
tao_k_inc = (rk_hat+u_hat+k_hat(-1)+tao_k);         // Capital tax
tao_w_inc = (w_hat+N_hat+tao_w);                    // Labor income tax

// Debt-to-output Ratio
b_hat_y = b_bar^(-1)*b_hat - y_hat; 

// 30-33. Fiscal policy rules:
// Change to (+/-) model base or shocks
g_hat  = rho_g*(g_hat(-1)) - (1-rho_g)*(psi_g*y_hat+ela_g*(b_hat(-1))) - e_g;             
Ig_hat = rho_ig*(Ig_hat(-1)) - (1-rho_ig)*(psi_ig*y_hat+ela_ig*(b_hat(-1))) - e_ig;    // 
trans  = rho_tr*(trans(-1)) - (1-rho_tr)*(psi_tr*y_hat+ela_tr*(b_hat(-1))) - e_tr;     // 

tao_c = rho_tc*(tao_c(-1)) + (1-rho_tc)*(psi_tc*(y_hat)+ela_lc*(b_hat(-1))) + e_tc;    // 
tao_k = rho_tk*(tao_k(-1)) + (1-rho_tk)*(psi_tk*(y_hat)+ela_lk*(b_hat(-1))) + e_tk;    // 
tao_w = rho_tw*(tao_w(-1)) + (1-rho_tw)*(psi_tw*(y_hat)+ela_lw*(b_hat(-1))) + e_tw;    //

// Net trade balance
nx_hat = y_hat - omega_x*deltaE_hat - c_hat - g_hat; % 

// -----------------------------------------------------------------------

//38-40 Foreign variables - VAR model
pistar_hat  = rho_pistar * pistar_hat(-1) + e_pistar;
ystar_hat   = rho_ystar  * ystar_hat(-1)  + e_ystar ;
Rstar_hat   = rho_Rstar  * Rstar_hat(-1)  + e_Rstar ;

//terms of trade:
// p_B = (infl_B/infl_A)*p_B(-1);

// 41-55. EXOGENOUS VARIABLES
e_tc = phi_tc*e_tc(-1)+tct;
e_tk = phi_tk*e_tk(-1)+tkt;
e_tw = phi_tw*e_tw(-1)+twt;
e_tr = phi_tr*e_tr(-1)+trt;
e_ig = phi_ig*e_ig(-1)+igt;
e_g  = phi_g*e_g(-1)+gt;
e_a  = rho_a*e_a(-1)+at;
e_i  = rho_i*e_i(-1)+it;
e_w  = rho_w*e_w(-1)+wt;
e_pi = rho_p*e_pi(-1)+pit;
e_pix= rho_px*e_pix(-1)+pixt;
e_pim= rho_pm*e_pim(-1)+pimt;
e_n  = rho_n*e_n(-1)+nt;
e_r  = rt;
e_dt = dt;
phitilde_hat = rho_phitilde*phitilde_hat(-1)+phitildet; //Risk premium
//e_oil = oilt;

// These denote various stochastic shocks to technology, investment, interest rate, inflation, export prices, 
// import prices, labor, wages, transfers, government investment, consumption tax, capital tax, labor tax, and
// government spending.

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// 58-69. Observables:
data_y     = y_hat;  % y_hat-y_hat(-1) + 100*log(yz_bar);
data_pi    = pic_hat;
data_con   = c_hat+(tao_cbar/(1+tao_cbar))*(tao_c);
data_inv   = I_hat;
data_gov   = g_hat;
data_r     = R_hat; 
data_rstar = Rstar_hat;
data_ystar = ystar_hat;
data_pistar= pistar_hat;
end;

// ---------------------------------------------------------------------
// -----------------------------------------------

initval;
%g_hat   = g_bar;      % 0.16
%Ig_hat  = Ig_bar;	   % 0.02
%trans   = tr_bar;	  
b_hat   =  b_bar;      %2.2
%c_hat   = c_bar; 
%I_hat   = I_bar;
%R_hat   = R_bar;       
%pic_hat = piT;         
%tao_c = tao_cbar;        
%tao_k = tao_kbar ;      
%tao_w = tao_wbar ;      
end;

// -----------------------------------------------


model_diagnostics;
%steady;
%check;

resid(1);
%Options_.dynatol.f=3e-4;
options.solve_tolf= 5e-3;
steady(solve_algo = 2, maxit = 1000000);
check;

// Shocks
shocks;
// Change values to policy shocks to 1pp
var tct; stderr std_tc;
var tkt; stderr std_tk;
var twt; stderr std_tw;
var gt; stderr std_g;
var igt; stderr std_ig;
var trt; stderr std_tr;
%var at; stderr std_at;
%var wt; stderr std_wt;
%var nt; stderr std_nt;
var pit; stderr std_pit;
%var it; stderr std_it;
var rt; stderr std_rt;
%var pixt; stderr std_pixt;
%var pimt; stderr std_pimt;
%var e_pistar; stderr std_pixt;
%var e_ystar; stderr std_pimt;
%var e_Rstar; stderr std_rt;
%var phitildet; stderr std_phitilde;
%var dt; stderr std_dt;
end;

estimated_params ;
    stderr at      ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr nt      ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr wt      ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr pit     ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr tct     ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr tkt     ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr twt     ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr it      ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr trt     ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr igt     ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr gt      ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr rt      ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr pixt    ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr pimt    ,        inv_gamma_pdf,    0.0125     ,        Inf;
    stderr phitildet,       inv_gamma_pdf,    0.0125     ,        Inf;
    //stderr oilt   ,        inv_gamma_pdf,    0.0125     ,        Inf;
    rho_tc         ,        beta_pdf,         0.80       ,        0.1;
    rho_tk         ,        beta_pdf,         0.80       ,        0.1;
    rho_tw         ,        beta_pdf,         0.80       ,        0.1;
    rho_g          ,        beta_pdf,         0.80       ,        0.1;
    rho_ig         ,        beta_pdf,         0.80       ,        0.1;
    rho_tr         ,        beta_pdf,         0.80       ,        0.1;
    rho_a          ,        beta_pdf,         0.80       ,        0.1;
    rho_n          ,        beta_pdf,         0.80       ,        0.1;
    rho_i          ,        beta_pdf,         0.80       ,        0.1;
    rho_pistar     ,        beta_pdf,         0.80       ,        0.1;
    rho_ystar      ,        beta_pdf,         0.80       ,        0.1;
    rho_Rstar      ,        beta_pdf,         0.80       ,        0.1;
    rho_phitilde   ,        beta_pdf,         0.80       ,        0.1;
    phi_g          ,        beta_pdf,         0.80       ,        0.1;
    phi_ig         ,        beta_pdf,         0.80       ,        0.1;
    phi_tr         ,        beta_pdf,         0.80       ,        0.1;
    phi_tk         ,        beta_pdf,         0.80       ,        0.1;
    phi_tw         ,        beta_pdf,         0.80       ,        0.1;
    phi_tc         ,        beta_pdf,         0.80       ,        0.1;
    //rho_oil        ,        beta_pdf,         0.70       ,        0.1;
    rho_w          ,        beta_pdf,         0.70       ,        0.1;
    rho_p          ,        beta_pdf,         0.70       ,        0.1;
    rho_px         ,        beta_pdf,         0.70       ,        0.1;
    rho_pm         ,        beta_pdf,         0.70       ,        0.1;
    phitildes      ,        beta_pdf,         0.70       ,        0.1;
    omega_d        ,        beta_pdf,         0.50       ,        0.05;
    omega_w        ,        beta_pdf,         0.50       ,        0.05;
    omega_im       ,        beta_pdf,         0.50       ,        0.05;
    omega_x        ,        beta_pdf,         0.50       ,        0.05;
    chi_h          ,        beta_pdf,         0.50       ,        0.05;
    chi_w          ,        beta_pdf,         0.50       ,        0.05;
    chi_x          ,        beta_pdf,         0.50       ,        0.05;
    chi_e          ,        beta_pdf,         0.50       ,        0.05;
    chi_im         ,        beta_pdf,         0.50       ,        0.05;
    r_pi           ,        normal_pdf,       1.5        ,        0.1;
    r_y            ,        normal_pdf,       0.125      ,        0.05;
    fc             ,        normal_pdf,       1.2        ,        0.1;
    phimc          ,        normal_pdf,       1.2        ,        0.1;
    phimi          ,        normal_pdf,       1.2        ,        0.1;
    etac           ,        normal_pdf,       1.2        ,        0.1;
    etai           ,        normal_pdf,       1.2        ,        0.1;
    kappa          ,        normal_pdf,       0.7        ,        0.1;
    psi_g          ,        normal_pdf,       0.5        ,        0.1;
    psi_ig         ,        normal_pdf,       0.5        ,        0.1;
    psi_tr         ,        normal_pdf,       0.5        ,        0.1;
    psi_tc         ,        normal_pdf,       0.5        ,        0.1;
    psi_tk         ,        normal_pdf,       0.5        ,        0.1;
    psi_tw         ,        normal_pdf,       0.5        ,        0.1;
    ela_g          ,        normal_pdf,       0.2        ,        0.1;
    ela_ig         ,        normal_pdf,       0.2        ,        0.1;
    ela_tr         ,        normal_pdf,       0.2        ,        0.1;
    ela_lw         ,        normal_pdf,       0.2        ,        0.1;  
    ela_lc         ,        normal_pdf,       0.2        ,        0.1;  
    ela_lk         ,        normal_pdf,       0.2        ,        0.1; 
    etax           ,        normal_pdf,      1.2        ,        0.1; 
end;


// ---------------------------------------------------------------------

varobs 
data_y data_pi data_r data_con data_gov data_inv data_pistar data_rstar data_ystar; 
// data_tr data_hours data_taxc data_taxk data_taxw data_wages data_inv_g data_pistar data_rstar data_ystar
// data_pet data_ex data_imp
// data_rev data_debt 

options_.order=1;
%options_.aim_solver=1;
options_.prior_trunc=0;


// Parameter Estimations
%estimation(datafile=datacol,plot_priors=0,mh_nblocks=4,mh_replic=100000,mh_jscale=0.29,lik_init=2,mode_compute=4); % mh_drop=.5
//estimation(datafile=datacol,plot_priors=0,mh_nblocks=2,mh_replic=75000,mh_jscale=0.29,lik_init=2,mode_compute=4) ;

// Schock simul..
stoch_simul (order=1,nofunctions,irf=40, nocorr, noprint, nograph);
//stoch_simul(order=1,nograph,irf=40) y_hat R_hat pid_hat c_hat cr_hat cnr_hat I_hat k_hat w_hat emp_hat;  //hp_filter=1600

// Simul
%simul(periods=40); %stack_solve_algo=0

%%Shock_decomposition;
//shock_decomposition y_hat pid_hat c_hat I_hat g_hat Ig_hat trans tao_c tao_k tao_w;

//estimation(datafile=datases,mode_compute=4,
//nobs=100,first_obs=1,mh_init_scale = 0.8, mh_replic=5000,
//mh_nblocks=2,mh_jscale=0.4,mode_check,
//conditional_variance_decomposition=[1,4,8,16,32]);
//write_latex_prior_table;  

%write_latex_parameter_table;
%write_latex_dynamic_model;
%write_latex_definitions;
%collect_latex_files;


//--------------------------------------------------------------------------------------------------------------

// Fiscal anaylisis
// g_hat, Ig_hat, trans || tao_c, tao_k, tao_w 

//define planner objective
%options_.ramsey.maxit = 10000;
%planner_objective beta*((cr_hat + (sigma_c/2)*(cr_hat)^2) - ((1-tao_wbar)/(1+tao_cbar)) * (N_hat+(sigma_n/2)*(N_hat)^2));
//set up Ramsey optimal policy problem with fiscal rules
%ramsey_model(planner_discount=beta,instruments=(tao_c, tao_w, tao_k, g_hat, Ig_hat, trans)); %, tao_c, tao_w, tao_k, g_hat, Ig_hat, trans
//conduct stochastic simulations of the Ramsey problem
%stoch_simul (order=2,nofunctions,irf=40, periods=500, nocorr, noprint, nograph);
//simulation
%steady(solve_algo = 2, maxit = 1000000);
// Simul. with perfect foresight
%perfect_foresight_setup(periods=500);
%perfect_foresight_solver;
%evaluate_planner_objective;
%----------------------------------------------------------------------------------------------


// End ----
// ----------------------------------------------------------------------------------------------
