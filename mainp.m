%%
%% -------------------------------------------------------------------
%% -------------------------------------------------------------------
% This .m file reproduces all output of the DSGE model and fiscal
% Multipliers 
% Copyright (C) 2024
% 
%% -------------------------------------------------------------------
%%-------------------------------------------------------------------
%
% Run for sections 
%
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

% This .m file uses Matlab Software  and Dyrane 5.4
clc
clear all;
warning off;
% Import data
datacol;
% Run Dynanre .mod code
addpath C:/dynare/5.4/matlab/     %% Adjust this line to your current Dynare version 

% Base line model
dynare code                       %% Run this line to the base line of the model

% Adjusted models: sections 7 and 8.
%dynare code_policy
%dynare code_oil

%% -------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%% Tax increase effect

y_hat_tct = oo_.irfs.y_hat_tct;
y_hat_twt = oo_.irfs.y_hat_twt;
y_hat_tkt = oo_.irfs.y_hat_tkt;
pid_hat_tct = oo_.irfs.pid_hat_tct;
pid_hat_twt = oo_.irfs.pid_hat_twt;
pid_hat_tkt = oo_.irfs.pid_hat_tkt;
R_hat_tct = oo_.irfs.R_hat_tct;
R_hat_twt = oo_.irfs.R_hat_twt;
R_hat_tkt = oo_.irfs.R_hat_tkt;
cr_hat_tct = oo_.irfs.cr_hat_tct;
cr_hat_twt = oo_.irfs.cr_hat_twt;
cr_hat_tkt = oo_.irfs.cr_hat_tkt;
cnr_hat_tct = oo_.irfs.cnr_hat_tct;
cnr_hat_twt = oo_.irfs.cnr_hat_twt;
cnr_hat_tkt = oo_.irfs.cnr_hat_tkt;
w_hat_tct = oo_.irfs.w_hat_tct;
w_hat_twt = oo_.irfs.w_hat_twt;
w_hat_tkt = oo_.irfs.w_hat_tkt;
emp_hat_tct = oo_.irfs.emp_hat_tct;
emp_hat_twt = oo_.irfs.emp_hat_twt;
emp_hat_tkt = oo_.irfs.emp_hat_tkt;
I_hat_tct = oo_.irfs.I_hat_tct;
I_hat_twt = oo_.irfs.I_hat_twt;
I_hat_tkt = oo_.irfs.I_hat_tkt;
k_hat_tct = oo_.irfs.k_hat_tct;
k_hat_twt = oo_.irfs.k_hat_twt;
k_hat_tkt = oo_.irfs.k_hat_tkt;
len = length(y_hat_tct);

% Figure 1
figure
subplot(3, 3, 1)
hold on
plot(y_hat_tct, "LineWidth", 1, "Color", 'k')
plot(y_hat_twt, "LineWidth", 1, "Color", 'b')
plot(y_hat_tkt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Output')
hold off
subplot(3,3, 2)
hold on
plot(pid_hat_tct, "LineWidth", 1, "Color", 'k')
plot(pid_hat_twt, "LineWidth", 1, "Color", 'b')
plot(pid_hat_tkt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Inflation')
hold off
subplot(3, 3, 3)
hold on
plot(cr_hat_tct, "LineWidth", 1, "Color", 'k')
plot(cr_hat_twt, "LineWidth", 1, "Color", 'b')
plot(cr_hat_tkt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Ric. Consumpt')
legend({'$\tau_c$','$\tau_w$','$\tau_k$'},'Location','best','Interpreter','latex'); %'Location','southeast'
legend('boxoff')
hold off
subplot(3, 3, 4)
hold on
plot(cnr_hat_tct, "LineWidth", 1, "Color", 'k')
plot(cnr_hat_twt, "LineWidth", 1, "Color", 'b')
plot(cnr_hat_tkt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('non-Ric. Consumpt.')
hold off
subplot(3, 3, 5)
hold on
plot(w_hat_tct, "LineWidth", 1, "Color", 'k')
plot(w_hat_twt, "LineWidth", 1, "Color", 'b')
plot(w_hat_tkt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Wage')
hold off
subplot(3,3,6)
hold on
plot(emp_hat_tct, "LineWidth", 1, "Color", 'k')
plot(emp_hat_twt, "LineWidth", 1, "Color", 'b')
plot(emp_hat_tkt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Employment')
hold off
subplot(3,3,7)
hold on
plot(I_hat_tct, "LineWidth", 1, "Color", 'k')
plot(I_hat_twt, "LineWidth", 1, "Color", 'b')
plot(I_hat_tkt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Investment')
hold off
subplot(3,3,8)
hold on
plot(k_hat_tct, "LineWidth", 1, "Color", 'k')
plot(k_hat_twt, "LineWidth", 1, "Color", 'b')
plot(k_hat_tkt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Capital')
hold off
subplot(3,3,9)
hold on
plot(b_hat_tct, "LineWidth", 1, "Color", 'k')
plot(b_hat_twt, "LineWidth", 1, "Color", 'b')
plot(b_hat_tkt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Debt')
hold off


%% -------------------------------------------------------------------
%% Government expendiature effect

y_hat_gt = oo_.irfs.y_hat_gt;
y_hat_igt = oo_.irfs.y_hat_igt;
y_hat_trt = oo_.irfs.y_hat_trt;
pid_hat_gt = oo_.irfs.pid_hat_gt;
pid_hat_igt = oo_.irfs.pid_hat_igt;
pid_hat_trt = oo_.irfs.pid_hat_trt;
R_hat_gt = oo_.irfs.R_hat_gt;
R_hat_igt = oo_.irfs.R_hat_igt;
R_hat_trt = oo_.irfs.R_hat_trt;
cr_hat_gt = oo_.irfs.cr_hat_gt;
cr_hat_igt = oo_.irfs.cr_hat_igt;
cr_hat_trt = oo_.irfs.cr_hat_trt;
cnr_hat_gt = oo_.irfs.cnr_hat_gt;
cnr_hat_igt = oo_.irfs.cnr_hat_igt;
cnr_hat_trt = oo_.irfs.cnr_hat_trt;
w_hat_gt = oo_.irfs.w_hat_gt;
w_hat_igt = oo_.irfs.w_hat_igt;
w_hat_trt = oo_.irfs.w_hat_trt;
emp_hat_gt = oo_.irfs.emp_hat_gt;
emp_hat_igt = oo_.irfs.emp_hat_igt;
emp_hat_trt = oo_.irfs.emp_hat_trt;
I_hat_gt = oo_.irfs.I_hat_gt;
I_hat_igt = oo_.irfs.I_hat_igt;
I_hat_trt = oo_.irfs.I_hat_trt;
k_hat_gt = oo_.irfs.k_hat_gt;
k_hat_igt = oo_.irfs.k_hat_igt;
k_hat_trt = oo_.irfs.k_hat_trt;

% Figure 1
figure
subplot(3, 3, 1)
hold on
plot(y_hat_gt, "LineWidth", 1, "Color", 'k')
plot(y_hat_igt, "LineWidth", 1, "Color", 'b')
plot(y_hat_trt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Output')
hold off
subplot(3,3, 2)
hold on
plot(pid_hat_gt, "LineWidth", 1, "Color", 'k')
plot(pid_hat_igt, "LineWidth", 1, "Color", 'b')
plot(pid_hat_trt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Inflation')
hold off
subplot(3, 3, 3)
hold on
plot(cr_hat_gt, "LineWidth", 1, "Color", 'k')
plot(cr_hat_igt, "LineWidth", 1, "Color", 'b')
plot(cr_hat_trt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Ric. Consumpt')
legend({'$g$','$Ig$','$tr$'},'Location','best','Interpreter','latex'); %'Location','southeast'
legend('boxoff')
hold off
subplot(3, 3, 4)
hold on
plot(cnr_hat_gt, "LineWidth", 1, "Color", 'k')
plot(cnr_hat_igt, "LineWidth", 1, "Color", 'b')
plot(cnr_hat_trt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('non-Ric. Consumpt.')
hold off
subplot(3, 3, 5)
hold on
plot(w_hat_gt, "LineWidth", 1, "Color", 'k')
plot(w_hat_igt, "LineWidth", 1, "Color", 'b')
plot(w_hat_trt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Wage')
hold off
subplot(3,3,6)
hold on
plot(emp_hat_gt, "LineWidth", 1, "Color", 'k')
plot(emp_hat_igt, "LineWidth", 1, "Color", 'b')
plot(emp_hat_trt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Employment')
hold off
subplot(3,3,7)
hold on
plot(I_hat_gt, "LineWidth", 1, "Color", 'k')
plot(I_hat_igt, "LineWidth", 1, "Color", 'b')
plot(I_hat_trt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Investment')
hold off
subplot(3,3,8)
hold on
plot(k_hat_gt, "LineWidth", 1, "Color", 'k')
plot(k_hat_igt, "LineWidth", 1, "Color", 'b')
plot(k_hat_trt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Capital')
hold off
subplot(3,3,9)
hold on
plot(b_hat_gt, "LineWidth", 1, "Color", 'k')
plot(b_hat_igt, "LineWidth", 1, "Color", 'b')
plot(b_hat_trt, "LineWidth", 1, "Color", 'r')
yline(0,'--')
xlabel('Quarters')
title('Debt')
hold off


%%
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% Monetary policy shock
 
figure
subplot(3, 3, 1)
hold on
plot(y_hat_rt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Output')
hold off
subplot(3, 3, 2)
hold on
plot(pid_hat_rt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Inflation')
hold off
subplot(3, 3, 3)
hold on
plot(R_hat_rt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Interest rate')
hold off
subplot(3, 3, 4)
hold on
plot(cr_hat_rt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Ric. Consumpt')
hold off
subplot(3, 3, 5)
hold on
plot(cnr_hat_rt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('non-Ric. Consumpt')
hold off
subplot(3, 3, 6)
hold on
plot(w_hat_rt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Wage')
hold off
subplot(3, 3, 7)
hold on
plot(emp_hat_rt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Employment')
hold off
subplot(3, 3, 8)
hold on
plot(I_hat_rt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Investment')
hold off
subplot(3, 3, 9)
hold on
plot(k_hat_rt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Capital')
hold off

%% 
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------


%% 
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% Observed variables 

figure
hold on
subplot(3, 3, 1)
plot(data_y, "LineWidth", 1, "Color", 'k')
title('Output')
xlabel('Quarters')
subplot(3, 3, 2)
plot(data_pi, "LineWidth", 1, "Color", 'k')
title('Inflation')
xlabel('Quarters')
subplot(3, 3, 3)
plot(data_r, "LineWidth", 1, "Color", 'k')
title('Interest Rate')
xlabel('Quarters')
subplot(3, 3, 4)
plot(data_con, "LineWidth", 1, "Color", 'k')
title('Consumption')
xlabel('Quarters')
subplot(3, 3, 5)
plot(data_inv, "LineWidth", 1, "Color", 'k')
title('Investment')
xlabel('Quarters')
subplot(3, 3, 6)
plot(data_gov, "LineWidth", 1, "Color", 'k')
title('Gov. compsumption')
xlabel('Quarters')
subplot(3, 3, 7)
plot(data_ystar, "LineWidth", 1, "Color", 'k')
title('Ext. Output')
xlabel('Quarters')
subplot(3, 3, 8)
plot(data_pistar, "LineWidth", 1, "Color", 'k')
title('Ext. Inflation')
xlabel('Quarters')
subplot(3, 3, 9)
plot(data_rstar, "LineWidth", 1, "Color", 'k')
title('Ext. Interest Rate')
xlabel('Quarters')
axis tight
hold off;


%%
%% IRF results and Fiscal multipliers
%% tao_c: 28 | tao_k: 29 | tao_w: 30 // g_hat: 25 | Ig_hat: 26 | trans: 27 //  y_hat: 10 | data_y: 49

%rate_gov = 1+R_bar/400;
%disc = 1/rate_gov;
clc
len = length(y_hat_gt);
%beta = 0.9922 ;
disc_vec = beta.^[0:1:len-1];  % beta | disc
%disc_mat=repmat(disc_vec,size(resp_y,1),1);

% Multipliers Gov. Expend.
resp_yg = disc_vec.*(oo_.irfs.y_hat_gt);
resp_rg  = disc_vec.*(oo_.irfs.g_hat_gt);
resp_ig = disc_vec.*(oo_.irfs.I_hat_gt);
resp_cg = disc_vec.*(oo_.irfs.c_hat_gt);
resp_rig  = disc_vec.*(oo_.irfs.Ig_hat_igt);
resp_yig = disc_vec.*(oo_.irfs.y_hat_igt);
resp_iig = disc_vec.*(oo_.irfs.I_hat_igt);
resp_cig = disc_vec.*(oo_.irfs.c_hat_igt);
resp_rtr  = disc_vec.*(oo_.irfs.trans_trt);
resp_ytr = disc_vec.*(oo_.irfs.y_hat_trt);
resp_ctr = disc_vec.*(oo_.irfs.c_hat_trt);
resp_itr = disc_vec.*(oo_.irfs.Ig_hat_trt);

% Gov. Expend.
% GDP
mult_yg  = [sum(resp_yg(1:4))./sum(resp_rg(1:4)), sum(resp_yg(1:8))./sum(resp_rg(1:8)), sum(resp_yg(1:12))./sum(resp_rg(1:12)), sum(resp_yg(1:20))./sum(resp_rg(1:20)), sum(resp_yg(1:end))./sum(resp_rg(1:end))];
mult_yig = [sum(resp_yig(1:4))./sum(resp_rig(1:4)), sum(resp_yig(1:8))./sum(resp_rig(1:8)), sum(resp_yig(1:12))./sum(resp_rig(1:12)), sum(resp_yig(1:20))./sum(resp_rig(1:20)), sum(resp_yig(1:end))./sum(resp_rig(1:end))];
mult_ytr = [sum(resp_ytr(1:4))./sum(resp_rtr(1:4)), sum(resp_ytr(1:8))./sum(resp_rtr(1:8)), sum(resp_ytr(1:12))./sum(resp_rtr(1:12)), sum(resp_ytr(1:20))./sum(resp_rtr(1:20)), sum(resp_ytr(1:end))./sum(resp_rtr(1:end))];

% Compsumption
mult_cg  = [sum(resp_cg(1:4))./sum(resp_rg(1:4)), sum(resp_cg(1:8))./sum(resp_rg(1:8)), sum(resp_cg(1:12))./sum(resp_rg(1:12)), sum(resp_cg(1:20))./sum(resp_rg(1:20)), sum(resp_cg(1:end))./sum(resp_rg(1:end))];
mult_cig = [sum(resp_cig(1:4))./sum(resp_rig(1:4)), sum(resp_cig(1:8))./sum(resp_rig(1:8)), sum(resp_cig(1:12))./sum(resp_rig(1:12)), sum(resp_cig(1:20))./sum(resp_rig(1:20)), sum(resp_cig(1:end))./sum(resp_rig(1:end))];
mult_ctr = [sum(resp_ctr(1:4))./sum(resp_rtr(1:4)), sum(resp_ctr(1:8))./sum(resp_rtr(1:8)), sum(resp_ctr(1:12))./sum(resp_rtr(1:12)), sum(resp_ctr(1:20))./sum(resp_rtr(1:20)), sum(resp_ctr(1:end))./sum(resp_rtr(1:end))];

% Investment
mult_invg  = [sum(resp_ig(1:4))./sum(resp_rg(1:4)), sum(resp_ig(1:8))./sum(resp_rg(1:8)), sum(resp_ig(1:12))./sum(resp_rg(1:12)), sum(resp_ig(1:20))./sum(resp_rg(1:20)), sum(resp_ig(1:end))./sum(resp_rg(1:end))];
mult_invig = [sum(resp_iig(1:4))./sum(resp_rig(1:4)), sum(resp_iig(1:8))./sum(resp_rig(1:8)), sum(resp_iig(1:12))./sum(resp_rig(1:12)), sum(resp_iig(1:20))./sum(resp_rig(1:20)), sum(resp_iig(1:end))./sum(resp_rig(1:end))];
mult_invtr = [sum(resp_itr(1:4))./sum(resp_rtr(1:4)), sum(resp_itr(1:8))./sum(resp_rtr(1:8)), sum(resp_itr(1:12))./sum(resp_rtr(1:12)), sum(resp_itr(1:20))./sum(resp_rtr(1:20)), sum(resp_itr(1:end))./sum(resp_rtr(1:end))];

% Fiscal multip.
% GDP
mult_fis_pib = [mult_yg; mult_yig; mult_ytr];
disp(' ');
disp('PV Fiscal Multipliers - GDP (1, 2, 3, 5 and 10 yrs):');
disp('------------------------------------------');
disp(mult_fis_pib);
disp(' ');

% Compsumpt.
mult_fis_cons = [mult_cg; mult_cig; mult_ctr];
disp(' ');
disp('PV Fiscal Multipliers - Consum. (1, 2, 3, 5 and 10 yrs):');
disp('------------------------------------------');
disp(mult_fis_cons);
disp(' ');

% Investment
mult_fis_inv = [mult_invg; mult_invig; mult_invtr];
disp(' ');
disp('PV Fiscal Multipliers - Invest. (1, 2, 3, 5 and 10 yrs):');
disp('------------------------------------------');
disp(mult_fis_inv);
disp(' ');

% Mean multip.
disp('Mean MF GDP.: ');
mean(mult_fis_pib,2)

disp('Mean MF Compsum.: ');
mean(mult_fis_cons,2)

disp('Mean MF Invest.: ');
mean(mult_fis_inv,2)

%% 
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% Tax multipliers

clc;
% Taxes
resp_rtaxc = disc_vec.*(oo_.irfs.tao_c_tct);
resp_ytaxc = disc_vec.*(oo_.irfs.y_hat_tct);
resp_ctaxc = disc_vec.*(oo_.irfs.c_hat_tct);
resp_itaxc = disc_vec.*(oo_.irfs.I_hat_tct);
resp_rtaxk = disc_vec.*(oo_.irfs.tao_k_tkt);
resp_ytaxk = disc_vec.*(oo_.irfs.y_hat_tkt);
resp_ctaxk = disc_vec.*(oo_.irfs.c_hat_tkt);
resp_itaxk = disc_vec.*(oo_.irfs.I_hat_tkt);
resp_rtaxw = disc_vec.*(oo_.irfs.tao_w_twt);
resp_ytaxw = disc_vec.*(oo_.irfs.y_hat_twt);
resp_ctaxw = disc_vec.*(oo_.irfs.c_hat_twt);
resp_itaxw = disc_vec.*(oo_.irfs.I_hat_twt);

% GDP
mult_ytaxc  = [sum(resp_ytaxc(1:4))./sum(resp_rtaxc(1:4)), sum(resp_ytaxc(1:8))./sum(resp_rtaxc(1:8)), sum(resp_ytaxc(1:12))./sum(resp_rtaxc(1:12)), sum(resp_ytaxc(1:20))./sum(resp_rtaxc(1:20)), sum(resp_ytaxc(1:end))./sum(resp_rtaxc(1:end))];
mult_ytaxk  = [sum(resp_ytaxk(1:4))./sum(resp_rtaxk(1:4)), sum(resp_ytaxk(1:8))./sum(resp_rtaxk(1:8)), sum(resp_ytaxk(1:12))./sum(resp_rtaxk(1:12)), sum(resp_ytaxk(1:20))./sum(resp_rtaxk(1:20)), sum(resp_ytaxk(1:end))./sum(resp_rtaxk(1:end))];
mult_ytaxw  = [sum(resp_ytaxw(1:4))./sum(resp_rtaxw(1:4)), sum(resp_ytaxw(1:8))./sum(resp_rtaxw(1:8)), sum(resp_ytaxw(1:12))./sum(resp_rtaxw(1:12)), sum(resp_ytaxw(1:20))./sum(resp_rtaxw(1:20)), sum(resp_ytaxw(1:end))./sum(resp_rtaxw(1:end))];

% Compsumpt.
mult_ctaxc  = [sum(resp_ctaxc(1:4))./sum(resp_rtaxc(1:4)), sum(resp_ctaxc(1:8))./sum(resp_rtaxc(1:8)), sum(resp_ctaxc(1:12))./sum(resp_rtaxc(1:12)), sum(resp_ctaxc(1:20))./sum(resp_rtaxc(1:20)), sum(resp_ctaxc(1:end))./sum(resp_rtaxc(1:end))];
mult_ctaxk  = [sum(resp_ctaxk(1:4))./sum(resp_rtaxk(1:4)), sum(resp_ctaxk(1:8))./sum(resp_rtaxk(1:8)), sum(resp_ctaxk(1:12))./sum(resp_rtaxk(1:12)), sum(resp_ctaxk(1:20))./sum(resp_rtaxk(1:20)), sum(resp_ctaxk(1:end))./sum(resp_rtaxk(1:end))];
mult_ctaxw  = [sum(resp_ctaxw(1:4))./sum(resp_rtaxw(1:4)), sum(resp_ctaxw(1:8))./sum(resp_rtaxw(1:8)), sum(resp_ctaxw(1:12))./sum(resp_rtaxw(1:12)), sum(resp_ctaxw(1:20))./sum(resp_rtaxw(1:20)), sum(resp_ctaxw(1:end))./sum(resp_rtaxw(1:end))];

% Investm.
mult_itaxc  = [sum(resp_itaxc(1:4))./sum(resp_rtaxc(1:4)), sum(resp_itaxc(1:8))./sum(resp_rtaxc(1:8)), sum(resp_itaxc(1:12))./sum(resp_rtaxc(1:12)), sum(resp_itaxc(1:20))./sum(resp_rtaxc(1:20)), sum(resp_itaxc(1:end))./sum(resp_rtaxc(1:end))];
mult_itaxk  = [sum(resp_itaxk(1:4))./sum(resp_rtaxk(1:4)), sum(resp_itaxk(1:8))./sum(resp_rtaxk(1:8)), sum(resp_itaxk(1:12))./sum(resp_rtaxk(1:12)), sum(resp_itaxk(1:20))./sum(resp_rtaxk(1:20)), sum(resp_itaxk(1:end))./sum(resp_rtaxk(1:end))];
mult_itaxw  = [sum(resp_itaxw(1:4))./sum(resp_rtaxw(1:4)), sum(resp_itaxw(1:8))./sum(resp_rtaxw(1:8)), sum(resp_itaxw(1:12))./sum(resp_rtaxw(1:12)), sum(resp_itaxw(1:20))./sum(resp_rtaxw(1:20)), sum(resp_itaxw(1:end))./sum(resp_rtaxw(1:end))];

% Multip.
% GDP
mult_fis_pib_tax = [mult_ytaxc; mult_ytaxk; mult_ytaxw];
disp(' ');
disp('PV Fiscal Multipliers - GDP (1, 2, 3, 5 and 10 yrs):');
disp('------------------------------------------');
disp(mult_fis_pib_tax);
disp(' ');

% Compsumpt.
mult_fis_cons_tax = [mult_ctaxc; mult_ctaxk; mult_ctaxw];
disp(' ');
disp('PV Fiscal Multipliers - Consum. (1, 2, 3, 5 and 10 yrs):');
disp('------------------------------------------');
disp(mult_fis_cons_tax);
disp(' ');

% Investm.
mult_fis_inv_tax = [mult_itaxc; mult_itaxk; mult_itaxw];
disp(' ');
disp('PV Fiscal Multipliers - Invest. (1, 2, 3, 5 and 10 yrs):');
disp('------------------------------------------');
disp(mult_fis_inv_tax);
disp(' ');

% Mean multip.
disp('Mean MF GDP.: ');
mean(mult_fis_pib_tax,2)

disp('Mean MF Compsum.: ');
mean(mult_fis_cons_tax,2)

disp('Mean MF Invest.: ');
mean(mult_fis_inv_tax,2)


%%
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% Fiscal consolidation
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

b_current = 2.24;
% 1 pp increase in taxes
% tao_c
y_hat_tct = oo_.irfs.y_hat_tct;
c_hat_tct = oo_.irfs.c_hat_tct;
cr_hat_tct = oo_.irfs.cr_hat_tct;
cnr_hat_tct = oo_.irfs.cnr_hat_tct;
w_hat_tct = oo_.irfs.w_hat_tct;
emp_hat_tct = oo_.irfs.emp_hat_tct;
I_hat_tct = oo_.irfs.I_hat_tct;
k_hat_tct = oo_.irfs.k_hat_tct;
b_hat_y_tct = oo_.irfs.b_hat_y_tct;

deltatauc = [y_hat_tct; c_hat_tct; cr_hat_tct; cnr_hat_tct; w_hat_tct; emp_hat_tct; I_hat_tct; k_hat_tct; b_hat_y_tct];
deltatauc = transpose(deltatauc);

% tao_w
y_hat_twt = oo_.irfs.y_hat_twt;
c_hat_twt = oo_.irfs.c_hat_twt;
cr_hat_twt = oo_.irfs.cr_hat_twt;
cnr_hat_twt = oo_.irfs.cnr_hat_twt;
w_hat_twt = oo_.irfs.w_hat_twt;
emp_hat_twt = oo_.irfs.emp_hat_twt;
I_hat_twt = oo_.irfs.I_hat_twt;
k_hat_twt = oo_.irfs.k_hat_twt;
b_hat_y_twt = oo_.irfs.b_hat_y_twt;

deltatauw = [y_hat_twt; c_hat_twt; cr_hat_twt; cnr_hat_twt; w_hat_twt; emp_hat_twt; I_hat_twt; k_hat_twt; b_hat_y_twt];
deltatauw = transpose(deltatauw);

% tao_k
y_hat_tkt = oo_.irfs.y_hat_tkt;
c_hat_tkt = oo_.irfs.c_hat_tkt;
cr_hat_tkt = oo_.irfs.cr_hat_tkt;
cnr_hat_tkt = oo_.irfs.cnr_hat_tkt;
w_hat_tkt = oo_.irfs.w_hat_tkt;
emp_hat_tkt = oo_.irfs.emp_hat_tkt;
I_hat_tkt = oo_.irfs.I_hat_tkt;
k_hat_tkt = oo_.irfs.k_hat_tkt;
b_hat_y_tkt = oo_.irfs.b_hat_y_tkt;

deltatauk = [y_hat_tkt; c_hat_tkt; cr_hat_tkt; cnr_hat_tkt; w_hat_tkt; emp_hat_tkt; I_hat_tkt; k_hat_tkt; b_hat_y_tkt];
deltatauk = transpose(deltatauk);

% Figure  

figure
hold on
subplot(3, 3, 1)
plot(deltatauc(:,1), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,1), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,1), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Output')
xlabel('Quarters')
subplot(3, 3, 2)
plot(deltatauc(:,2), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,2), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,2), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Compsumption')
xlabel('Quarters')
subplot(3, 3, 3)
plot(deltatauc(:,3), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,3), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,3), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Compsumption R')
legend({'$\tau_c$','$\tau_w$','$\tau_k$'},'Location','best','Interpreter','latex'); %'Location','southeast'
xlabel('Quarters')
subplot(3, 3, 4)
plot(deltatauc(:,4), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,4), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,4), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Compsumption NR')
xlabel('Quarters')
subplot(3, 3, 5)
plot(deltatauc(:,5), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,5), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,5), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Wage')
xlabel('Quarters')
subplot(3, 3, 6)
plot(deltatauc(:,6), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,6), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,6), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Employment')
xlabel('Quarters')
subplot(3, 3, 7)
plot(deltatauc(:,7), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,7), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,7), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Investment')
xlabel('Quarters')
subplot(3, 3, 8)
plot(deltatauc(:,8), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,8), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,8), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Capital')
xlabel('Quarters')
subplot(3, 3, 9)
plot(deltatauc(:,9), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,9), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,9), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Debt')
xlabel('Quarters')
axis tight
hold off;

% Debt reduction
d_tct = cumsum([b_current,disc_vec.*b_hat_y_tct])./4;
d_twt = cumsum([b_current,disc_vec.*b_hat_y_twt])./4;
d_tkt = cumsum([b_current,disc_vec.*b_hat_y_tkt])./4;

figure
hold on
plot(d_tct, "LineWidth", 1, "Color", 'k')
hold on 
plot(d_twt, "LineWidth", 1, "Color", 'b')
hold on 
plot(d_tkt, "LineWidth", 1, "Color", 'r')
xlabel('Quarters')
legend({'$\tau_c$','$\tau_w$','$\tau_k$'},'Location','best','Interpreter','latex'); %'Location','southeast'
hold off;

%% ----------------------------------------------------------------------
%  1 pp reduction in Gov. Expend.
%% ----------------------------------------------------------------------

% 1 pp reduction in expend.
% g
y_hat_gt = oo_.irfs.y_hat_gt;
c_hat_gt = oo_.irfs.c_hat_gt;
cr_hat_gt = oo_.irfs.cr_hat_gt;
cnr_hat_gt = oo_.irfs.cnr_hat_gt;
w_hat_gt = oo_.irfs.w_hat_gt;
emp_hat_gt = oo_.irfs.emp_hat_gt;
I_hat_gt = oo_.irfs.I_hat_gt;
k_hat_gt = oo_.irfs.k_hat_gt;
b_hat_y_gt = oo_.irfs.b_hat_y_gt;

deltag1 = [y_hat_gt; c_hat_gt; cr_hat_gt; cnr_hat_gt; w_hat_gt; emp_hat_gt; I_hat_gt; k_hat_gt; b_hat_y_gt];
deltag1 = transpose(deltag1);

% Ig
y_hat_igt = oo_.irfs.y_hat_igt;
c_hat_igt = oo_.irfs.c_hat_igt;
cr_hat_igt = oo_.irfs.cr_hat_igt;
cnr_hat_igt = oo_.irfs.cnr_hat_igt;
w_hat_igt = oo_.irfs.w_hat_igt;
emp_hat_igt = oo_.irfs.emp_hat_igt;
I_hat_igt = oo_.irfs.I_hat_igt;
k_hat_igt = oo_.irfs.k_hat_igt;
b_hat_y_igt = oo_.irfs.b_hat_y_igt;

deltaig1 = [y_hat_igt; c_hat_igt; cr_hat_igt; cnr_hat_igt; w_hat_igt; emp_hat_igt; I_hat_igt; k_hat_igt; b_hat_y_igt];
deltaig1 = transpose(deltaig1);

% trans
y_hat_trt = oo_.irfs.y_hat_trt;
c_hat_trt = oo_.irfs.c_hat_trt;
cr_hat_trt = oo_.irfs.cr_hat_trt;
cnr_hat_trt = oo_.irfs.cnr_hat_trt;
w_hat_trt = oo_.irfs.w_hat_trt;
emp_hat_trt = oo_.irfs.emp_hat_trt;
I_hat_trt = oo_.irfs.I_hat_trt;
k_hat_trt = oo_.irfs.k_hat_trt;
b_hat_y_trt = oo_.irfs.b_hat_y_trt;

deltatr1 = [y_hat_trt; c_hat_trt; cr_hat_trt; cnr_hat_trt; w_hat_trt; emp_hat_trt; I_hat_trt; k_hat_trt; b_hat_y_trt];
deltatr1 = transpose(deltatr1);

% Figure  

figure
hold on
subplot(3, 3, 1)
plot(deltag1(:,1), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,1), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,1), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Output')
xlabel('Quarters')
subplot(3, 3, 2)
plot(deltag1(:,2), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,2), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,2), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Compsumption')
xlabel('Quarters')
subplot(3, 3, 3)
plot(deltag1(:,3), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,3), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,3), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Compsumption R')
legend({'$g$','$Ig$','$tr$'},'Location','best','Interpreter','latex'); %'Location','southeast'
xlabel('Quarters')
subplot(3, 3, 4)
plot(deltag1(:,4), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,4), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,4), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Compsumption NR')
xlabel('Quarters')
subplot(3, 3, 5)
plot(deltag1(:,5), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,5), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,5), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Wage')
xlabel('Quarters')
subplot(3, 3, 6)
plot(deltag1(:,6), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,6), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,6), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Employment')
xlabel('Quarters')
subplot(3, 3, 7)
plot(deltag1(:,7), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,7), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,7), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Investment')
xlabel('Quarters')
subplot(3, 3, 8)
plot(deltag1(:,8), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,8), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,8), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Capital')
xlabel('Quarters')
subplot(3, 3, 9)
plot(deltag1(:,9), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,9), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,9), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Debt')
xlabel('Quarters')
axis tight
hold off;

% Reduction in debt
d_gt = cumsum([b_current,disc_vec.*b_hat_y_gt])./4;
d_igt = cumsum([b_current,disc_vec.*b_hat_y_igt])./4;
d_trt = cumsum([b_current,disc_vec.*b_hat_y_trt])./4;

figure
hold on
plot(d_gt, "LineWidth", 1, "Color", 'k')
hold on 
plot(d_igt, "LineWidth", 1, "Color", 'b')
hold on 
plot(d_trt, "LineWidth", 1, "Color", 'r')
xlabel('Quarters')
legend({'$g$','$Ig$','$tr$'},'Location','best','Interpreter','latex'); %'Location','southeast'
hold off;


%%
% Permanent reduction of debt

clc
[sum(disc_vec*b_hat_gt)*100, sum(disc_vec*b_hat_igt)*100, sum(disc_vec*b_hat_trt)*100]
[sum(disc_vec*b_hat_tct)*100,sum(disc_vec*b_hat_twt)*100,sum(disc_vec*b_hat_tkt)*100]


%%

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%% Ramsey policy model
%M_.endo_names
%global M_ oo_ options_ ys0_ ex0_ estimation_info

clc
% Run Dynanre .mod code
%addpath C:/dynare/5.4/matlab/
dynare code_policy
%oo_.planner_objective_value.conditional.steady_initial_multiplier

%%
%%

% options_.ramsey.maxit = 10000;
% planner_objective beta*((cr_hat + (sigma_c/2)*(cr_hat)^2) - ((1-tao_wbar)/(1+tao_cbar)) * (N_hat+(sigma_n/2)*(N_hat)^2));
% ramsey_model(planner_discount=beta,instruments=(tao_c, tao_w, tao_k, g_hat, Ig_hat, trans)); %, tao_c, tao_w, tao_k, g_hat, Ig_hat, trans
% 
% compute consumption equivalent
% options_old=options_;
% options_.nocorr=1;
% options_.noprint=1;
% lambda_unconditional_shock=csolve('get_consumption_equivalent_unconditional_welfare',0,[],1e-8,1000,M_,oo_,options_)
% lambda_conditional_shock=csolve('get_consumption_equivalent_conditional_welfare',lambda_unconditional_shock,[],1e-8,1000,M_,oo_,options_)
% options_=options_old;
% 
% steady(solve_algo = 2, maxit = 1000000);
% // Simul. with perfect foresight
% perfect_foresight_setup(periods=500);
% perfect_foresight_solver;
% 
% evaluate_planner_objective;
 
c_pos= 2; %%
variance.c_hat=oo_.var(c_pos,c_pos);

display results
labels={'sigma(c_hat)';'L unc';'L cond'};
headers={'Utility shock'};
values_shock= [sqrt(variance.c_hat)];
options_.noprint=0;
dyntable(options_,[],headers,labels,100*values_shock,size(labels,2)+2,4,3)


%%
clc;
% G_t
disc_vec = beta.^[0:1:len-1];  % beta | disc
c_gt = (c_hat_gt + 0.5*(1-sigma_c)*(c_hat_gt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_gt +0.5*(1+sigma_n)*(N_hat_gt).^2);
cr_gt = (cr_hat_gt + 0.5*(1-sigma_c)*(cr_hat_gt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_gt +0.5*(1+sigma_n)*(N_hat_gt).^2);
cnr_gt = (cnr_hat_gt + 0.5*(1-sigma_c)*(cnr_hat_gt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_gt +0.5*(1+sigma_n)*(N_hat_gt).^2);
% Ig_t
c_igt = (c_hat_igt + 0.5*(1-sigma_c)*(c_hat_igt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_igt +0.5*(1+sigma_n)*(N_hat_igt).^2);
cr_igt = (cr_hat_igt + 0.5*(1-sigma_c)*(cr_hat_igt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_igt +0.5*(1+sigma_n)*(N_hat_igt).^2);
cnr_igt = (cnr_hat_igt + 0.5*(1-sigma_c)*(cnr_hat_igt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_igt +0.5*(1+sigma_n)*(N_hat_igt).^2);
%TR_t
c_tr = (c_hat_trt + 0.5*(1-sigma_c)*(c_hat_trt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_trt +0.5*(1+sigma_n)*(N_hat_trt).^2);
cr_tr = (cr_hat_trt + 0.5*(1-sigma_c)*(cr_hat_trt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_trt +0.5*(1+sigma_n)*(N_hat_trt).^2);
cnr_tr = (cnr_hat_trt + 0.5*(1-sigma_c)*(cnr_hat_trt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_trt +0.5*(1+sigma_n)*(N_hat_trt).^2);

% Tao_c
c_tct = (c_hat_tct + 0.5*(1-sigma_c)*(c_hat_tct).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_tct +0.5*(1+sigma_n)*(N_hat_tct).^2);
cr_tct = (cr_hat_tct + 0.5*(1-sigma_c)*(cr_hat_tct).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_tct +0.5*(1+sigma_n)*(N_hat_tct).^2);
cnr_tct = (cnr_hat_tct + 0.5*(1-sigma_c)*(cnr_hat_tct).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_tct +0.5*(1+sigma_n)*(N_hat_tct).^2);
% Tao_w
c_twt = (c_hat_twt + 0.5*(1-sigma_c)*(c_hat_twt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_twt +0.5*(1+sigma_n)*(N_hat_twt).^2);
cr_twt = (cr_hat_twt + 0.5*(1-sigma_c)*(cr_hat_twt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_twt +0.5*(1+sigma_n)*(N_hat_twt).^2);
cnr_twt = (cnr_hat_twt + 0.5*(1-sigma_c)*(cnr_hat_twt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_twt +0.5*(1+sigma_n)*(N_hat_twt).^2);
% Tao_k
c_tkt = (c_hat_tkt + 0.5*(1-sigma_c)*(c_hat_tkt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_tkt +0.5*(1+sigma_n)*(N_hat_tkt).^2);
cr_tkt = (cr_hat_tkt + 0.5*(1-sigma_c)*(cr_hat_tkt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_tkt +0.5*(1+sigma_n)*(N_hat_tkt).^2);
cnr_tkt = (cnr_hat_tkt + 0.5*(1-sigma_c)*(cnr_hat_tkt).^2) - (1-tao_wbar)/(1+tao_cbar)*(N_hat_tkt +0.5*(1+sigma_n)*(N_hat_tkt).^2);

% Expend.
[sum(disc_vec*c_gt),sum(disc_vec*cr_gt),sum(disc_vec*cnr_gt)]
[sum(disc_vec*c_igt),sum(disc_vec*cr_igt),sum(disc_vec*cnr_igt)]
[sum(disc_vec*c_tr),sum(disc_vec*cr_tr),sum(disc_vec*cnr_tr)]

% Taxes
[sum(disc_vec*c_tct),sum(disc_vec*cr_tct),sum(disc_vec*cnr_tct)]
[sum(disc_vec*c_twt),sum(disc_vec*cr_twt),sum(disc_vec*cnr_twt)]
[sum(disc_vec*c_tkt),sum(disc_vec*cr_tkt),sum(disc_vec*cnr_tkt)]

%% Ramsey analysis - Var decomp. 

c_pos = 2;
cr_pos = 3;
cnr_pos = 4;
n_pos = 11;

var_c = sqrt(oo_.var(c_pos,c_pos))*100;
var_cr = sqrt(oo_.var(cr_pos,cr_pos))*100;
var_cnr = sqrt(oo_.var(cnr_pos,cnr_pos))*100;
var_n = sqrt(oo_.var(n_pos,n_pos))*100;

[var_c, var_cr, var_cnr]

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% Oil shock

y_hat_oil = oo_.irfs.y_hat_oilt;
b_hat_oil = oo_.irfs.b_hat_oilt;
oil_shock = [y_hat_oil;b_hat_oil];

figure
subplot(1, 2, 1) 
hold on
plot(oil_shock(1,:), "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Output')
hold off
subplot(1, 2, 2)
hold on
plot(oil_shock(2,:), "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Debt')
hold off


%% For all variables
figure
subplot(3, 3, 1)
hold on
plot(y_hat_oilt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('GDP')
hold off
subplot(3, 3, 2)
hold on
plot(b_hat_oilt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Debt')
hold off
subplot(3, 3, 3)
hold on
plot(pic_hat_oilt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Inflation')
hold off
subplot(3, 3, 4)
hold on
plot(cr_hat_oilt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Compsumpt. R')
hold off
subplot(3, 3, 5)
hold on
plot(cnr_hat_oilt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Compsumpt. NR')
hold off
subplot(3, 3, 6)
hold on
plot(w_hat_oilt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Wage')
hold off
subplot(3, 3, 7)
hold on
plot(emp_hat_oilt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Employment')
hold off
subplot(3, 3, 8)
hold on
plot(I_hat_oilt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Investment')
hold off
subplot(3, 3, 9)
hold on
plot(k_hat_oilt, "LineWidth", 1, "Color", 'k')
yline(0,'--')
xlabel('Quarters')
title('Capital')
hold off


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%% Oil shocks

% - 1pp tao_c
y_hat_tct = oo_.irfs.y_hat_tct;
c_hat_tct = oo_.irfs.c_hat_tct;
cr_hat_tct = oo_.irfs.cr_hat_tct;
cnr_hat_tct = oo_.irfs.cnr_hat_tct;
w_hat_tct = oo_.irfs.w_hat_tct;
emp_hat_tct = oo_.irfs.emp_hat_tct;
I_hat_tct = oo_.irfs.I_hat_tct;
k_hat_tct = oo_.irfs.k_hat_tct;
b_hat_y_tct = oo_.irfs.b_hat_y_tct;

%% Baseline model - run code.mod
deltatauc = [y_hat_tct; c_hat_tct; cr_hat_tct; cnr_hat_tct; w_hat_tct; emp_hat_tct; I_hat_tct; k_hat_tct; b_hat_y_tct];
deltatauc = transpose(deltatauc);

%% Oil shock - run codeoil.mod
deltatauc1 = [y_hat_tct; c_hat_tct; cr_hat_tct; cnr_hat_tct; w_hat_tct; emp_hat_tct; I_hat_tct; k_hat_tct; b_hat_y_tct];
deltatauc1 = transpose(deltatauc1);
%%

% - 1pp tao_w
y_hat_twt = oo_.irfs.y_hat_twt;
c_hat_twt = oo_.irfs.c_hat_twt;
cr_hat_twt = oo_.irfs.cr_hat_twt;
cnr_hat_twt = oo_.irfs.cnr_hat_twt;
w_hat_twt = oo_.irfs.w_hat_twt;
emp_hat_twt = oo_.irfs.emp_hat_twt;
I_hat_twt = oo_.irfs.I_hat_twt;
k_hat_twt = oo_.irfs.k_hat_twt;
b_hat_y_twt = oo_.irfs.b_hat_y_twt;

%% Baseline model - run code.mod
deltatauw = [y_hat_twt; c_hat_twt; cr_hat_twt; cnr_hat_twt; w_hat_twt; emp_hat_twt; I_hat_twt; k_hat_twt; b_hat_y_twt];
deltatauw = transpose(deltatauw);
%% oil shock - run codeoil.mod
deltatauw1 = [y_hat_twt; c_hat_twt; cr_hat_twt; cnr_hat_twt; w_hat_twt; emp_hat_twt; I_hat_twt; k_hat_twt; b_hat_y_twt];
deltatauw1 = transpose(deltatauw1);
%%

% - 1pp tao_k
y_hat_tkt = oo_.irfs.y_hat_tkt;
c_hat_tkt = oo_.irfs.c_hat_tkt;
cr_hat_tkt = oo_.irfs.cr_hat_tkt;
cnr_hat_tkt = oo_.irfs.cnr_hat_tkt;
w_hat_tkt = oo_.irfs.w_hat_tkt;
emp_hat_tkt = oo_.irfs.emp_hat_tkt;
I_hat_tkt = oo_.irfs.I_hat_tkt;
k_hat_tkt = oo_.irfs.k_hat_tkt;
b_hat_y_tkt = oo_.irfs.b_hat_y_tkt;

%% Baseline model - run code.mod
deltatauk = [y_hat_tkt; c_hat_tkt; cr_hat_tkt; cnr_hat_tkt; w_hat_tkt; emp_hat_tkt; I_hat_tkt; k_hat_tkt; b_hat_y_tkt];
deltatauk = transpose(deltatauk);

%% oil shock - run codeoil.mod
deltatauk1 = [y_hat_tkt; c_hat_tkt; cr_hat_tkt; cnr_hat_tkt; w_hat_tkt; emp_hat_tkt; I_hat_tkt; k_hat_tkt; b_hat_y_tkt];
deltatauk1 = transpose(deltatauk1);
%% 
% Figure  

figure
hold on
subplot(3, 3, 1)
plot(deltatauc(:,1), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauc1(:,1), "--", "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,1), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauw1(:,1), "--", "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,1), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatauk1(:,1), "--", "LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('GDP')
xlabel('Quarters')
subplot(3, 3, 2)
plot(deltatauc(:,2), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauc1(:,2),"--", "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,2), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauw1(:,2),"--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,2), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatauk1(:,2),"--","LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Compsumption')
xlabel('Quarters')
subplot(3, 3, 3)
plot(deltatauc(:,3), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,3), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,3), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatauc1(:,3), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw1(:,3), "--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk1(:,3),"--", "LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Compsumpt. R')
legend({'$\tau_c$','$\tau_w$','$\tau_k$'},'Location','best','Interpreter','latex'); %'Location','southeast'
xlabel('Quarters')
subplot(3, 3, 4)
plot(deltatauc(:,4), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauc1(:,4), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,4), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauw1(:,4), "--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,4), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatauk1(:,4),"--" ,"LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Compsumpt. NR')
xlabel('Quarters')
subplot(3, 3, 5)
plot(deltatauc(:,5), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauc1(:,5),"--" ,"LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,5), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauw1(:,5),"--", "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,5), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatauk1(:,5), "--","LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Wage')
xlabel('Quarters')
subplot(3, 3, 6)
plot(deltatauc(:,6), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauc1(:,6), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,6), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauw1(:,6), "--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,6), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatauk1(:,6), "--","LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Employment')
xlabel('Quarters')
subplot(3, 3, 7)
plot(deltatauc(:,7), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauc1(:,7), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,7), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauw1(:,7), "--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,7), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatauk1(:,7), "--","LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Investment')
xlabel('Quarters')
subplot(3, 3, 8)
plot(deltatauc(:,8), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauc1(:,8), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,8), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauw1(:,8),"--", "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,8), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatauk1(:,8), "--","LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Capital')
xlabel('Quarters')
subplot(3, 3, 9)
plot(deltatauc(:,9), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauc1(:,9), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltatauw(:,9), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauw1(:,9),"--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatauk(:,9), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatauk1(:,9),"--", "LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Debt')
xlabel('Quarters')
axis tight
hold off;

%% ----------------------------------------------------------------------
%  1 pp reduction in Gov. Expend.
%% ----------------------------------------------------------------------

% -1 pp g
y_hat_gt = oo_.irfs.y_hat_gt;
c_hat_gt = oo_.irfs.c_hat_gt;
cr_hat_gt = oo_.irfs.cr_hat_gt;
cnr_hat_gt = oo_.irfs.cnr_hat_gt;
w_hat_gt = oo_.irfs.w_hat_gt;
emp_hat_gt = oo_.irfs.emp_hat_gt;
I_hat_gt = oo_.irfs.I_hat_gt;
k_hat_gt = oo_.irfs.k_hat_gt;
b_hat_y_gt = oo_.irfs.b_hat_y_gt;

%% Baseline model - run code.mod
deltag = [y_hat_gt; c_hat_gt; cr_hat_gt; cnr_hat_gt; w_hat_gt; emp_hat_gt; I_hat_gt; k_hat_gt; b_hat_y_gt];
deltag = transpose(deltag);
%% oil shock - run codeoil.mod
deltag1 = [y_hat_gt; c_hat_gt; cr_hat_gt; cnr_hat_gt; w_hat_gt; emp_hat_gt; I_hat_gt; k_hat_gt; b_hat_y_gt];
deltag1 = transpose(deltag1);
%%

%-1 pp Ig
y_hat_igt = oo_.irfs.y_hat_igt;
c_hat_igt = oo_.irfs.c_hat_igt;
cr_hat_igt = oo_.irfs.cr_hat_igt;
cnr_hat_igt = oo_.irfs.cnr_hat_igt;
w_hat_igt = oo_.irfs.w_hat_igt;
emp_hat_igt = oo_.irfs.emp_hat_igt;
I_hat_igt = oo_.irfs.I_hat_igt;
k_hat_igt = oo_.irfs.k_hat_igt;
b_hat_y_igt = oo_.irfs.b_hat_y_igt;

%% Baseline model - run code.mod
deltaig = [y_hat_igt; c_hat_igt; cr_hat_igt; cnr_hat_igt; w_hat_igt; emp_hat_igt; I_hat_igt; k_hat_igt; b_hat_y_igt];
deltaig = transpose(deltaig);
%% oil shock - run codeoil.mod
deltaig1 = [y_hat_igt; c_hat_igt; cr_hat_igt; cnr_hat_igt; w_hat_igt; emp_hat_igt; I_hat_igt; k_hat_igt; b_hat_y_igt];
deltaig1 = transpose(deltaig1);
%% 

% -1 pp tr
y_hat_trt = oo_.irfs.y_hat_trt;
c_hat_trt = oo_.irfs.c_hat_trt;
cr_hat_trt = oo_.irfs.cr_hat_trt;
cnr_hat_trt = oo_.irfs.cnr_hat_trt;
w_hat_trt = oo_.irfs.w_hat_trt;
emp_hat_trt = oo_.irfs.emp_hat_trt;
I_hat_trt = oo_.irfs.I_hat_trt;
k_hat_trt = oo_.irfs.k_hat_trt;
b_hat_y_trt = oo_.irfs.b_hat_y_trt;

%% Baseline model - run code.mod
deltatr = [y_hat_trt; c_hat_trt; cr_hat_trt; cnr_hat_trt; w_hat_trt; emp_hat_trt; I_hat_trt; k_hat_trt; b_hat_y_trt];
deltatr = transpose(deltatr);

%% oil shock - run codeoil.mod
deltatr1 = [y_hat_trt; c_hat_trt; cr_hat_trt; cnr_hat_trt; w_hat_trt; emp_hat_trt; I_hat_trt; k_hat_trt; b_hat_y_trt];
deltatr1 = transpose(deltatr1);
%%
% Figure  

figure
hold on
subplot(3, 3, 1)
plot(deltag(:,1), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltag1(:,1), "--", "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig(:,1), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltaig1(:,1), "--", "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr(:,1), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatr1(:,1), "--", "LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Output')
xlabel('Quarters')
subplot(3, 3, 2)
plot(deltag(:,2), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltag1(:,2),"--", "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig(:,2), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltaig1(:,2),"--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr(:,2), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatr1(:,2),"--","LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Compsumption')
xlabel('Quarters')
subplot(3, 3, 3)
plot(deltag(:,3), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig(:,3), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr(:,3), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltag1(:,3), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,3), "--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,3),"--", "LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Compsumt. R')
legend({'$g$','$Ig$','$tr$'},'Location','best','Interpreter','latex'); %'Location','southeast'
xlabel('Quarters')
subplot(3, 3, 4)
plot(deltag(:,4), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltag1(:,4), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig(:,4), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltaig1(:,4), "--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr(:,4), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatr1(:,4),"--" ,"LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Compsumt. NR')
xlabel('Quarters')
subplot(3, 3, 5)
plot(deltag(:,5), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltag1(:,5),"--" ,"LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig(:,5), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltaig1(:,5),"--", "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr(:,5), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatr1(:,5), "--","LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Wage')
xlabel('Quarters')
subplot(3, 3, 6)
plot(deltag(:,6), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltag1(:,6), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig(:,6), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltaig1(:,6), "--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr(:,6), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatr1(:,6), "--","LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Employment')
xlabel('Quarters')
subplot(3, 3, 7)
plot(deltag(:,7), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltag1(:,7), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig(:,7), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltaig1(:,7), "--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr(:,7), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatr1(:,7), "--","LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Investment')
xlabel('Quarters')
subplot(3, 3, 8)
plot(deltag(:,8), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltag1(:,8), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig(:,8), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltaig1(:,8),"--", "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr(:,8), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatr1(:,8), "--","LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Capital')
xlabel('Quarters')
subplot(3, 3, 9)
plot(deltag(:,9), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltag1(:,9), "--","LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig(:,9), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltaig1(:,9),"--","LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr(:,9), "LineWidth", 1, "Color", 'r')
hold on 
plot(deltatr1(:,9),"--", "LineWidth", 1, "Color", 'r')
%yline(0,'--')
title('Debt')
xlabel('Quarters')
axis tight
hold off;

%%
% Figure  

figure
hold on
subplot(3, 3, 1)
plot(deltag1(:,1), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,1), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,1), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Output')
xlabel('Quarters')
subplot(3, 3, 2)
plot(deltag1(:,2), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,2), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,2), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Compsumption')
xlabel('Quarters')
subplot(3, 3, 3)
plot(deltag1(:,3), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,3), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,3), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Compsumt. R')
xlabel('Quarters')
subplot(3, 3, 4)
plot(deltag1(:,4), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,4), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,4), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Compsumt. NR')
xlabel('Quarters')
subplot(3, 3, 5)
plot(deltag1(:,5), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,5), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,5), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Wage')
xlabel('Quarters')
subplot(3, 3, 6)
plot(deltag1(:,6), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,6), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,6), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Employment')
xlabel('Quarters')
subplot(3, 3, 7)
plot(deltag1(:,7), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,7), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,7), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Investment')
xlabel('Quarters')
subplot(3, 3, 8)
plot(deltag1(:,8), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,8), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,8), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Capital')
xlabel('Quarters')
subplot(3, 3, 9)
plot(deltag1(:,9), "LineWidth", 1, "Color", 'k')
hold on 
plot(deltaig1(:,9), "LineWidth", 1, "Color", 'b')
hold on 
plot(deltatr1(:,9), "LineWidth", 1, "Color", 'r')
yline(0,'--')
title('Debt')
xlabel('Quarters')
axis tight
hold off;


%%

% End -------------------

%%
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------