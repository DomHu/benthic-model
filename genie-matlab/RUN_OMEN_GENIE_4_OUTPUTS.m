function RUN_OMEN_GENIE_4_OUTPUTS(PEXP1, PKPARAM, PK2SCALING)
% 
%   PKPARAM [STRING] (e.g., 'boudreau1997')
%   --> the string for the k parameterisation scheme
%   -->  options are 
%   -->  'boudreau1997'
%   -->  'tromp1995'
%   -->  'stolpovsky2016'
%   -->  'boudreau1997fPOC'
%   -->  'invariant' here, global values must be specified in
%         benthic_test.OMEN_with_GENIE_input(...)
%   --> PK2SCALING [REAL] scaling factor for calculating k2 = PK2SCALING * k1

% set passed parameters
exp_1 = PEXP1;
k_parametr = PKPARAM;
k2scaling = PK2SCALING;

% run_and_plot_OMEN_with_GENIE_data(exp_1,'','Mean OM wt% in upper 5cm','',0.0,0.0,0,'',1.0,0.0,3.0,20,'','','',17, k_parametr, k2scaling)

% run_and_plot_OMEN_with_GENIE_data(exp_1,'','oxygen penetration depth','',0.0,0.0,0,'',1.0,0.0,10.0,20,'','','',1, k_parametr, k2scaling)
% % %run_and_plot_OMEN_with_GENIE_data(exp_1,'','frac_of_aerobic_Cox','',0.0,0.0,0,'',1.0,0.0,80.0,20,'','','',6, k_parametr, k2scaling)
% run_and_plot_OMEN_with_GENIE_data(exp_1,'','frac_of_aerobic_Cox_upper_xcm','',0.0,0.0,0,'',1.0,0.0,80.0,20,'','','',19, k_parametr, k2scaling)
%run_and_plot_OMEN_with_GENIE_data(exp_1,'','Total_Cox_rat','',0.0,0.0,0,'',1e-6,0.0,100.0,20,'','','',3, k_parametr, k2scaling)
% run_and_plot_OMEN_with_GENIE_data(exp_1,'','Total_Cox_rate upper 5cm_till_100','',0.0,0.0,0,'',1e-6,0.0,100.0,20,'','','',18, k_parametr, k2scaling)
% % %run_and_plot_OMEN_with_GENIE_data(exp_1,'','TOC wtpc at 100cm','',0.0,0.0,0,'',1.0,0.0,1.0,10,'','','',8, k_parametr, k2scaling)
% run_and_plot_OMEN_with_GENIE_data(exp_1,'','O2_SWI-flux_until100','',0.0,0.0,0,'',1e-6,0.0,-100.0,20,'','','',11, k_parametr, k2scaling)
% run_and_plot_OMEN_with_GENIE_data(exp_1,'','SO4_SWI-flux_until5','',0.0,0.0,0,'',1e-6,0.0,-5.0,20,'','','',12, k_parametr, k2scaling)
% run_and_plot_OMEN_with_GENIE_data(exp_1,'','H2S_SWI-flux_until5','',0.0,0.0,0,'',1e-6,0.0,5.0,20,'','','',13, k_parametr, k2scaling)
run_and_plot_OMEN_with_GENIE_data(exp_1,'','PO4_SWI','',0.0,0.0,0,'',1e-9,-20.0,20.0,20,'','','',14, k_parametr, k2scaling)

% run_and_plot_OMEN_with_GENIE_data(exp_1,'','Mean_PO4_biot_layer','',0.0,0.0,0,'',1e-9,0.0,1,20,'','','',20, k_parametr, k2scaling)
%run_and_plot_OMEN_with_GENIE_data(exp_1,'','Mean_FeP_biot_layer','',0.0,0.0,0,'',1e-5,0.0,10,20,'','','',21, k_parametr, k2scaling)
end