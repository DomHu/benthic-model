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

% run_and_plot_OMEN_with_GENIE_data(exp_1,'','Mean OM wt% in upper 10cm','',0.0,0.0,0,'',1.0,0.0,3.0,20,'','','',17, k_parametr, k2scaling)
run_and_plot_OMEN_with_GENIE_data(exp_1,'','oxygen penetration depth','',0.0,0.0,0,'',1.0,0.0,10.0,10,'','','',1, k_parametr, k2scaling)
run_and_plot_OMEN_with_GENIE_data(exp_1,'','frac_of_aerobic_Cox','',0.0,0.0,0,'',1.0,0.0,80.0,10,'','','',6, k_parametr, k2scaling)
run_and_plot_OMEN_with_GENIE_data(exp_1,'','Total_Cox_rate','',0.0,0.0,0,'',1e-6,0.0,400.0,10,'','','',3, k_parametr, k2scaling)
run_and_plot_OMEN_with_GENIE_data(exp_1,'','TOC wtpc at 100cm','',0.0,0.0,0,'',1.0,0.0,1.0,10,'','','',8, k_parametr, k2scaling)
end