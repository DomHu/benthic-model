function RUN_OMEN_GENIE_4_OUTPUTS(PEXP1)
% set passed parameters
exp_1 = PEXP1;

run_and_plot_OMEN_with_GENIE_data(exp_1,'','oxygen penetration depth','',0.0,0.0,0,'',1.0,0.0,10.0,10,'','','',1)
run_and_plot_OMEN_with_GENIE_data(exp_1,'','frac_of_aerobic_Cox','',0.0,0.0,0,'',1.0,0.0,80.0,10,'','','',6)
run_and_plot_OMEN_with_GENIE_data(exp_1,'','Total_Cox_rate','',0.0,0.0,0,'',1e-6,0.0,400.0,10,'','','',3)
run_and_plot_OMEN_with_GENIE_data(exp_1,'','TOC wtpc at 100cm','',0.0,0.0,0,'',1.0,0.0,1.0,10,'','','',8)
end