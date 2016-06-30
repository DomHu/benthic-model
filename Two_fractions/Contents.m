% BENTHIC-MODEL-2308-SD20150406
%
% Examples:
%
% Calculate a single sediment column with default parameters
%
% Execute now with 
% >> benthic_test.run_OMEN
% was:
% >> swi=benthic_test.default_swi()
% >> res=benthic_test.test_benthic(1,swi);
% >> benthic_test.plot_column(res, false)
%


% Calculate a single sediment column with default parameters and compare
% with BRNS results saved in OMEN-BRNS - save output *.ps in OMEN-BRNS
% >> benthic_test.run_OMEN_BRNS


% OMEN_sensitivity.m
% generate Latin Hypercube for various parameters and calculate SWI fluxes
% and penetration depths with OMEN 
% save resulting *.ps in folder Sensitivity






%%%%% STUART

% Calculate and plot 1000 sed columns vs O2 gradient 
% >> res=benthic_test.test_benthic(1000,swi);
% >> benthic_test.plot_summary(res)
%
% Files
%   benthic_main  - Global properties for benthic model
%   benthic_test  - test cases for benthic layer model
%   benthic_utils - Utility functions for solution matching etc
%   benthic_zH2S  - Solve SO4
%   benthic_zNH4  - Solve NH4
%   benthic_zNO3  - Solve NO3
%   benthic_zO2   - Solve O2
%   benthic_zSO4  - Solve SO4
%   benthic_zTOC  - Organic matter- 2 Fractions
%   fzero_vec     - Vectorized single-variable nonlinear zero finding based on fzero
