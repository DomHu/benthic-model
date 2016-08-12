% Sensitivity study for OMEN-SED
%    
% http://www2.mae.ufl.edu/mdo/Papers/5176.pdf

% just change one Parameter e.h. k1

clear

% Make/Save Latin hypercube for parameters

% X = lhsdesign(n,p) returns an n-by-p matrix, X, 
% containing a latin hypercube sample of n values on each of p variables. 
% For each column of X, the n values are randomly distributed with one from each interval (0,1/n), (1/n,2/n), ..., (1-1/n,1), and they are randomly permuted.

% Parameters:
% k1	OM degradation frac 1 (labile)
% k2 = 1.0e-4 (as this is smallest for k1 in Sandras review study)	OM degradation frac 2 (refractory)
% KIPO4	P adsorption coeff. oxic
% KIIPO4	P adsorption coeff. anoxic
par = 1;    % parameters
n = 300;     % n values are randomly distributed with one from each interval (0,1/n), (1/n,2/n), ..., (1-1/n,1)
% calculate Latin Hypercube: each row is one representations of the variables for an Experiment

Latin_Cube_k_TOC_f = lhsdesign(n,par);
save Latin_Cube_k_TOC_f.mat Latin_Cube_k_TOC_f
%load('Latin_Cube.mat')

% Set the parameter Params.ranges
Params.range_k1 = [log10(1e-4), log10(20)];  % [-4, 1.3010] OM degradation frac 1 (labile) 
%Params.range_f1 = [0.05, 0.95];              % fraction of labile OM
%Params.range_wtpc = [log10(0.01), log10(20)];              % POC wtpc adsorption

Para1 = [Params.range_k1(1)]; %[Params.range_wtpc(1), Params.range_f1(1),Params.range_wtpc(1)];

% Calculate length of Intervals
Int = [diff(Params.range_k1)]; %diff(Params.range_wtpc), diff(Params.range_f1), 
% Use Latin hypercube and Int to calculate which values to chose; multiple element by element
WhichV = bsxfun(@times, Latin_Cube_k_TOC_f, Int);

% Now get actual values, add element by element
V = bsxfun(@plus, WhichV, Para1);
% unlog k1 and kaPO4

% p = [0.01:0.002:0.99];
% q = [0.1:0.1:20];
% Params.wtpc = [p q];
%Params.wtpc = 10.^V(:,1);
Params.k1 = 10.^V(:,1); %ones(1,length(Params.wtpc))/10; %
Params.f1 = 0.5;
Params.wtpc = ones(1,length(Params.k1));

% initialize SWI concentrationsn and other parameters
swi=benthic_test.default_swi();
% 

% % set date-time
str_date = '0307_Shallow_k1'; %datestr(now,'ddmmyy_HH_MM_SS');
swi=benthic_test.sensitivity_swi_singleParameter(swi, Params, str_date);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%                       PLOT the Results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_results_singlePar = true;
if(plot_results_singlePar)
%    Results_NaN = load('Results_NaN.txt','ascii');  % TODO: check how to create this file automaticly
% with created results    Plot_sensitivity_singleParameter(swi.Results, log10(Params.k1), Params.range_k1, 'k1', 'log(k1) [yr^{-1}]')
Plot_sensitivity_penetrationdepths_SWIflux(swi.Results, log10(Params.k1), Params.range_k1, 'k1', 'log(k1) [yr^{-1}]', str_date)
% Dom for TOC load     Plot_sensitivity_penetrationdepths_SWIflux(swi.Results, log10(Params.wtpc), [log10(0.01) log10(20)], 'TOC_load', 'log(TOC wtpc)', str_date)
   
%%    Plot_sensitivity_singleParameter(swi.Results, Params.f1, [0, 1], 'f1', 'labile fraction')
%%    Plot_sensitivity_singleParameter(swi.Results, log10(Params.wtpc), Params.range_wtpc,'TOCwtpc', 'TOC wtpc [%]')
end

plot_results_singleOut = false;
if(plot_results_singleOut)    
%    Results_NaN = load('Results_NaN.txt','ascii'); % TODO: check how to create this file automaticly

    Plot_sensitivity_singleOutput(Params, swi.Results, 2, 'FO2', 'F_{O_2}')
    Plot_sensitivity_singleOutput(Params, swi.Results, 3, 'FNO3', 'F_{NO_3}')
    Plot_sensitivity_singleOutput(Params, swi.Results, 4, 'FSO4', 'F_{SO_4}')
    Plot_sensitivity_singleOutput(Params, swi.Results, 5, 'FNH4', 'F_{NH_4}')
    Plot_sensitivity_singleOutput(Params, swi.Results, 6, 'FH2S', 'F_{H_2S}')
    Plot_sensitivity_singleOutput(Params, swi.Results, 7, 'FPO4', 'F_{PO_4}')
    Plot_sensitivity_singleOutput(Params, swi.Results, 8, 'zox', 'z_{O_2}')
    Plot_sensitivity_singleOutput(Params, swi.Results, 9, 'zNO3', 'z_{NO_3}')
    Plot_sensitivity_singleOutput(Params, swi.Results, 10, 'zSO4', 'z_{SO_4}')    
end

%%%%% Other techniques as proposed to use different ones: http://www2.mae.ufl.edu/mdo/Papers/5176.pdf

% % % % Generating Quasi-Random Numbers


% % % % http://uk.mathworks.com/help/stats/generating-quasi-random-rng default  % For reproducibility
% % % p = haltonset(4,'Skip',1e3,'Leap',1e2)      % Use haltonset to generate a 11-D Halton point set, skip the first 1000 values, and then retain every 101st point
% % % 
% % % p = scramble(p,'RR2')                       % Use scramble to apply reverse-radix scrambling.
% % % X0 = net(p,200);                            % Use net to generate the first 500 points.
% % % figure;
% % % scatter(X0(:,1),X0(:,2),5,'r')
% % % axis square
% % % title('{\bf Quasi-Random Scatter}')