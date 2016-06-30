% Sensitivity study for OMEN-SED
%    
% http://www2.mae.ufl.edu/mdo/Papers/5176.pdf

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
par = 10;    % parameters
n = 200;     % n values are randomly distributed with one from each interval (0,1/n), (1/n,2/n), ..., (1-1/n,1)
% calculate Latin Hypercube: each row is one representations of the variables for an Experiment

%Latin_Cube = lhsdesign(n,par);
%save Latin_Cube.mat Latin_Cube
load('Latin_Cube.mat')

% Set the parameter Params.ranges
Params.range_k1 = [log10(1e-4), log10(20)];  % [-4, 1.3010] OM degradation frac 1 (labile) 
Params.range_f1 = [0.05, 0.95];              % fraction of labile OM
Params.range_KNH4 = [0.8, 1.7];              % NH4 adsorption
Params.range_KPO4ox = [100.0, 400.0];        % P Adsorption coefficient in oxic layer
Params.range_KPO4anox = [1.3, 2.0];          % P Adsorption coefficient in anoxic layer
Params.range_ksPO4 = [40.0, 900.0];          % Rate constant for kinetic P sorption
Params.range_kmPO4 = [0.015, 0.02];          % Rate constant for Fe-bound P release
Params.range_kaPO4 = [log10(0.001), log10(10.0)];          % [0.001, 10.0] Rate constant for authigenic P formation
Params.range_gammaNH4 = [0.5, 1];            % fraction of NH4 that is oxidised in oxic layer
Params.range_gammaH2S = [0.5, 1];            % fraction of H2S that is oxidised in oxic layer

Para1 = [Params.range_k1(1), Params.range_f1(1),Params.range_KNH4(1), Params.range_KPO4ox(1), Params.range_KPO4anox(1), Params.range_ksPO4(1), Params.range_kmPO4(1), Params.range_kaPO4(1), Params.range_gammaNH4(1), Params.range_gammaH2S(1)];

% Calculate length of Intervals
Int = [diff(Params.range_k1), diff(Params.range_f1), diff(Params.range_KNH4), diff(Params.range_KPO4ox), diff(Params.range_KPO4anox), diff(Params.range_ksPO4), diff(Params.range_kmPO4), diff(Params.range_kaPO4), diff(Params.range_gammaNH4), diff(Params.range_gammaH2S)]; 
% Use Latin hypercube and Int to calculate which values to chose; multiple element by element
WhichV = bsxfun(@times, Latin_Cube, Int);

% Now get actual values, add element by element
V = bsxfun(@plus, WhichV, Para1);
% unlog k1 and kaPO4
Params.k1 = 10.^V(:,1);
Params.f1 = V(:,2);
Params.KNH4 = V(:,3);
Params.KPO4ox = V(:,4); 
Params.KPO4anox= V(:,5);
Params.ksPO4 = V(:,6);
Params.kmPO4 = V(:,7);
Params.kaPO4 = 10.^V(:,8);
Params.gammaNH4 = V(:,9);
Params.gammaH2S = V(:,10);

% initialize SWI concentrationsn and other parameters
swi=benthic_test.default_swi();
% 
Params.wtpc = ones(1,n)./10;     % Overall POC wt\% reaching the SWI

swi=benthic_test.sensitivity_swi(swi, Params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%                       PLOT the Results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_results_singlePar = false;
if(plot_results_singlePar)
    Results_NaN = load('Results_NaN.txt','ascii');  % TODO: check how to create this file automaticly
% with created results    Plot_sensitivity_singleParameter(swi.Results, log10(Params.k1), Params.range_k1, 'k1', 'log(k1) [yr^{-1}]')
    Plot_sensitivity_singleParameter(Results_NaN, log10(Params.k1), Params.range_k1, 'k1', 'log(k1) [yr^{-1}]')
    Plot_sensitivity_singleParameter(swi.Results, Params.f1, [0, 1], 'f1', 'labile fraction')
    Plot_sensitivity_singleParameter(swi.Results, Params.KNH4, Params.range_KNH4,'KNH4', 'K_{NH_4} [-]')
    Plot_sensitivity_singleParameter(swi.Results, Params.KPO4ox, Params.range_KPO4ox,'KPO4ox', 'KPO4ox [-]')
    Plot_sensitivity_singleParameter(swi.Results, Params.KPO4anox, Params.range_KPO4anox,'KPO4anox', 'KPO4anox [-]')
    Plot_sensitivity_singleParameter(swi.Results, Params.ksPO4, Params.range_ksPO4,'ksPO4', 'ksPO4 [yr^{-1}]')
    Plot_sensitivity_singleParameter(swi.Results, Params.kmPO4, Params.range_kmPO4,'kmPO4', 'kmPO4 [yr^{-1}]')
    Plot_sensitivity_singleParameter(swi.Results, log10(Params.kaPO4), Params.range_kaPO4,'kaPO4', 'log(kaPO4) [yr^{-1}]')
    Plot_sensitivity_singleParameter(swi.Results, Params.gammaNH4, Params.range_gammaNH4,'gammaNH4', 'gammaNH4 [-]')
    Plot_sensitivity_singleParameter(swi.Results, Params.gammaH2S, Params.range_gammaH2S,'gammaH2S', 'gammaH2S [-]')
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