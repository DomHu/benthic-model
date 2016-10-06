% This script provides an application example
% of the PAWN sensitivity analysis approach (Pianosi and Wagener, 2015)
%
% MODEL AND STUDY AREA
%
% The model under study is the rainfall-runoff model Hymod
% (see help of function hymod_sim.m for more details) 
% applied to the Leaf catchment in Mississipi, USA
% (see header of file LeafCatch.txt for more details).
% The inputs subject to SA are the 5 model parameters, and the scalar 
% output for SA is a statistic of the simulated time series
% (e.g. the maximum flow over the simulation horizon)
% 
% REFERENCES
%
% Pianosi, F. and Wagener, T. (2015), A simple and efficient method 
% for global sensitivity analysis based on cumulative distribution 
% functions, Env. Mod. & Soft., 67, 1-11.

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2015
% mail to: francesca.pianosi@bristol.ac.uk

%% Step 1: Add paths

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath([ my_dir '/sampling'])
%addpath([ my_dir '/example/hymod'])
addpath([ my_dir '/../Two_fractions'])
addpath([ my_dir '/visualization'])
addpath([ my_dir '/util'])
addpath([ my_dir '/PAWN'])

%% Step 2: setup the OMEN model

% Initialze the model
% initialize SWI concentrationsn and other parameters
swi=benthic_test.default_swi();

% Dominik: could give Corg wt% here...
% wtpc = 1.0;

res.bsd = benthic_main(1);
res.bsd.usescalarcode = true;

res.swi = swi;

% Set default values 
res.zTOC = benthic_zTOC(res.bsd);
res.zO2 = benthic_zO2(res.bsd, res.swi);           
res.zNO3 = benthic_zNO3(res.bsd, res.swi);
res.zSO4 = benthic_zSO4(res.bsd, res.swi);
res.zNH4 = benthic_zNH4(res.bsd, res.swi);
res.zH2S = benthic_zH2S(res.bsd, res.swi);
res.zPO4_M = benthic_zPO4_M(res.bsd, res.swi);


% Define uncertain inputs (parameters):

% O2, NO3, SO4, NH4, H2S
M = 11 ; % number of inputs
labelparams={ 'k1','k2ord','f1','KNH4','gamma NH4','gamma H2S','K I','K II','ks' 'km' 'ka' } ; % input names
xmin = [  1e-4 1e-6 0.02 0.8 0.5 0.5  100 1.3  0.1 0.015  0.001];
xmax = [    5  1e-1 0.98 1.7 1.0 1.0  400 2.0 100.0 0.02  10.0 ];

% PO4 
% M = 7 ; % number of inputs
% labelparams={ 'k1','f1','K I','K II','ks' 'km' 'ka'} ; % input names
% xmin = [  1e-4 0.02 100 1.3  50.0 0.015  0.001];
% xmax = [    5  0.98 400 2.0 850.0 0.02  10.0];
% % M = 5 ; % number of inputs
% % labelparams={ 'k1','f1','ks' 'km' 'ka'} ; % input names
% % xmin = [  1e-4 0.02  50.0 0.015  0.001];
% % xmax = [    5  0.98  850.0 0.02  10.0];

distrpar=cell(M,1); for i=1:M; distrpar{i}=[xmin(i) xmax(i)]; end

% Define model output:
fun_test = 'OMEN_SED_vbsa';
Titles = {'O_2', 'NO_3', 'SO_4', 'NH_4', 'H_2S', 'PO_4'};

%% Step 3: Apply PAWN

NU = 200 ; % number of samples to estimate unconditional CDF
NC = 100 ; % number of samples to estimate conditional CDFs
n  = 10 ; % number of conditioning points
out = 6; % number of output from OMEN (e.g. SWI-flux O2, NO3, ...)

% Create input/output samples to estimate the unconditional output CDF:
Xu = AAT_sampling('lhs',M,'unif',distrpar,NU); % matrix (NU,M)
Yu = model_evaluation(fun_test,Xu,res)  ; % matrix (NU,out) % was vector (1,M)
save('RESULTS_PAWN_ALL_OUTPUT/Xu_400m.mat','Xu')
save('RESULTS_PAWN_ALL_OUTPUT/Yu_400m.mat','Yu')
% Yu lines: input param sets
% Yu columns: outputs
%      y(:,1) = O2 SWI flux
%      y(:,2) = NO3 SWI flux
%      y(:,3) = SO4 SWI flux
%      y(:,4) = NH4 SWI flux
%      y(:,5) = H2S SWI flux
%      y(:,6)   = P SWI flux

% Create input/output samples to estimate the conditional output CDFs:
[ XX, xc ] = pawn_sampling('lhs',M,'unif',distrpar,n,NC);
save('RESULTS_PAWN_ALL_OUTPUT/XX_400m.mat','XX')
save('RESULTS_PAWN_ALL_OUTPUT/xc_400m.mat','xc')

YY = pawn_model_evaluation(fun_test,XX,res) ;
save('RESULTS_PAWN_ALL_OUTPUT/YY_400m.mat','YY')

for j=1:out
% Estimate unconditional and conditional CDFs:
[ YF, Fu, Fc  ] = pawn_cdfs(Yu,YY, j) ;


% % % Plot CDFs:
% % figure
% % for i=1:M
% %    subplot(1,M,i)
% %    pawn_plot_cdf(YF, Fu, Fc(i,:),[],'SWI-flux')
% % end
% % title(Titles(j))
% % 
% % % Further analyze CDF of one input:
% % i = 1 ;
% % figure;
% % pawn_plot_cdf(YF, Fu, Fc(i,:),xc{i},'SWI-flux',labelparams{i}) % same
% % title(Titles(j))
% % % function as before but exploiting more optional input arguments
% % 
% Compute KS statistics:
KS = pawn_ks(YF,Fu,Fc) ;

% Plot KS statistics:
figure
for i=1:M
   subplot(1,M,i)
   pawn_plot_kstest(KS(:,i),NC,NU,0.05,xc{i},labelparams{i})
end
title(Titles(j))
print('-depsc2', ['RESULTS_PAWN_ALL_OUTPUT/2_KS_400m_' char(Titles(j)) '.eps']);

% %  HERE without confidence intervals
% % % Compute PAWN index by taking a statistic of KSs (e.g. max):
% % Pi = max(KS);
% % % Plot:
% % figure 
% % boxplot1(Pi,labelparams)
% % title(Titles(j))

% Use bootstrapping to assess robustness of PAWN indices:
stat = 'max' ; % statistic to be applied to KSs
Nboot = 100  ; % number of boostrap resamples
[ T_m, T_lb, T_ub ] = pawn_indices(Yu,YY,stat,j,Nboot);
SIndex(j,:) = T_m;

% Plot:
figure; boxplot1(T_m,labelparams,[],T_lb,T_ub)
title(Titles(j))
print('-depsc2', ['RESULTS_PAWN_ALL_OUTPUT/1_SIndex_400m_' char(Titles(j)) '.eps']);

% Convergence analysis:
stat = 'max' ; % statistic to be applied to KSs
NCb = [ NC/10 NC/2 NC ] ;
NUb = [ NU/10 NU/2 NU ] ;

[ T_m_n, T_lb_n, T_ub_n ] = pawn_convergence( Yu, YY, stat, NUb, NCb,j,Nboot );
NN = NUb+n*NCb ;
figure; plot_convergence(T_m_n,NN,T_lb_n,T_ub_n,[],'no of evals',[],labelparams)
title(Titles(j))
print('-depsc2', ['RESULTS_PAWN_ALL_OUTPUT/3_Conv_400m_' char(Titles(j)) '.eps']);

% % 
% % %% Step 4: Apply PAWN to sub-region of the output range
% % 
% % % Compute the PAWN index over a sub-range of the output distribution, for
% % % instance only output values above a given threshold
% % thres = 50 ;
% % [ T_m2, T_lb2, T_ub2 ]= pawn_indices( Yu, YY, stat,[], Nboot,[],'above',thres ) ;
% % 
% % % Plot:
% % figure; boxplot1(T_m2,labelparams,[],T_lb2,T_ub2)

if(false)
    % % %% Step 4: create coloured scatter plot
    % k1 vs f1 -> flux 
    figure
    scatter_plots_col(Xu,Yu(:,j),1,2,16,labelparams)
    title(Titles(j))
    print('-depsc2', ['RESULTS_PAWN_ALL_OUTPUT/4000m/k1_vs_f1_SWIflux_' char(Titles(j)) '.eps']);

    %scatter_plots_col(XD,YC,1,2,16,X_Labels)

    % k1 vs SWI flux 
    figure
    scatter_plots(Xu(:,1),Yu(:,j),1,'SWI fluxes',{'k1'})
    title(Titles(j))
    print('-depsc2', ['RESULTS_PAWN_ALL_OUTPUT/4000m/k1_SWIflux_' char(Titles(j)) '.eps']);
end

end
SIndex400m = SIndex;
save('RESULTS_PAWN_ALL_OUTPUT/SIndex_400m.mat','SIndex400m')
