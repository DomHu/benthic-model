% This script provides an application example of Variance Based Sensitivity
% Analysis (VBSA)
%
% METHODS
%
% We use two well established variance-based sensitivity indices: 
% - the first-order sensitivity index (or 'main effects')
% - the total-order sensitivity index (or 'total effects')
% (see help of 'vbsa_indices' for more details and references)
%
% MODEL AND STUDY AREA
%
% The routine under study is the Organic Matter burial scheme by Sandra
% Arndt to find in GitHub/Burial_Sandra
% INDEX
%
% Steps:
% 1. Add paths to required directories
% 2. Load data, set-up the Hymod model and define input ranges
% 3. Compute first-order (main effects) and total-order (total effects)
%    variance-based indices.
% 4: Example of how to repeat computions after adding up new 
%    input/output samples.
% 5. Example of how to compute indices when dealing with multiple outputs. 

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2014
% mail to: francesca.pianosi@bristol.ac.uk


clear

%% Step 1: set paths

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath([ my_dir '/sampling'])
addpath([ my_dir '/util'])
addpath([ my_dir '/visualization'])
%addpath([ my_dir '/example/hymod'])
addpath([ my_dir '/../../Burial_Sandra'])
addpath([ my_dir '/VBSA'])

%% Step 2: setup the model and define input ranges

% Define input distribution and ranges:

M  = 7 ; % number of uncertain parameters [ a nu w0 Db0 OC0 beta por0 ]
DistrFun  = {'unif', 'unif', 'unif', 'unif', 'unif', 'unif', 'unif'}; %'unif'  ; % Parameter distribution
DistrPar  = { [1e-1 1e5]; [0.1 2]; [1e-4 10];  [0.1 30]; [0.01 30]; [1e-7 1e-5]; [0.5 0.9]} ; % Parameter ranges

%% Step 3: Compute first-order and total-order variance-based indices

myfun = 'burial_Sandra' ;

% Sample parameter space using the resampling strategy proposed by 
% (Saltelli, 2008; for reference and more details, see help of functions
% vbsa_resampling and vbsa_indices) 
SampStrategy = 'lhs' ;
N = 2000 ; % Base sample size.
% Comment: the base sample size N is not the actual number of input 
% samples that will be evaluated. In fact, because of the resampling
% strategy, the total number of model evaluations to compute the two
% variance-based indices is equal to N*(M+2) 
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,2*N);

%X(:,1)=10.^X(:,1); % other Parameter range

[ XA, XB, XC ] = vbsa_resampling(X) ;
%save('RESULTS/X_5000.mat','X')
%save('RESULTS/XA_5000.mat','XA')
%save('RESULTS/XB_5000.mat','XB')
%save('RESULTS/XC_5000.mat','XC')


% Run the model and compute selected model output at sampled parameter
% sets: 
% lines: input param sets
% columns: outputs (see below)
YA_all = model_evaluation(myfun,XA); % size (N,1)
YB_all = model_evaluation(myfun,XB); % size (N,1)
YC_all = model_evaluation(myfun,XC); % size (N*M,1)
%save('RESULTS/YA_all_5000.mat','YA_all')
%save('RESULTS/YB_all_5000.mat','YB_all')
%save('RESULTS/YC_all_5000.mat','YC_all')

% select the j-th model output:
%      y(:,1) = O2 SWI flux
%      y(:,2) = NO3 SWI flux
%      y(:,3) = SO4 SWI flux
%      y(:,4) = NH4 SWI flux
%      y(:,5) = H2S SWI flux
%      y(:,6)   = P SWI flux

for j=8:16; 
YA = YA_all(:,j);
YB = YB_all(:,j);
YC = YC_all(:,j);

% Compute main (first-order) and total effects:
%was [ Si, STi ] = vbsa_indices(YA,YB,YC);
[ Si_tmp, STi_tmp ] = vbsa_indices(YA,YB,YC);
Si(j,:)=Si_tmp;
STi(j,:)=STi_tmp;
% Plot results:
X_Labels = {'a', 'nu', 'w0', 'Db0', 'OC0', 'beta', 'por0'} ;
Titles = {'F1', 'FF10', 'F2', 'FF100', 'F3', 'FF1000', 'F4', 'FF10000', 'dF1', 'dFF10', 'dF2', 'dFF100', 'dF3', 'dFF1000', 'dF4', 'dFF10000'};

% % figure % plot main and total separately
% % subplot(121); boxplot1(Si,X_Labels,'main effects')
% % subplot(122); boxplot1(STi,X_Labels,'total effects')

if(false)
    figure % plot both in one plot:
    boxplot2([Si_tmp; STi_tmp],X_Labels)
    legend('main effects','total effects')% add legend
    title(Titles(j))

    % Check the model output distribution (if multi-modal or highly skewed, the
    % variance-based approach may not be adequate):
    Y = [ YA; YC ] ;
    figure; plot_cdf(Y,'NSE') ;
    title(Titles(j))
end    
    % Check the model output distribution (if multi-modal or highly skewed, the
    % variance-based approach may not be adequate):
    Y = [ YA; YC ] ;
    figure; plot_pdf(Y,'NSE') ;
    title(Titles(j));
    print('-depsc2', ['../../Burial_Sandra/RESULTS/PDF_' char(Titles(j)) '.eps']);

% Compute confidence bounds:
Nboot = 500 ;
[ Si_tmp, STi_tmp, Si_sd_tmp, STi_sd_tmp, Si_lb_tmp, STi_lb_tmp, Si_ub_tmp, STi_ub_tmp ] = vbsa_indices(YA,YB,YC,Nboot);
Si(j,:)=Si_tmp;
STi(j,:)=STi_tmp;
Si_sd(j,:)=Si_sd_tmp;
STi_sd(j,:)=STi_sd_tmp;
Si_lb(j,:)=Si_lb_tmp;
STi_lb(j,:)=STi_lb_tmp;
Si_ub(j,:)=Si_ub_tmp;
STi_ub(j,:)=STi_ub_tmp;
% Plot:
% % figure % plot main and total separately
% % subplot(121); boxplot1(Si,X_Labels,'main effects',Si_lb,Si_ub)
% % subplot(122); boxplot1(STi,X_Labels,'total effects',STi_lb,STi_ub)
figure % plot both in one plot:
boxplot2([Si_tmp; STi_tmp],X_Labels,[ Si_lb_tmp; STi_lb_tmp ],[ Si_ub_tmp; STi_ub_tmp ])
legend('main effects','total effects')
title(Titles(j))
print('-depsc2', ['../../Burial_Sandra/RESULTS/SIndex_' char(Titles(j)) '.eps']);

if(false)
    % Analyze convergence of sensitivity indices:
    NN = [N/10:N/10:N] ;
    [ Sic, STic ] = vbsa_convergence([YA;YB;YC],M,NN);
    % % figure
    % % subplot(121); plot_convergence(Sic,NN*(M+2),[],[],[],'model evals','main effect',X_Labels)
    % % subplot(122); plot_convergence(STic,NN*(M+2),[],[],[],'model evals','total effect',X_Labels);
    % With confidence bounds:
    [ Sic, STic, Si_sdc, STi_sdc, Si_lbc, STi_lbc, Si_ubc, STi_ubc  ] = vbsa_convergence([YA;YB;YC],M,NN,Nboot);
    figure
    subplot(121); plot_convergence(Sic,NN*(M+2),Si_lbc,Si_ubc,[],'model evals','main effect',X_Labels)
    subplot(122); plot_convergence(STic,NN*(M+2),STi_lbc,STi_ubc,[],'model evals','total effect',X_Labels);
    title(X_Labels(j))
end

%% Step 4: create coloured scatter plot
if(false)
% % XD=XC;
% % XD(:,1)=log10(XC(:,1));
% k1 vs f1 -> flux 
figure
scatter_plots_col(XC,YC_all(:,j),1,2,16,X_Labels)
title(Titles(j))
print('-depsc2', ['RESULTS/uniform_k1_5_4000m/k1_vs_f1_SWIflux_' char(Titles(j)) '.eps']);

%scatter_plots_col(XD,YC,1,2,16,X_Labels)

% k1 vs SWI flux 
figure
scatter_plots(XC(:,1),YC_all(:,j),1,'SWI fluxes',{'k1'})
title(Titles(j))
print('-depsc2', ['RESULTS/uniform_k1_5_4000m/k1_SWIflux_' char(Titles(j)) '.eps']);
end
end

%save('RESULTS/STi_all_5000.mat','STi')

if(false)
%% Step 5: Adding up new samples

N2 = 500 ; % increase of base sample size
% (that means: N2*(M+2) new samples that will need to be evaluated)

Xext = AAT_sampling_extend(X,DistrFun,DistrPar,2*(N+N2)) ; % extended sample 
% (it includes the already evaluated samples 'X' and the new ones
Xnew = Xext(2*N+1:end,:) ; % extract the new input samples that need to be evaluated
% Resampling strategy:
[ XA2, XB2,XC2 ] = vbsa_resampling(Xnew);
% Evaluate model against new samples:
YA2 =model_evaluation(myfun,XA2); % should have size (N2,1)
YB2 =model_evaluation(myfun,XB2); % should have size (N2,1)
YC2 =model_evaluation(myfun,XC2); % should have size (N2*M,1)

% Put new and old results toghether:
YAn = [ YA; YA2 ] ; % should have size (N+N2,1)
YBn = [ YB; YB2 ] ; % should have size (N+N2,1)
YCn = [ reshape(YC,N,M); reshape(YC2,N2,M) ] ; % should have size (N+N2,M)
YCn = YCn(:); % should have size ((N+N2)*M,1)

% Recompute indices:
Nboot = 500 ;
[ Sin, STin, Si_sdn, STi_sdn, Si_lbn, STi_lbn, Si_ubn, STi_ubn ] = vbsa_indices(YAn,YBn,YCn,Nboot);

figure
subplot(121)
boxplot2([Si; STi],X_Labels,[ Si_lb; STi_lb ],[ Si_ub; STi_ub ]); title([ num2str(N*(M+2)) ' model eval.'])
subplot(122)
boxplot2([Sin; STin],X_Labels,[ Si_lbn; STi_lbn ],[ Si_ubn; STi_ubn ]); title([ num2str((N+N2)*(M+2)) ' model eval.'])
legend('main effects','total effects')

%% Step 6: case of multiple outputs 
% (In this example: RMSE and AME)


myfun = '' ;
YA = model_evaluation(myfun,XA) ; % size (N,P)
YB = model_evaluation(myfun,XB) ; % size (N,P)
YC = model_evaluation(myfun,XC) ; % size (N*M,P)

% select the j-th model output:
j = 1 ; 
[ Si1, STi1 ] = vbsa_indices(YA(:,j),YB(:,j),YC(:,j));
j = 2 ;
[ Si2, STi2 ] = vbsa_indices(YA(:,j),YB(:,j),YC(:,j));

% Compare boxplots:
figure
subplot(121)
boxplot2([Si1; STi1],X_Labels)
title('RMSE')
subplot(122)
boxplot2([Si2; STi2],X_Labels)
legend('main effects','total effects')
title('BIAS')

% Use stacked bar to put all outputs on one plot:
figure
subplot(121)
stackedbar([Si1; Si2],[],'main effects',[],{'RMSE','BIAS'})
subplot(122)
stackedbar([STi1; STi2],X_Labels,'total effects',[],{'RMSE','BIAS'})

% If you want to add samples in this case:
N2 = 500 ; % increase of base sample size (see previous Step)
Xext = AAT_sampling_extend(X,DistrFun,DistrPar,2*(N+N2)) ; % extended sample 
Xnew = Xext(2*N+1:end,:) ; % extract the new input samples that need to be evaluated
% Resampling strategy:
[ XA2, XB2,XC2 ] = vbsa_resampling(Xnew);
% Evaluate model against new samples:
YA2 =model_evaluation(myfun,XA2,rain,evap,flow); % should have size (N2,2)
YB2 =model_evaluation(myfun,XB2,rain,evap,flow); % should have size (N2,2)
YC2 =model_evaluation(myfun,XC2,rain,evap,flow); % should have size (N2*M,2)
% Select the j-th model output:
j = 1 ; 
% Put new and old results together:
YAn = [ YA(:,j); YA2(:,j) ] ; % should have size (N+N2,1)
YBn = [ YB(:,j); YB2(:,j) ] ; % should have size (N+N2,1)
YCn = [ reshape(YC(:,j),N,M); reshape(YC2(:,j),N2,M) ] ; % should have size (N+N2,M)
YCn = YCn(:); % should have size ((N+N2)*M,1)

[ Si1n, STi1n ] = vbsa_indices(YAn,YBn,YCn);
figure
subplot(121)
boxplot2([Si1; STi1],X_Labels); title([ num2str(N*(M+2)) ' model eval.'])
subplot(122)
boxplot2([Si1n; STi1n],X_Labels); title([ num2str((N+N2)*(M+2)) ' model eval.'])


end

