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
% The model under study is the Organic Matter ENabled SEDiment model
% 2G model to find in GitHub/benthic-model/Two_fractions
%
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
addpath([ my_dir '/../Two_fractions'])
addpath([ my_dir '/VBSA'])

%% Step 2: setup the model and define input ranges

% % % Load data:
% % load -ascii LeafCatch.txt
% % rain = LeafCatch(1:365,1)   ;
% % evap = LeafCatch(1:365,2)   ;
% % flow = LeafCatch(1:365,3)   ;

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
res.zH2S = benthic_zH2S(res.bsd, res.swi);
res.zPO4_M = benthic_zPO4_M(res.bsd, res.swi);


% Define input distribution and ranges:
M  = 5 ; % number of uncertain parameters [ k1 f1 KNH4 gammaNH4 gammaH2S ]
DistrFun  = {'unif', 'unif', 'unif', 'unif', 'unif'}; %'unif'  ; % Parameter distribution
DistrPar  = { [ 1e-4 5 ]; [ 0.05 0.95 ]; [ 0.8 1.7 ]; [ 0.5 1] ; [ 0.5 1 ] } ; % Parameter ranges
%DistrPar  = { [ log10(1e-4) log10(5) ]; [ 0.05 0.95 ]; [ 0.8 1.7 ]; [ 0.5 1] ; [ 0.5 1 ] } ; % Parameter ranges

%% Step 3: Compute first-order and total-order variance-based indices

myfun = 'OMEN_SED_vbsa' ;

% Sample parameter space using the resampling strategy proposed by 
% (Saltelli, 2008; for reference and more details, see help of functions
% vbsa_resampling and vbsa_indices) 
SampStrategy = 'lhs' ;
N = 1000 ; % Base sample size.
% Comment: the base sample size N is not the actual number of input 
% samples that will be evaluated. In fact, because of the resampling
% strategy, the total number of model evaluations to compute the two
% variance-based indices is equal to N*(M+2) 
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,2*N);

%X(:,1)=10.^X(:,1);
[ XA, XB, XC ] = vbsa_resampling(X) ;

% Run the model and compute selected model output at sampled parameter
% sets:
YA = model_evaluation(myfun,XA,res) ; % size (N,1)

if(false)
YB = model_evaluation(myfun,XB) ; % size (N,1)
YC = model_evaluation(myfun,XC) ; % size (N*M,1)

% Compute main (first-order) and total effects:
[ Si, STi ] = vbsa_indices(YA,YB,YC);

% Plot results:
X_Labels = {'Sm','beta','alfa','Rs','Rf'} ;
figure % plot main and total separately
subplot(121); boxplot1(Si,X_Labels,'main effects')
subplot(122); boxplot1(STi,X_Labels,'total effects')
figure % plot both in one plot:
boxplot2([Si; STi],X_Labels)
legend('main effects','total effects')% add legend

% Check the model output distribution (if multi-modal or highly skewed, the
% variance-based approach may not be adequate):
Y = [ YA; YC ] ;
figure; plot_cdf(Y,'NSE')
figure; plot_pdf(Y,'NSE');

% Compute confidence bounds:
Nboot = 500 ;
[ Si, STi, Si_sd, STi_sd, Si_lb, STi_lb, Si_ub, STi_ub ] = vbsa_indices(YA,YB,YC,Nboot);
% Plot:
figure % plot main and total separately
subplot(121); boxplot1(Si,X_Labels,'main effects',Si_lb,Si_ub)
subplot(122); boxplot1(STi,X_Labels,'total effects',STi_lb,STi_ub)
figure % plot both in one plot:
boxplot2([Si; STi],X_Labels,[ Si_lb; STi_lb ],[ Si_ub; STi_ub ])
legend('main effects','total effects')

% Analyze convergence of sensitivity indices:
NN = [N/10:N/10:N] ;
[ Sic, STic ] = vbsa_convergence([YA;YB;YC],M,NN);
figure
subplot(121); plot_convergence(Sic,NN*(M+2),[],[],[],'model evals','main effect',X_Labels)
subplot(122); plot_convergence(STic,NN*(M+2),[],[],[],'model evals','total effect',X_Labels);
% With confidence bounds:
[ Sic, STic, Si_sdc, STi_sdc, Si_lbc, STi_lbc, Si_ubc, STi_ubc  ] = vbsa_convergence([YA;YB;YC],M,NN,Nboot);
figure
subplot(121); plot_convergence(Sic,NN*(M+2),Si_lbc,Si_ubc,[],'model evals','main effect',X_Labels)
subplot(122); plot_convergence(STic,NN*(M+2),STi_lbc,STi_ubc,[],'model evals','total effect',X_Labels);

%% Step 4: Adding up new samples

N2 = 500 ; % increase of base sample size
% (that means: N2*(M+2) new samples that will need to be evaluated)

Xext = AAT_sampling_extend(X,DistrFun,DistrPar,2*(N+N2)) ; % extended sample 
% (it includes the already evaluated samples 'X' and the new ones
Xnew = Xext(2*N+1:end,:) ; % extract the new input samples that need to be evaluated
% Resampling strategy:
[ XA2, XB2,XC2 ] = vbsa_resampling(Xnew);
% Evaluate model against new samples:
YA2 =model_evaluation(myfun,XA2,rain,evap,flow); % should have size (N2,1)
YB2 =model_evaluation(myfun,XB2,rain,evap,flow); % should have size (N2,1)
YC2 =model_evaluation(myfun,XC2,rain,evap,flow); % should have size (N2*M,1)

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

%% Step 5: case of multiple outputs 
% (In this example: RMSE and AME)

myfun = 'hymod_MulObj' ;
YA = model_evaluation(myfun,XA,rain,evap,flow) ; % size (N,P)
YB = model_evaluation(myfun,XB,rain,evap,flow) ; % size (N,P)
YC = model_evaluation(myfun,XC,rain,evap,flow) ; % size (N*M,P)

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

