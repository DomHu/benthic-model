function [y] = OMEN_SED_vbsa(x,res)
%
% This function runs the OMEN-SED sediment model
% and returns the associated SWI interface flux (TODO: include other
% variables as well)
%
% [y,Q_sim,STATES,FLUXES] = hymod_nse(param,rain,evap,flow)
% 
% Input:
% x = vector of model parameters ('k1','k2ord','f1','KNH4','gamma NH4','gamma H2S','K I','K II','ks' 'km' 'ka' )  - vector (1,11)
%
% Output:
%      y(1) = O2 SWI flux                     - scalar
%      y(2) = NO3 SWI flux                     - scalar
%      y(3) = SO4 SWI flux                     - scalar
%      y(4) = NH4 SWI flux                     - scalar
%      y(5) = H2S SWI flux                     - scalar
%      y(6)   = P SWI flux                     - scalar
%

% how much Corg wtpc at top of sedments:
wtpc = 2.0;

M = 11 ; % number of model parameters
x = x(:);
if ~isnumeric(x); error('input argument ''param'' must be numeric'); end
if length(x)~=M; error('input argument ''param'' must have %d components',M); end


res.swi.C01=x(3)*wtpc*1e-2/12*res.bsd.rho_sed;       %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
res.swi.C02=(1-x(3))*wtpc*1e-2/12*res.bsd.rho_sed;
res.zTOC.k1 = x(1);
res.zTOC.k2 = x(1)*x(2);                
% O2, NO3, SO4, NH4, H2S
res.zNO3.KNH4 = x(4);
res.zNH4.KNH4 = x(4);
res.bsd.gamma = x(5);
res.bsd.gammaH2S = x(6);

% PO4 use just anoxic 400m setup
res.zPO4_M.KPO41 = x(7);
res.zPO4_M.KPO42 = x(8);
res.zPO4_M.ksPO4 = x(9);
res.zPO4_M.kmPO4 = x(10);
res.zPO4_M.kaPO4 = x(11);


res = res.zTOC.calc(res.bsd,res.swi, res);
res = res.zO2.calc(res.bsd, res.swi, res);
%            if(swi.Nitrogen)
res = res.zNO3.calc(res.bsd, res.swi, res);
%            else
%            res.zno3=res.zox;
%            end
res = res.zSO4.calc(res.bsd, res.swi, res);
%            if(swi.Nitrogen)
res = res.zNH4.calc(res.bsd, res.swi, res);
%            end
res = res.zH2S.calc(res.bsd, res.swi, res);
res = res.zPO4_M.calc(res.bsd, res.swi, res);
% % res = res.zDIC.calc(res.bsd, res.swi, res);
% % res = res.zALK.calc(res.bsd, res.swi, res);

% % y = res.flxswiO2;
y(1) = res.flxswiO2;
y(2) = res.flxswiNO3; 
y(3) = res.flxswiSO4; 
y(4) = res.flxswiNH4; 
y(5) = res.flxswiH2S; 
y(6) = res.flxswi_P;

end