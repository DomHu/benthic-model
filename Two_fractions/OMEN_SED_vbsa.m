function [y] = OMEN_SED_vbsa(x,res)
%
% This function runs the OMEN-SED sediment model
% and returns the associated SWI interface flux (TODO: include other
% variables as well)
%
% [y,Q_sim,STATES,FLUXES] = hymod_nse(param,rain,evap,flow)
% 
% Input:
% x = vector of model parameters (k1, f1,  KNH4, gammaNH4, gammaH2S)  - vector (1,5)
%
% Output:
%      y = O2 SWI flux                     - scalar
%

% how much Corg wtpc at top of sedments:
wtpc = 10.0;

M = 5 ; % number of model parameters
x = x(:);
if ~isnumeric(x); error('input argument ''param'' must be numeric'); end
if length(x)~=M; error('input argument ''param'' must have %d components',M); end


res.swi.C01=x(2)*wtpc*1e-2/12*res.bsd.rho_sed;       %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
res.swi.C02=(1-x(2))*wtpc*1e-2/12*res.bsd.rho_sed;
res.zTOC.k1 = x(1);
res.zTOC.k2 = x(1)*0.01;                
res.zNO3.KNH4 = x(3);
res.bsd.gamma = x(4);
res.bsd.gammaH2S = x(5);


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

y = res.flxswiO2

end