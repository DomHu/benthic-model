% OMEN-SED 1.0 BENTHIC-MODEL Stand-alone matlab code
% Hülse et al (2017) GMD paper

% benthic_test.m
% functions to run, plot and evaluate the model

% Command to run the model: benthic_test.run_OMEN
% 1) default sediment-water interface boundary conditions are set as
%    prescribed in default_swi()
% 2) the subroutines for the different properties are called in test_benthic(1,swi):
% 3) results saved in res and are plotted with plot_column()

classdef benthic_test
    % test cases for benthic layer model
    
    properties
    end
    
    methods(Static)
        
        function swi = default_swi()
            % set default SWI conditions
            
            format longEng
            
            bsd = benthic_main();
            %bottom water concentrations
            swi.T = 8.0;                                        % temperature (degree C)
            swi.C01_nonbio= 1.0*1e-2/12*bsd.rho_sed;            % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            swi.C02_nonbio= 1.0*1e-2/12*bsd.rho_sed;            % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            swi.Fnonbio1 = swi.C01_nonbio*(1-bsd.por)*bsd.w;    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
            swi.Fnonbio2 = swi.C02_nonbio*(1-bsd.por)*bsd.w;    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
            swi.C01 = swi.C01_nonbio;                           % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
            swi.C02 = swi.C02_nonbio;                           % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
            swi.O20=300.0E-009;                                 % O2  concentration at SWI (mol/cm^3)
            swi.NO30=40.0e-9;                                   % NO3 concentration at SWI (mol/cm^3)
            swi.Nitrogen=true;                                  % calculate N (true/false)
            swi.NH40=10.0e-9;                                 	% NH4 concentration at SWI (mol/cm^3)
            swi.SO40=2.8E-005;                                	% SO4 concentration at SWI (mol/cm^3)
            swi.H2S0=0.0;                                       % H2S concentration at SWI (mol/cm^3)
            swi.PO40=40.0e-9;                                   % PO4 concentration at SWI (mol/cm^3)
            swi.Mflux0=365*0.2e-10*1/(1-bsd.por)*1/bsd.w;       % actually CONCENTRATION of M at the sediment [mol/cm3] : from flux input  365*0.2e-10 (mol/(cm2*yr))
            swi.DIC0=2.4E-006;                                 	% DIC concentration at SWI (mol/cm^3)
            swi.ALK0=2.4E-006;                                 	% ALK concentration at SWI (mol/cm^3)
            swi.S0=35;                                         	% Salinity at SWI (not used at the moment)
            swi.plot_PO4_DIC_ALK=true;
        end
        
        function run_OMEN()
            % run OMEN-SED with default SWI conditions as in default_swi()
            clear
            tic;
            swi=benthic_test.default_swi()
            %            % set date-time
            %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
            res=benthic_test.test_benthic(1,swi);
            toc;
            benthic_test.plot_column(res, false, swi, 'FULL_OMEN')
            
            % calculate depth integrated OM degradation rates
            Cox_rate.Cox_total = res.zTOC.calcReac(0.0, res.bsd.zinf, 1, 1, res.bsd, swi, res);
            Cox_rate.Cox_aerobic = res.zTOC.calcReac(0.0, res.zox, 1, 1, res.bsd, swi, res);
            if(swi.Nitrogen)
                Cox_rate.Cox_denitr = res.zTOC.calcReac(res.zox, res.zno3, 1, 1, res.bsd, swi, res);
            end
            Cox_rate.Cox_sulfred = res.zTOC.calcReac(res.zno3, res.bsd.zinf, 1, 1, res.bsd, swi, res)
            
            % calculate mean OM concentration in upper x cm
            [C_10, C1_10, C2_10] = res.zTOC.calcC( 10, res.bsd, res.swi, res);
            OM_10=C_10* 100*12/res.bsd.rho_sed
            x = 10;
            Mean_OM = 1/x * 100*12/res.bsd.rho_sed*res.zTOC.calcOM(0.0, x, 1, 1, res.bsd, swi, res)
        end
        
        
        function res = test_benthic( ncl, swi )
            
            if nargin < 1
                ncl = 1;
            end
            
            res.bsd = benthic_main(ncl);
            res.bsd.usescalarcode = ncl==1;
            
            
            if nargin < 2 || isempty(swi)
                swi = benthic_test.default_swi();
            end
            
            if ncl > 1  % set up O2 gradient for testing
                O20 = swi.O20;
                for i = 1:ncl
                    swi.O20(i) = 10*(i-1)/(ncl-1)*O20;
                end
            end
            
            res.swi = swi;
            
            % calculate
            res.zTOC = benthic_zTOC(res.bsd);
            res.zO2 = benthic_zO2(res.bsd, res.swi);
            res.zNO3 = benthic_zNO3(res.bsd, res.swi);
            res.zSO4 = benthic_zSO4(res.bsd, res.swi);
            res.zNH4 = benthic_zNH4(res.bsd, res.swi);
            res.zH2S = benthic_zH2S(res.bsd, res.swi);
            res.zPO4_M = benthic_zPO4_M(res.bsd, res.swi);
            res.zDIC = benthic_zDIC(res.bsd, res.swi);
            res.zALK = benthic_zALK(res.bsd, res.swi);
            
            %            tic;
            res = res.zTOC.calc(res.bsd,res.swi, res);
            % check O2 demand using O2 to C ratio and (convert POC concentr. to flux analog to fortran)
            % POC_flux*OC = POC_conc * w * 1/(1 - por) * OC
            O2_demand_flux = -(res.swi.Fnonbio1+res.swi.Fnonbio2)*res.bsd.OC/((1-res.bsd.por)./res.bsd.por)
            %            O2_demand = (res.swi.C01+res.swi.C02)*res.bsd.OC
            %            O2_demand = (res.swi.C01+res.swi.C02)*res.bsd.w*res.bsd.OC
            if(res.swi.O20<=0.0)
                res.zox=0.0;
                res.flxzox = 0.0;
                res.conczox = 0.0;
                res.flxswiO2=0.0;
                res.zxf=0.0;
            else
                res = res.zO2.calc(res.bsd, res.swi, res);
            end
            if(swi.Nitrogen)
                res = res.zNO3.calc(res.bsd, res.swi, res);
            else
                res.zno3=res.zox;
                %                res.zso4=res.zox;   % for test-case with just TOC & O2
            end
            res = res.zSO4.calc(res.bsd, res.swi, res);
            if(swi.Nitrogen)
                res = res.zNH4.calc(res.bsd, res.swi, res);
            end
            res = res.zH2S.calc(res.bsd, res.swi, res);
            res = res.zPO4_M.calc(res.bsd, res.swi, res);
            res = res.zDIC.calc(res.bsd, res.swi, res);
            res = res.zALK.calc(res.bsd, res.swi, res);
            %            toc;
            
            %%%%% WRITE OUTPUT:
            answ = res
            [Cinf, C1inf, C2inf] = res.zTOC.calcC( 100, res.bsd, res.swi, res);
            [Cswi, C1swi, C2swi] = res.zTOC.calcC( 0, res.bsd, res.swi, res);
            fprintf('frac1 concentration at zinf %g \n',  C1inf);
            fprintf('frac2 concentration at zinf %g \n',  C2inf);
            fprintf('both concentration at zinf %g \n',  Cinf);
            fprintf('frac1 concentration at swi %g \n',  C1swi);
            fprintf('frac2 concentration at swi %g \n',  C2swi);
            fprintf('both concentration at swi %g \n',  Cswi);
            
            fprintf('sed preservation of POC %g \n',  Cinf/Cswi);
            %             %%% WRITE EXACT FLUX
            %             FO2_exact=res.zO2.calcFO2_exact(res.zox,res.bsd, res.swi, res);
            %             fprintf('exact F_O2 flux (mol cm^{-2} yr^{-1}) %g \n',  FO2_exact);
            
        end
        
        function plot_column(res, debug, swi, str_date)
            % plot single sediment column vs depth
            
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            bsd = res.bsd;
            zgrid = 0:0.1:bsd.zinf;
            
            
            % CONCENTRATIONS WITHOUT PO4
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            % Nitrogen included?
            if(swi.Nitrogen)
                figure;
                % TOC
                subplot(3,2,1)
                for i=1:length(zgrid)
                    [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                
                %%% TOC wt %
                plot(100*C1*12/bsd.rho_sed, -zgrid, 'b')
                hold on
                plot(100*C2*12/bsd.rho_sed, -zgrid, 'g')
                plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                
                %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                %                ylim([-50 0.0])
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
                %            title('Total TOC (wt%)')
                
                %%% O2
                if(res.zox>0.0)
                    for i=1:length(zgrid)
                        [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                    end
                else
                    O2(i) = 0.0;
                end
                subplot(3,2,3)
                plot(O2, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                ylim([-20 0.0])
                xlabel ('O_2 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('O2 (mol/cm^3)')
                
                %%% NO3
                for i=1:length(zgrid)
                    [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,2,5)
                plot(NO3, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                ylim([-20 0.0])
                xlabel ('NO_3 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('NO3 (mol/cm^3)')
                
                %%% NH4
                for i=1:length(zgrid)
                    [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,2,4)
                plot(NH4, -zgrid, 'b')
                hold on
%                xlim([0 1e-7])
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                hold off
                xlabel ('NH_4 (mol/cm^3)')
                %                ylabel('Depth (cm)')
                %            title ('NH4 (mol/cm^3)')
                
                %%% SO4
                subplot(3,2,2)
                for i=1:length(zgrid)
                    [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                end
                plot(SO4, -zgrid, 'b')
                hold on
                %                xlim([2.7e-5 swi.SO40])
                %                xlim([2.7e-5 swi.SO40])
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                hold off
                xlabel ('SO_4 (mol/cm^3)')
                %                ylabel('Depth (cm)')
                %            title ('SO4 (mol/cm^3)')
                
                %%% H2S
                subplot(3,2,6)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                end
                plot(H2S, -zgrid, 'b')
                hold on
                %                xlim([0 4e-7])
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('H_2S (mol/cm^3)')
                %               ylabel('Depth (cm)')
                %            title ('H2S (mol/cm^3)')
                
                % save Figure
                print('-depsc2', ['0_ALL_PROFILES_' str_date '.eps']);
                
            else % no Nitrogen
                
                figure;
                % TOC
                subplot(2,2,1)
                for i=1:length(zgrid)
                    [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                %%% TOC wt %
                plot(100*C1*12/bsd.rho_sed, -zgrid, 'b')
                hold on
                plot(100*C2*12/bsd.rho_sed, -zgrid, 'g')
                plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                
                %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
                %            title('Total TOC (wt%)')
                
                %%% O2
                for i=1:length(zgrid)
                    [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                end
                subplot(2,2,3)
                plot(O2, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('O_2 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('O2 (mol/cm^3)')
                
                %%% SO4
                subplot(2,2,2)
                for i=1:length(zgrid)
                    [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                end
                plot(SO4, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                hold off
                %            xlim([2.7e-5 swi.SO40])
                xlabel ('SO_4 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('SO4 (mol/cm^3)')
                
                %%% H2S
                subplot(2,2,4)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                end
                plot(H2S, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('H_2S (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('H2S (mol/cm^3)')
                
                
                print('-dpsc2', ['0_PROFILES_' str_date '.eps']);
            end
            
            % Also plot the 2nd figure for PO4, DIC and ALK?
            if(swi.plot_PO4_DIC_ALK)
                figure
                
                %%% PO4
                subplot(3,2,1)
                for i=1:length(zgrid)
                    [PO4(i), flxPO4(i), M(i), flxM(i), e_M(i), f_M(i), p_M(i), q_M(i), g_M(i), dedz_M(i), dfdz_M(i), dpdz_M(i), dqdz_M(i), dgdz_M(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
                end
                plot(PO4, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                %          axis([0 1.5*10^(-9) -100 0])
                xlabel ('PO_4 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('PO_4 (mol/cm^3)')
                
                %%% Fe-bound P (M)
                subplot(3,2,2)
                %for i=1:length(zgrid)
                %    [PO4(i), flxPO4(i), M(i), flxM(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
                %end
                plot(M, -zgrid, 'b')
                hold on
                %            plot([0,max(M)], [-bsd.zbio,-bsd.zbio], 'k--')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('Fe-bound P (mol/cm^3)')
                %            ylabel('Depth (cm)')
                %            title ('Fe-bound P (mol/cm^3)')
                
                %%% DIC
                subplot(3,2,3)
                for i=1:length(zgrid)
                    [DIC(i), flxDIC(i)] = res.zDIC.calcDIC(zgrid(i), bsd, res.swi, res);
                end
                plot(DIC, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('DIC (mol/cm^3)')
                ylabel('Depth (cm)')
                
                %%% ALK
                subplot(3,2,4)
                for i=1:length(zgrid)
                    [ALK(i), flxALK(i)] = res.zALK.calcALK(zgrid(i), bsd, res.swi, res);
                end
                plot(ALK, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('ALK (mol/cm^3)')
                ylabel('Depth (cm)')
                
                
                print('-depsc2', ['0_PO4_PROFILES_' str_date '.eps']);
            end
            
            
            %%% some plots for debugging!
            if debug
                figure
                
                subplot(3,2,1)
                hold on
                plot(e_M, -zgrid, 'b')
                title ('Fe-P - ODE solution')
                
                subplot(3,2,2)
                hold on
                plot(f_M, -zgrid, 'b')
                
                subplot(3,2,3)
                hold on
                plot(p_M, -zgrid, 'b')
                
                subplot(3,2,4)
                hold on
                plot(q_M, -zgrid, 'b')
                
                subplot(3,2,5)
                hold on
                plot(g_M, -zgrid, 'b')
                
                figure
                hold on
                subplot(3,2,1)
                plot(dedz_M, -zgrid, 'b')
                title ('Fe-P - ODE derivations')
                subplot(3,2,2)
                plot(dfdz_M, -zgrid, 'b')
                subplot(3,2,3)
                plot(dpdz_M, -zgrid, 'b')
                subplot(3,2,4)
                plot(dqdz_M, -zgrid, 'b')
                subplot(3,2,5)
                plot(dgdz_M, -zgrid, 'b')
                
                
                %%%%%%%%%%%%%%%%%%%%%
                
                %         H2S
                
                %%%%%%%%%%%%%%%%%%%%%
                figure
                subplot(3,3,2)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i), e_H2S(i), dedz_H2S(i), f_H2S(i), dfdz_H2S(i), g_H2S(i), dgdz_H2S(i)] = res.zH2S.calcH2S_debug(zgrid(i), bsd, res.swi, res);
                end
                plot(H2S, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('H_2S (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('H2S (mol/cm^3)')
                
                subplot(3,3,4)
                hold on
                plot(e_H2S, -zgrid, 'b')
                title ('H_2S - ODE solution')
                
                subplot(3,3,5)
                hold on
                plot(f_H2S, -zgrid, 'b')
                
                subplot(3,3,6)
                hold on
                plot(g_H2S, -zgrid, 'b')
                
                subplot(3,3,7)
                hold on
                plot(dedz_H2S, -zgrid, 'b')
                title ('H_2S - ODE derivations')
                
                subplot(3,3,8)
                hold on
                plot(dfdz_H2S, -zgrid, 'b')
                
                subplot(3,3,9)
                hold on
                plot(dgdz_H2S, -zgrid, 'b')
                
                
                
                % CONCENTRATION + Vertical Trransport
                figure;
                % TOC
                subplot(3,4,1)
                for i=1:length(zgrid)
                    [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                % TOC wt %
                plot(100*C1*12/bsd.rho_sed, -zgrid, 'b')
                hold on
                plot(100*C2*12/bsd.rho_sed, -zgrid, 'g')
                plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
                title('Total TOC (wt%)')
                % TOC vertical transport flux
                subplot(3,4,2);
                plot(C1flx, -zgrid, 'b')
                hold on
                plot(C2flx, -zgrid, 'g')
                plot(Cflx, -zgrid, 'k')
                xlabel ('TOC trspt (mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title('TOC vert transport')
                
                
                % O2
                for i=1:length(zgrid)
                    [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,4,3)
                plot(O2, -zgrid, 'b')
                hold on
                plot([0,res.swi.O20], [-bsd.zbio,-bsd.zbio], 'k--')
                xlabel ('O2 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('O2 (mol/cm^3)')
                subplot(3,4,4);
                plot(flxO2, -zgrid, 'b', flxO2D,-zgrid,'b--',flxO2adv,-zgrid,'c--');%,flxO2D+flxO2adv,-zgrid,'r--');
                legend('tot','diff','adv','diff+adv');
                legend boxoff;
                xlabel ('O2 trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('O2 vert transport')
                
                
                %             % NO3
                %
                %             for i=1:length(zgrid)
                %                 [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                %             end
                %             subplot(3,4,5)
                %             plot(NO3, -zgrid, 'b')
                %             hold on
                %             plot([0,res.swi.NO30], [-bsd.zbio,-bsd.zbio], 'k--')
                %             xlabel ('NO3 (mol/cm^3)')
                %             ylabel('Depth (cm)')
                %             title ('NO3 (mol/cm^3)')
                %             subplot(3,4,6)
                %             plot(flxNO3, -zgrid, 'b')
                %             xlabel ('NO3 trsp(mol cm^{-2}yr^{-1})')
                %             ylabel('Depth (cm)')
                %             title ('NO3 vert transport');
                %
                %
                %
                %             for i=1:length(zgrid)
                %                 [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                %             end
                %             subplot(3,4,7)
                %             plot(NH4, -zgrid, 'b')
                %             hold on
                %             plot([0,res.swi.NH40], [-bsd.zbio,-bsd.zbio], 'k--')
                %             hold off
                %             xlabel ('NH4 (mol/cm^3)')
                %             ylabel('Depth (cm)')
                %             title ('NH4 (mol/cm^3)')
                %             subplot(3,4,8)
                %             plot(flxNH4, -zgrid, 'b');
                %             xlabel ('NH4 trsp(mol cm^{-2}yr^{-1})')
                %             ylabel('Depth (cm)')
                %             title ('NH4 vert transport')
                
                subplot(3,4,9)
                for i=1:length(zgrid)
                    [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                end
                plot(SO4, -zgrid, 'b')
                hold on
                plot([0,res.swi.SO40], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                %xlim([0 SO40])
                xlabel ('SO4 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('SO4 (mol/cm^3)')
                subplot(3,4,10)
                plot(flxSO4, -zgrid, 'b');
                xlabel ('SO4 trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('SO4 vert transport')
                
                subplot(3,4,11)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                end
                plot(H2S, -zgrid, 'b')
                hold on
                plot([0,res.swi.H2S0], [-bsd.zbio,-bsd.zbio], 'k--')
                xlabel ('H2S (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('H2S (mol/cm^3)')
                subplot(3,4,12)
                plot(flxH2S, -zgrid, 'b');
                xlabel ('H2S trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('H2S vert transport')
                
            end
            
        end
        
        function test_k_SWIflux()
            clear
            swi=benthic_test.default_swi()
            %            % set date-time
            %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
            res=benthic_test.test_benthic(1,swi);
            benthic_test.plot_column(res, false, swi, 'k_0.1')
        end
        
        
        function run_OMEN_BRNS()
            % compare OMEN results with BRNS
            clear
            swi=benthic_test.default_swi()
            %             % set date-time
            %             str_date = datestr(now,'ddmmyy_HH_MM_SS');
            res=benthic_test.test_benthic(1,swi);
            benthic_test.plot_OMEN_BRNS(res, swi, '0107')
        end
        
        
        function res = OMEN_with_GENIE_input(bc, Nitrogen, k_parametr, k2scaling)
            
            swi.Nitrogen = Nitrogen;
            
            % Conversions after GENIE
            conv_POC_mol_cm3=12.0;
            conv_cal_mol_cm3=100.0/2.7;
            conv_det_mol_cm3=60.0/3;
            conv_cm3_kg = 1000;
            
            % set wdepth
            res.bsd.wdepth = -bc(end-2);
            
            res.bsd = benthic_main(1, res.bsd.wdepth);
            res.bsd.usescalarcode = true;
            
            %            if(Nitrogen)
            
            %            else
            
            % calculate sediment burial rate using w from GENIE from
            % previous time-step
            res.bsd.w = 1/(1-res.bsd.por)*bc(end);
            % calculate sediment accumulation rate using POC, CaCO3 and
            % OLD detrital burial/rain flux (convert from mol to cm3)
            %                res.bsd.w = 1/(1-res.bsd.por)*(conv_POC_mol_cm3*bc(end)+conv_cal_mol_cm3*bc(8) + conv_det_mol_cm3*bc(9));   % + bc(10))
            % without POC               res.bsd.w = 1/(1-res.bsd.por)*(conv_cal_mol_cm3*bc(8) + conv_det_mol_cm3*bc(9));   % + bc(10))
            if(bc(9)<=0.0)  % check detrital flux < 0.0
                % remineralise everything manually
                res.zox=0.0;
                res.zso4=0.0;
                res.Cox_rate_total=0.0; %res.zTOC.calcReac(0.0, res.bsd.zinf, 1, 1, res.bsd, res.swi, res);;
                res.Cox_rate_aerobic=0.0;
                res.Cox_rate_sulfred=0.0;
                res.Cox_perc_aerobic=0.0;
                res.Cox_perc_sulfred=0.0;
                res.C_zinf_wtpc=0.0;
                res.C1_zinf_wtpc=0.0;
                res.C2_zinf_wtpc=0.0;
                res.flxswiO2=0.0;
                res.flxswiSO4=0.0;
                res.flxswiH2S=0.0;
                res.flxswi_P=0.0;
                res.flxswiDIC=0.0;
                res.flxswiALK=0.0;
                res.Mean_OM=0.0;
                res.Cox_rate_total_xcm = 0.0;
                res.Cox_perc_aerobic_xcm = 0.0;
                
            else % USE OMEN-SED
                % check lower than detrital
                if(res.bsd.w<=bc(9))
                    res.bsd.w=bc(9);
                end
                if(res.bsd.w<=5.0e-4)
                    res.bsd.w=5.0e-4;
                end
                %bottom water concentrations
                swi.T = bc(11); %20.0;                         %temperature (degree C)
                
                % set fraction of labile and refractory POC
                f1=1-bc(end-1);
                f2=bc(end-1);
                swi.Fnonbio1 = f1*bc(1);    % [mol/(cm2 yr)] according non-bioturbated flux
                swi.Fnonbio2 = f2*bc(1);
                
                
                swi.C01 = 0.0;  % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
                swi.C02 = 0.0;
                
                % bottom water concentrations: GENIE mol/kg -> SEDIMENT MODEL needs mol/cm^3
                swi.O20= 1/conv_cm3_kg*bc(2);                   %O2  concentration at SWI (mol/cm^3)
                swi.SO40=1/conv_cm3_kg*bc(3);               	%SO4 concentration at SWI (mol/cm^3)
                swi.H2S0=1/conv_cm3_kg*bc(4);                 	%H2S concentration at SWI (mol/cm^3)
                swi.PO40=1/conv_cm3_kg*bc(5);                   %PO4 concentration at SWI (mol/cm^3)
                swi.Mflux0=365*0.2e-10;                         %flux of M to the sediment (mol/(cm2*yr))
                swi.DIC0=1/conv_cm3_kg*bc(6);               	%DIC concentration at SWI (mol/cm^3)
                swi.ALK0=1/conv_cm3_kg*bc(7);               	%ALK concentration at SWI (mol/cm^3)
                
                swi.plot_PO4_DIC_ALK=true;
                
                %            end
                % Set default values
                res.zTOC = benthic_zTOC(res.bsd);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%     Choose k parameterisation    %%%%%%%%%%
                
                % use oxic degradation rates
                switch k_parametr
                    case 'boudreau1997'
                        % use parameterisation of Boudreau 1997 dependent on sediment accumulation rate (w)
                        loc_k_apparent = 0.38*res.bsd.w^0.59;
                        if(res.bsd.wdepth < 1000)
                            k2scaling=5.0;
                        elseif(res.bsd.wdepth < 2000)
                            k2scaling=8.0;
                        elseif(res.bsd.wdepth < 3000)
                            k2scaling=12.0;
                            %elseif(res.bsd.wdepth < 4000)
                            %    k2scaling=50.0;
                            %elseif(res.bsd.wdepth < 5000)
                            %    k2scaling=100.0;
                        else
                            k2scaling=25.0;
                        end
                        
                        res.zTOC.k1=loc_k_apparent/(f1+f2*k2scaling);
                        res.zTOC.k2=res.zTOC.k1*k2scaling;
                        
                    case 'tromp1995'
                        OMEN_result(j,i) = test.zso4;
                    case 'stolpovsky2016'
                        OMEN_result(j,i) = test.Cox_rate_total;
                    case 'boudreau1997fPOC'
                        OMEN_result(j,i) = test.Cox_rate_total;
                    case 'invariant'
                        res.zTOC.k1=0.0065;
                        res.zTOC.k2=0.005;
                    otherwise
                        error('Error. Unknown k parameterization.')
                end
                
                % MIN oxic from Arndt et al. 2013
                %              res.zTOC.k1=1.0e-4;
                %              res.zTOC.k2=1.0e-6;
                
                %       after Boudreau 1997 - k dependent on OM flux (in micromol/(cm^2yr):
                %        res.zTOC.k1 = 2.2*1e-5*(bc(1)*10^6)^2.1;
                
                % if anoxic, decrease zbio and use anoxic degradation rate
                if(swi.O20 < 5.0e-9 )
                    res.bsd.zbio=0.01;
                    res.zTOC.k2 = res.zTOC.k1/650;
                    %                 switch k_parametr
                    %                     case 'boudreau1997'
                    %                         % use parameterisation of Boudreau 1997 dependent on sediment accumulation rate (w)
                    %                         % which is actually Toth and Lerman (1977) - as no anoxic rate in Boudreau 1997:
                    %                         loc_k_apparent = 0.04*res.bsd.w^2;
                    %                         res.zTOC.k1=loc_k_apparent/(f1+f2/100);
                    %                         res.zTOC.k2=res.zTOC.k1/100;
                    %                     case 'tromp1995'
                    %                         OMEN_result(j,i) = test.zso4;
                    %                     case 'stolpovsky2016'
                    %                         OMEN_result(j,i) = test.Cox_rate_total;
                    %                     case 'boudreau1997fPOC'
                    %                         OMEN_result(j,i) = test.Cox_rate_total;
                    %                     case 'invariant'
                    %                         OMEN_result(j,i) = test.Cox_rate_total;
                    %                    otherwise
                    %                     error('Error. Unknown k parameterization.')
                    %                 end
                end
                
                %            res.zTOC.k2 = res.zTOC.k1/100;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%        RUN OMEN       %%%%%%%%%%%%%%%%%
                
                
                res.swi = swi;
                
                res.zO2 = benthic_zO2(res.bsd, res.swi);
                res.zxf = 0.0;                              % roll off oxidation at low zox
                if(Nitrogen)
                    res.zNO3 = benthic_zNO3(res.bsd, res.swi);
                    res.zNH4 = benthic_zNH4(res.bsd, res.swi);
                end
                res.zSO4 = benthic_zSO4(res.bsd, res.swi);
                res.zH2S = benthic_zH2S(res.bsd, res.swi);
                res.zH2S = benthic_zH2S(res.bsd, res.swi);
                res.zPO4_M = benthic_zPO4_M(res.bsd, res.swi);
                res.zDIC = benthic_zDIC(res.bsd, res.swi);
                res.zALK = benthic_zALK(res.bsd, res.swi);
                
                
                % Calculate OMEN profiles and fluxes
                res = res.zTOC.calc(res.bsd,res.swi, res);
                if(swi.O20<= 1.0e-18)
                    res.zox = 0.0;
                    res.zno3=res.zox;
                    res.flxswiO2=0.0;
                else
                    res = res.zO2.calc(res.bsd, res.swi, res);
                end
                if(Nitrogen)
                    res = res.zNO3.calc(res.bsd, res.swi, res);
                else
                    res.zno3=res.zox;
                end
                
                res = res.zSO4.calc(res.bsd, res.swi, res);
                
                if(Nitrogen)
                    res = res.zNH4.calc(res.bsd, res.swi, res);
                end
                res = res.zH2S.calc(res.bsd, res.swi, res);
                res = res.zPO4_M.calc(res.bsd, res.swi, res);
                res = res.zDIC.calc(res.bsd, res.swi, res);
                res = res.zALK.calc(res.bsd, res.swi, res);
                
                %%%%% WRITE OUTPUT:
                sed_depth=100.0;
                answ = res;
                [Cinf, C1inf, C2inf] = res.zTOC.calcC( sed_depth, res.bsd, res.swi, res);
                [Cswi, C1swi, C2swi] = res.zTOC.calcC( 0, res.bsd, res.swi, res);
                
                %%%% TOC wt %%%%
                res.C1_zinf_wtpc=100*C1inf*12/res.bsd.rho_sed;
                res.C2_zinf_wtpc=100*C2inf*12/res.bsd.rho_sed;
                res.C_zinf_wtpc=100*Cinf*12/res.bsd.rho_sed;
                
                
                %             fprintf('frac1 concentration at zinf %g \n',  C1inf);
                %             fprintf('frac2 concentration at zinf %g \n',  C2inf);
                %             fprintf('both concentration at zinf %g \n',  Cinf);
                %             fprintf('frac1 concentration at swi %g \n',  C1swi);
                %             fprintf('frac2 concentration at swi %g \n',  C2swi);
                %             fprintf('both concentration at swi %g \n',  Cswi);
                %
                %             fprintf('sed preservation of POC %g \n',  Cinf/Cswi);
                
                x = 5;
                % calculate depth integrated OM degradation rates [mol cm-2 yr-1]
                res.Cox_rate_total_xcm = res.zTOC.calcReac(0.0, x, 1, 1, res.bsd, res.swi, res);
                res.Cox_rate_total = res.zTOC.calcReac(0.0, res.bsd.zinf, 1, 1, res.bsd, res.swi, res);
                res.Cox_rate_aerobic = res.zTOC.calcReac(0.0, res.zox, 1, 1, res.bsd, res.swi, res);
                if(swi.Nitrogen)
                    res.Cox_rate_denitr = res.zTOC.calcReac(res.zox, res.zno3, 1, 1, res.bsd, res.swi, res);
                end
                res.Cox_rate_sulfred = res.zTOC.calcReac(res.zno3, res.bsd.zinf, 1, 1, res.bsd, res.swi, res);
                
                res.Cox_perc_aerobic= res.Cox_rate_aerobic/res.Cox_rate_total*100;
                res.Cox_perc_sulfred= res.Cox_rate_sulfred/res.Cox_rate_total*100;
                
                if(res.zox>=x)
                    res.Cox_perc_aerobic_xcm = 100;
                else
                    res.Cox_perc_aerobic_xcm = res.Cox_rate_aerobic/res.Cox_rate_total_xcm *100;
                end
                
                % calculate mean OM concentration in upper x cm
                
                res.Mean_OM = 1/x * 100*12/res.bsd.rho_sed*res.zTOC.calcOM(0.0, x, 1, 1, res.bsd, res.swi, res);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%  TEST PROFILES  %%%%%%%%%%%%%%%%%%%%
                
                
                %   if(bsd.wdepth < 500)
                %               benthic_test.plot_column(res, false, res.swi, '_shallow_2804_GENIEw')
                %   end
            end
            
        end
        
        function swi = sensitivity_swi(swi, Params, str_date)
            % calculate sensitivity indices
            
            res.bsd = benthic_main(1);
            res.bsd.usescalarcode = true;
            
            if nargin < 2 || isempty(swi)
                swi = benthic_test.default_swi();
            end
            
            
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
            res.zDIC = benthic_zDIC(res.bsd, res.swi);
            res.zALK = benthic_zALK(res.bsd, res.swi);
            
            fileID = fopen(['./Sensitivity/Results_' str_date '.txt'],'w');
            fprintf(fileID,'%1s %8s %12s %12s %12s %12s %12s %12s %8s %8s\n','% Exp','F_O2','F_NO3','F_SO4','F_NH4','F_H2S','F_PO4','zox','zNO3','zSO4');
            fclose(fileID);
            
            for i=1:length(Params.k1)
                i
                res.swi.C01=Params.f1(i)*Params.wtpc(i)*1e-2/12*res.bsd.rho_sed;       %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                res.swi.C02=(1-Params.f1(i))*Params.wtpc(i)*1e-2/12*res.bsd.rho_sed;
                res.zTOC.k1 = Params.k1(i);
                res.zTOC.k2 = Params.k1(i)*0.01;
                %                 res.zNO3.KNH4 = Params.KNH4(i);
                %                 res.zPO4_M.KPO4_ox = Params.KPO4ox(i);
                %                 res.zPO4_M.KPO4_anox = Params.KPO4anox(i);
                %                 res.zPO4_M.ksPO4 = Params.ksPO4(i);
                %                 res.zPO4_M.kmPO4 = Params.kmPO4(i);
                %                 res.zPO4_M.kaPO4 = Params.kaPO4(i);
                %                 res.bsd.gamma = Params.gammaNH4(i);
                %                 res.bsd.gammaH2S = Params.gammaH2S(i);
                
                
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
                res = res.zDIC.calc(res.bsd, res.swi, res);
                res = res.zALK.calc(res.bsd, res.swi, res);
                
                swi.Results(i,:) = [i res.flxswiO2 res.flxswiNO3 res.flxswiSO4 res.flxswiNH4 res.flxswiH2S res.flxswi_P res.zox res.zno3 res.zso4];
                
                fileID = fopen(['./Sensitivity/Results_' str_date '.txt'],'a');
                fprintf(fileID,'%3d %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %8.5f %8.5f %8.5f\n',i, res.flxswiO2, res.flxswiNO3, res.flxswiSO4, res.flxswiNH4, res.flxswiH2S, res.flxswi_P, res.zox, res.zno3, res.zso4);
                fclose(fileID);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%  TEST PROFILES
                
                %           benthic_test.plot_column(res, false, swi, '0107')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             %%%%% WRITE OUTPUT:
                %             answ = res
                %             [Cinf, C1inf, C2inf] = res.zTOC.calcC( 100, res.bsd, res.swi, res);
                %             [Cswi, C1swi, C2swi] = res.zTOC.calcC( 0, res.bsd, res.swi, res);
                %             fprintf('frac1 concentration at zinf %g \n',  C1inf);
                %             fprintf('frac2 concentration at zinf %g \n',  C2inf);
                %             fprintf('both concentration at zinf %g \n',  Cinf);
                %             fprintf('frac1 concentration at swi %g \n',  C1swi);
                %             fprintf('frac2 concentration at swi %g \n',  C2swi);
                %             fprintf('both concentration at swi %g \n',  Cswi);
                %
                %             fprintf('sed preservation of POC %g \n',  Cinf/Cswi);
                %             %%% WRITE EXACT FLUX
                %             FO2_exact=res.zO2.calcFO2_exact(res.zox,res.bsd, res.swi, res);
                %             fprintf('exact F_O2 flux (mol cm^{-2} yr^{-1}) %g \n',  FO2_exact);
                %
                %
            end
            
        end
        
        function swi = sensitivity_swi_singleParameter(swi, Params, str_date)
            
            res.bsd = benthic_main(1);
            res.bsd.usescalarcode = true;
            
            if nargin < 2 || isempty(swi)
                swi = benthic_test.default_swi();
            end
            
            
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
            res.zDIC = benthic_zDIC(res.bsd, res.swi);
            res.zALK = benthic_zALK(res.bsd, res.swi);
            
            fileID = fopen(['./Sensitivity/Results_' str_date '.txt'],'w');
            %            fprintf(fileID,'%1s %8s %12s %12s %12s %12s %12s %12s %8s %8s\n','% Exp','F_O2','F_NO3','F_SO4','F_NH4','F_H2S','F_PO4','zox','zNO3','zSO4');
            fprintf(fileID,'%1s %8s %8s\n','% Exp','F_O2','zox');
            fclose(fileID);
            
            for i=1:length(Params.k1)
                i
                res.swi.C01=Params.f1*Params.wtpc(i)*1e-2/12*res.bsd.rho_sed;       %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                res.swi.C02=(1-Params.f1)*Params.wtpc(i)*1e-2/12*res.bsd.rho_sed;
                res.zTOC.k1 = Params.k1(i);
                res.zTOC.k2 = Params.k1(i)*0.01;
                %                 res.zNO3.KNH4 = Params.KNH4(i);
                %                 res.zPO4_M.KPO41 = Params.KPO4ox(i);
                %                 res.zPO4_M.KPO42 = Params.KPO4anox(i);
                %                 res.zPO4_M.ksPO4 = Params.ksPO4(i);
                %                 res.zPO4_M.kmPO4 = Params.kmPO4(i);
                %                 res.zPO4_M.kaPO4 = Params.kaPO4(i);
                %                 res.bsd.gamma = Params.gammaNH4(i);
                %                 res.bsd.gammaH2S = Params.gammaH2S(i);
                
                
                res = res.zTOC.calc(res.bsd,res.swi, res);
                res = res.zO2.calc(res.bsd, res.swi, res);
                % % % %
                % % % % %            if(swi.Nitrogen)
                % % % %             res = res.zNO3.calc(res.bsd, res.swi, res);
                % % % % %            else
                % % % % %            res.zno3=res.zox;
                % % % % %            end
                % % % %             res = res.zSO4.calc(res.bsd, res.swi, res);
                % % % % %            if(swi.Nitrogen)
                % % % %             res = res.zNH4.calc(res.bsd, res.swi, res);
                % % % % %            end
                % % % %             res = res.zH2S.calc(res.bsd, res.swi, res);
                % % % %             res = res.zPO4_M.calc(res.bsd, res.swi, res);
                % % % %             res = res.zDIC.calc(res.bsd, res.swi, res);
                % % % %             res = res.zALK.calc(res.bsd, res.swi, res);
                
                %           swi.Results(i,:) = [i res.flxswiO2 res.flxswiNO3 res.flxswiSO4 res.flxswiNH4 res.flxswiH2S res.flxswi_P res.zox res.zno3 res.zso4];
                swi.Results(i,:) = [i res.flxswiO2 res.zox];
                
                fileID = fopen(['./Sensitivity/Results_' str_date '.txt'],'a');
                %            fprintf(fileID,'%3d %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %8.5f %8.5f %8.5f\n',i, res.flxswiO2, res.flxswiNO3, res.flxswiSO4, res.flxswiNH4, res.flxswiH2S, res.flxswi_P, res.zox, res.zno3, res.zso4);
                fprintf(fileID,'%3d %7.6e %8.5f\n',i, res.flxswiO2, res.zox);
                fclose(fileID);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%  TEST PROFILES
                
                %           benthic_test.plot_column(res, false, swi, '0107')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             %%%%% WRITE OUTPUT:
                %             answ = res
                %             [Cinf, C1inf, C2inf] = res.zTOC.calcC( 100, res.bsd, res.swi, res);
                %             [Cswi, C1swi, C2swi] = res.zTOC.calcC( 0, res.bsd, res.swi, res);
                %             fprintf('frac1 concentration at zinf %g \n',  C1inf);
                %             fprintf('frac2 concentration at zinf %g \n',  C2inf);
                %             fprintf('both concentration at zinf %g \n',  Cinf);
                %             fprintf('frac1 concentration at swi %g \n',  C1swi);
                %             fprintf('frac2 concentration at swi %g \n',  C2swi);
                %             fprintf('both concentration at swi %g \n',  Cswi);
                %
                %             fprintf('sed preservation of POC %g \n',  Cinf/Cswi);
                %             %%% WRITE EXACT FLUX
                %             FO2_exact=res.zO2.calcFO2_exact(res.zox,res.bsd, res.swi, res);
                %             fprintf('exact F_O2 flux (mol cm^{-2} yr^{-1}) %g \n',  FO2_exact);
                %
                %
            end
            
        end
        
        function test_w()
            wdepth = 0:5000;
            w = benthic_main.sedrate(wdepth);
            
            figure;
            plot(w,-wdepth);
            xlabel('w, cm/yr');
            ylabel('depth, m');
        end
        
        function res = test_TOC_load( ncl, swi )
            clear
            swi=benthic_test.default_swi()
            a=500;
            
            for n = 1:a
                swi.C01=0.5*n/10*1e-2/12*2.5;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                swi.C02=0.5*n/10*1e-2/12*2.5;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                
                if nargin < 1
                    ncl = 1;
                end
                
                res.bsd = benthic_main(ncl);
                res.bsd.usescalarcode = ncl==1;
                
                if nargin < 2 || isempty(swi)
                    swi = benthic_test.default_swi();
                end
                
                if ncl > 1  % set up O2 gradient for testing
                    O20 = swi.O20;
                    for i = 1:ncl
                        swi.O20(i) = 10*(i-1)/(ncl-1)*O20;
                    end
                end
                
                res.swi = swi;
                
                
                % calculate
                res.zTOC = benthic_zTOC(res.bsd);
                res.zO2 = benthic_zO2(res.bsd, res.swi);
                res.zNO3 = benthic_zNO3(res.bsd, res.swi);
                res.zSO4 = benthic_zSO4(res.bsd, res.swi);
                res.zNH4 = benthic_zNH4(res.bsd, res.swi);
                res.zH2S = benthic_zH2S(res.bsd, res.swi);
                res.zH2S = benthic_zH2S(res.bsd, res.swi);
                res.zPO4_M = benthic_zPO4_M(res.bsd, res.swi);
                
                %           tic;
                res = res.zTOC.calc(res.bsd,res.swi, res);
                res = res.zO2.calc(res.bsd, res.swi, res);
                if(swi.Nitrogen)
                    res = res.zNO3.calc(res.bsd, res.swi, res);
                else
                    res.zno3=res.zox;
                end
                res = res.zSO4.calc(res.bsd, res.swi, res);
                if(swi.Nitrogen)
                    res = res.zNH4.calc(res.bsd, res.swi, res);
                end
                res = res.zH2S.calc(res.bsd, res.swi, res);
                res = res.zPO4_M.calc(res.bsd, res.swi, res);
                %            toc;
                
                res.O2depth(n)=res.zox;
                res.NO3depth(n)=res.zno3;
                res.SO4depth(n)=res.zso4;
            end
            
            %             %%%%% WRITE OUTPUT:
            %             answ = res
            %             [Cinf, C1inf, C2inf] = res.zTOC.calcC( 100, res.bsd, res.swi, res);
            %             [Cswi, C1swi, C2swi] = res.zTOC.calcC( 0, res.bsd, res.swi, res);
            %             fprintf('frac1 concentration at zinf %g \n',  C1inf);
            %             fprintf('frac2 concentration at zinf %g \n',  C2inf);
            %             fprintf('both concentration at zinf %g \n',  Cinf);
            %             fprintf('frac1 concentration at swi %g \n',  C1swi);
            %             fprintf('frac2 concentration at swi %g \n',  C2swi);
            %             fprintf('both concentration at swi %g \n',  Cswi);
            %
            %             fprintf('sed preservation of POC %g \n',  Cinf/Cswi);
            %             %%% WRITE EXACT FLUX
            %             FO2_exact=res.zO2.calcFO2_exact(res.zox,res.bsd, res.swi, res);
            %             fprintf('exact F_O2 flux (mol cm^{-2} yr^{-1}) %g \n',  FO2_exact);
            
            figure
            hold on;
            box on;
            plot((1:a)/10, -res.O2depth, 'b',(1:a)/10, -res.NO3depth, 'g', (1:a)/10, -res.SO4depth, 'r')
            xlabel('TOC (wt%)');
            ylabel('TEA penetration depth (cm)');
            hleg=legend('Oxygen', 'Nitrate', 'Sulfate');
            set(hleg,'Location','SouthWest')
            print('-depsc2', ['./Sensitivity/TOC_LOAD.eps']);
            
            
        end
        
        
        function res=test_budgets(gamma, zbio)
            % check conservation of fluxes at swi
            res.bsd = benthic_main(1);
            res.bsd.gamma = gamma;
            res.bsd.zbio = zbio;
            res.swi = benthic_test.default_swi();
            
            % fully oxic sediment: check swiO2 flux = net C + 2*swiNO3flx
            %            res.swi.O20 = 5000e-9;
            % calculate
            res.zTOC = benthic_zTOC(res.bsd);
            res.zO2 = benthic_zO2(res.bsd, res.swi);
            res.zNO3 = benthic_zNO3(res.bsd, res.swi);
            res = res.zTOC.calc(res.bsd,res.swi, res);
            res = res.zO2.calc(res.bsd, res.swi, res);
            res = res.zNO3.calc(res.bsd, res.swi, res);
            
            netC = res.Fswi_TOC - res.F_TOC;
            %           -*-ve    +ve NO3 from sed
            O2consump = -netC + 2*res.flxswiNO3;
            
            [swiO2, fswiO2, fDswiO2, fadvswiO2] = res.zO2.calcO2(0, res.bsd, res.swi, res);
            fprintf('swi O2 %g (mol cm^{-3}) flux (mol cm^{-2} yr^{-1}) tot %g D %g adv %g\n',  swiO2, fswiO2, fDswiO2, fadvswiO2);
            
            [zinfO2, fzinfO2, fDzinfO2, fadvzinfO2] = res.zO2.calcO2(res.bsd.zinf, res.bsd, res.swi, res);
            fprintf('zinf O2 %g (mol cm^{-3}) flux (mol cm^{-2} yr^{-1}) tot %g D %g adv %g\n',  zinfO2, fzinfO2, fDzinfO2, fadvzinfO2);
            
            netO2 = fswiO2 - fzinfO2;
            fprintf('netC %g net O2 %g swi NO3 %g (mol cm^{-2} yr^{-1}\n',netC, netO2, res.flxswiNO3);
            fprintf('O2 consumed %g  consumed - net input %g (mol cm^{-2} yr^{-1}\n', O2consump,O2consump + netO2);
            
        end
        
        function plot_OMEN_BRNS(res, swi, str_date)
            % plot single sediment column vs depth and compare with BRNS
            
            g1=load('OMEN-BRNS/FortranFiles_OLDADVECTION_4_1/g1.dat','ascii');
            g2=load('OMEN-BRNS/FortranFiles_OLDADVECTION_4_1/g2.dat','ascii');
            zzo2=load('OMEN-BRNS/FortranFiles_OLDADVECTION_4_1/zzo2.dat','ascii');
            zno3=load('OMEN-BRNS/FortranFiles_OLDADVECTION_4_1/zno3.dat','ascii');
            zso4=load('OMEN-BRNS/FortranFiles_OLDADVECTION_4_1/zso4.dat','ascii');
            zpo4=load('OMEN-BRNS/FortranFiles_OLDADVECTION_4_1/zpo4.dat','ascii');
            znh4=load('OMEN-BRNS/FortranFiles_OLDADVECTION_4_1/znh4.dat','ascii');
            zh2s=load('OMEN-BRNS/FortranFiles_OLDADVECTION_4_1/zh2s.dat','ascii');
            
            % Porosity and bioturbation coeff. depth invariant
            g1_di=load('OMEN-BRNS/FortranFilesDominik_4_1Test_depthinvariant/g1.dat','ascii');
            g2_di=load('OMEN-BRNS/FortranFilesDominik_4_1Test_depthinvariant/g2.dat','ascii');
            zzo2_di=load('OMEN-BRNS/FortranFilesDominik_4_1Test_depthinvariant/zzo2.dat','ascii');
            zno3_di=load('OMEN-BRNS/FortranFilesDominik_4_1Test_depthinvariant/zno3.dat','ascii');
            zso4_di=load('OMEN-BRNS/FortranFilesDominik_4_1Test_depthinvariant/zso4.dat','ascii');
            zpo4_di=load('OMEN-BRNS/FortranFilesDominik_4_1Test_depthinvariant/zpo4.dat','ascii');
            znh4_di=load('OMEN-BRNS/FortranFilesDominik_4_1Test_depthinvariant/znh4.dat','ascii');
            zh2s_di=load('OMEN-BRNS/FortranFilesDominik_4_1Test_depthinvariant/zh2s.dat','ascii');
            
            
            [l1, l2]=size(zzo2);
            l=122;  % depth steps
            %l=159;
            ylim=100; %zzo2(l,2); %39; %
            %            xlim = 10*10^(-10);
            mm=l1/(l)-1;
            
            set(0,'defaultLineLineWidth', 1)
            set(0,'DefaultAxesFontSize',12) % plots 18
            
            bsd = res.bsd;
            zgrid = 0:0.1:bsd.zinf;
            
            if(true)
                figure
                % PO4
                subplot(2,2,1)
                for i=1:length(zgrid)
                    [PO4(i), flxPO4(i), M(i), flxM(i), e_M(i), f_M(i), p_M(i), q_M(i), g_M(i), dedz_M(i), dfdz_M(i), dpdz_M(i), dqdz_M(i), dgdz_M(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
                end
                plot(PO4, -zgrid, 'b')
                hold on
                % BRNS
                plot(zpo4(mm*l+(1:l),1), -zpo4(mm*l+(1:l),2), 'r')
                plot(zpo4_di(mm*l+(1:l),1), -zpo4_di(mm*l+(1:l),2), 'g')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                axis ([-Inf Inf -ylim 0.0])
                %          axis([0 1.5*10^(-9) -100 0])
                xlabel ('PO_4 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('PO_4 (mol/cm^3)')
                
                % Fe-bound P (M)
                subplot(2,2,2)
                %for i=1:length(zgrid)
                %    [PO4(i), flxPO4(i), M(i), flxM(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
                %end
                plot(M, -zgrid, 'b')
                hold on
                %            plot([0,max(M)], [-bsd.zbio,-bsd.zbio], 'k--')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('Fe-bound P (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('Fe-bound P (mol/cm^3)')
                
                % DIC
                subplot(2,2,3)
                for i=1:length(zgrid)
                    [DIC(i), flxDIC(i)] = res.zDIC.calcDIC(zgrid(i), bsd, res.swi, res);
                end
                plot(DIC, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('DIC (mol/cm^3)')
                ylabel('Depth (cm)')
                
                % ALK
                subplot(2,2,4)
                for i=1:length(zgrid)
                    [ALK(i), flxALK(i)] = res.zALK.calcALK(zgrid(i), bsd, res.swi, res);
                end
                plot(ALK, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('ALK (mol/cm^3)')
                ylabel('Depth (cm)')
                
                
                %            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
                %            text(0.5, 1,'\bf Test 4\_1: 500m anoxic (no Mn, Fe)','HorizontalAlignment','center','VerticalAlignment', 'top')
                print('-dpsc2', ['./OMEN-BRNS/PO4_PROFILES_old_advection_' str_date '.ps']);
            end
            
            
            if(true)
                % CONCENTRATIONS WITHOUT PO4
                set(0,'defaultLineLineWidth', 1)
                set(0,'DefaultAxesFontSize',12)
                figure;
                % TOC
                for i=1:length(zgrid)
                    [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                % TOC wt %
                plot(100*C1*12/bsd.rho_sed, -zgrid, 'b-.')
                hold on
                plot(100*C2*12/bsd.rho_sed, -zgrid, 'b--')
                plot(100*C*12/bsd.rho_sed, -zgrid, 'b')
                % BRNS
                plot(g1(mm*l+(1:l),1)*12/2.5*100, -g1(mm*l+(1:l),2), 'r-.')
                plot(g2(mm*l+(1:l),1)*12/2.5*100, -g2(mm*l+(1:l),2), 'r--')
                plot(g1_di(mm*l+(1:l),1)*12/2.5*100, -g1_di(mm*l+(1:l),2), 'g-.')
                plot(g2_di(mm*l+(1:l),1)*12/2.5*100, -g2_di(mm*l+(1:l),2), 'g--')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                axis ([-Inf Inf -ylim 0.0])
                %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
                ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
                text(0.5, 1,'\bf OMEN (blue) vs BRNS depthinv. (green) vs old advection (red) - Test 4\_1: 500m anoxic (no Mn, Fe)','HorizontalAlignment','center','VerticalAlignment', 'top')
                print('-depsc2', ['./OMEN-BRNS/TOC_PROFILES_old_advection_' str_date '.eps']);
                
                if(swi.Nitrogen)
                    figure;
                    % TOC
                    subplot(3,2,1)
                    for i=1:length(zgrid)
                        [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                        [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                    end
                    % TOC wt %
                    plot(100*C1*12/bsd.rho_sed, -zgrid, 'b-.')
                    hold on
                    plot(100*C2*12/bsd.rho_sed, -zgrid, 'b--')
                    plot(100*C*12/bsd.rho_sed, -zgrid, 'b')
                    % BRNS
                    plot(g1(mm*l+(1:l),1)*12/2.5*100, -g1(mm*l+(1:l),2), 'r-.')
                    plot(g2(mm*l+(1:l),1)*12/2.5*100, -g2(mm*l+(1:l),2), 'r--')
                    plot(g1_di(mm*l+(1:l),1)*12/2.5*100, -g1_di(mm*l+(1:l),2), 'g-.')
                    plot(g2_di(mm*l+(1:l),1)*12/2.5*100, -g2_di(mm*l+(1:l),2), 'g--')
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    axis ([-Inf Inf -ylim 0.0])
                    %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                    hold off
                    xlabel ('TOC (wt%)')
                    ylabel('Depth (cm)')
                    %            title('Total TOC (wt%)')
                    
                    % O2
                    for i=1:length(zgrid)
                        [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                    end
                    subplot(3,2,3)
                    plot(O2, -zgrid, 'b')
                    hold on
                    % BRNS
                    plot(zzo2(mm*l+(1:l),1), -zzo2(mm*l+(1:l),2), 'r')
                    plot(zzo2_di(mm*l+(1:l),1), -zzo2_di(mm*l+(1:l),2), 'g')
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    axis ([-Inf Inf -5.0 0.0])
                    xlabel ('O_2 (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('O2 (mol/cm^3)')
                    
                    % NO3
                    
                    for i=1:length(zgrid)
                        [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                    end
                    subplot(3,2,5)
                    plot(NO3, -zgrid, 'b')
                    hold on
                    % BRNS
                    plot(zno3(mm*l+(1:l),1), -zno3(mm*l+(1:l),2), 'r')
                    plot(zno3_di(mm*l+(1:l),1), -zno3_di(mm*l+(1:l),2), 'g')
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    axis ([0.0 4e-8 -40.0 0.0])
                    xlabel ('NO_3 (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('NO3 (mol/cm^3)')
                    
                    
                    for i=1:length(zgrid)
                        [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                    end
                    subplot(3,2,4)
                    plot(NH4, -zgrid, 'b')
                    hold on
                    % BRNS
                    plot(znh4(mm*l+(1:l),1), -znh4(mm*l+(1:l),2), 'r')
                    plot(znh4_di(mm*l+(1:l),1), -znh4_di(mm*l+(1:l),2), 'g')
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    hold off
                    axis ([-Inf Inf -ylim 0.0])
                    xlabel ('NH_4 (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('NH4 (mol/cm^3)')
                    
                    subplot(3,2,2)
                    for i=1:length(zgrid)
                        [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                    end
                    plot(SO4, -zgrid, 'b')
                    hold on
                    % BRNS
                    plot(zso4(mm*l+(1:l),1), -zso4(mm*l+(1:l),2), 'r--')
                    plot(zso4_di(mm*l+(1:l),1), -zso4_di(mm*l+(1:l),2), 'g--')
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    hold off
                    axis ([-Inf Inf -ylim 0.0])
                    %xlim([0 SO40])
                    xlabel ('SO_4 (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('SO4 (mol/cm^3)')
                    
                    subplot(3,2,6)
                    for i=1:length(zgrid)
                        [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                    end
                    plot(H2S, -zgrid, 'b')
                    hold on
                    % BRNS
                    plot(zh2s(mm*l+(1:l),1), -zh2s(mm*l+(1:l),2), 'r--')
                    plot(zh2s_di(mm*l+(1:l),1), -zh2s_di(mm*l+(1:l),2), 'g--')
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    axis ([-Inf Inf -ylim 0.0])
                    xlabel ('H_2S (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('H2S (mol/cm^3)')
                    
                    %                ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
                    %                text(0.5, 1,'\bf OMEN (blue) vs BRNS (red);  depthinv. (green) - Test 4\_1: 500m anoxic (no Mn, Fe)','HorizontalAlignment','center','VerticalAlignment', 'top')
                    print('-dpsc2', ['./OMEN-BRNS/ALL_PROFILES_old_advection_' str_date '.ps']);
                    
                else
                    figure;
                    % TOC
                    subplot(2,2,1)
                    for i=1:length(zgrid)
                        [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                        [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                    end
                    % TOC wt %
                    plot(100*C1*12/bsd.rho_sed, -zgrid, 'b')
                    hold on
                    plot(100*C2*12/bsd.rho_sed, -zgrid, 'g')
                    plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    
                    %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                    hold off
                    xlabel ('TOC (wt%)')
                    ylabel('Depth (cm)')
                    %            title('Total TOC (wt%)')
                    
                    % O2
                    for i=1:length(zgrid)
                        [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                    end
                    subplot(2,2,3)
                    plot(O2, -zgrid, 'b')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    xlabel ('O_2 (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('O2 (mol/cm^3)')
                    
                    subplot(2,2,2)
                    for i=1:length(zgrid)
                        [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                    end
                    plot(SO4, -zgrid, 'b')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    hold off
                    %xlim([0 SO40])
                    xlabel ('SO_4 (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('SO4 (mol/cm^3)')
                    
                    subplot(2,2,4)
                    for i=1:length(zgrid)
                        [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                    end
                    plot(H2S, -zgrid, 'b')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    xlabel ('H_2S (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('H2S (mol/cm^3)')
                    
                    
                    print('-dpsc2', ['./OMEN-BRNS/ALL_PROFILES_' str_date '.ps']);
                end
                
            end
        end
        
        
        
        function plot_summary(res)
            % plot summary statistics from a set of runs
            
            figure;
            subplot(2,2,1);
            plot(res.swi.O20,-res.zox,res.swi.O20,-res.zno3,res.swi.O20,-res.zso4);
            xlabel('swi O2');
            ylabel('front depth (cm)');
            legend('zox','zno3','zso4');
            legend boxoff;
            
            subplot(2,2,2);
            plot(res.swi.O20,res.Fswi_TOC,'--');
            hold on;
            plot(res.swi.O20,res.Fswi_TOC-res.F_TOC);
            hold all;
            plot(res.swi.O20,res.flxswiO2);
            plot(res.swi.O20,res.flxswiNO3);
            plot(res.swi.O20,res.flxswiSO4);
            plot(res.swi.O20,res.flxswiNH4);
            plot(res.swi.O20,res.flxswiH2S);
            
            xlabel('swi O2 (\mu mol / L)');
            ylabel('flux (mol / cm^2 water column)');
            legend('C in','C net','O2 swi','NO3 swi','SO4 swi','NH4 swi','H2S swi');
            %legend boxoff;
        end
        
        
        function res = test_zO2(ncl)
            res.bsd = benthic_main(ncl);
            
            %bottom water concentrations
            swi.T=20.0;                                                     %temperature (degree C)
            swi.C01=0.01*1e-2/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm^3)
            swi.C02=0.01*1e-2/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm^3)
            
            O20=300.0e-9;                                                 %O2  concentration at SWI (mol/cm^3)
            
            if ncl == 1
                swi.O20 = O20;
            else
                for i=1:ncl
                    swi.O20(i) = O20*(i-1)/(ncl-1);
                end
            end
            
            npl = min(10,ncl);
            
            res.swi = swi;
            
            % calculate
            res.zTOC = benthic_zTOC(res.bsd);
            res = res.zTOC.calc(res.bsd,res.swi, res);
            
            res.zO2 = benthic_zO2(res.bsd, res.swi);
            
            figure;
            zgrid = 0:0.1:res.bsd.zinf;
            
            
            res.zox = res.bsd.zbio/2;
            
            [flxzox, conczox, flxswi,res] = res.zO2.calcbc(res.zox, res.bsd, res.swi, res, 1);
            for i=1:length(zgrid)
                [O2(i,:), flxO2(i,:)] = res.zO2.calcO2(zgrid(i), res.bsd, res.swi, res);
            end
            subplot(2,2,1)
            for i=1:npl
                plot(O2(:,i),-zgrid);
                hold all;
            end
            xlabel('O2');
            subplot(2,2,2);
            for i=1:npl
                plot(flxO2(:,i),-zgrid);
                hold all;
            end
            xlabel('flx O2');
            
            res.zox = res.bsd.zbio*2;
            
            [flxzox, conczox, flxswi, res] = res.zO2.calcbc(res.zox, res.bsd, res.swi, res, 1);
            for i=1:length(zgrid)
                [O2(i,:), flxO2(i,:)] = res.zO2.calcO2(zgrid(i), res.bsd, res.swi, res);
            end
            subplot(2,2,3)
            for i=1:npl
                plot(O2(:,i),-zgrid);
                hold all;
            end
            xlabel('O2');
            subplot(2,2,4);
            for i=1:npl
                plot(flxO2(:,i),-zgrid);
                hold all;
            end
            xlabel('flx O2');
            
            if res.bsd.ncl == 1
                figure;
                
                flxzox(1)=NaN;
                flxswi(1)=NaN;
                for i=2:length(zgrid)
                    [flxzoxz, conczox, flxswiz, res] = res.zO2.calcbc(zgrid(i), res.bsd, res.swi, res, 1);
                    flxzox(i)=flxzoxz;
                    flxswi(i)=flxswiz;
                    FO2(i) = res.zO2.calcFO2(zgrid(i), res.bsd, res.swi, res);
                end
                subplot(2,2,1)
                plot(flxzox,-zgrid,-FO2,-zgrid);
                subplot(2,2,2);
                plot(flxswi,-zgrid);
            end
            
        end
        
        function res = test_zNO3(ncl)
            res.bsd = benthic_main(ncl);
            
            %bottom water concentrations
            swi.T=20.0;                                                     %temperature (degree C)
            swi.C01=0.1/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm^3)
            swi.C02=0.1/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm^3)
            NO30=50.0e-9;                                               %NO3 concentration at SWI (mol/cm^3)
            
            if ncl == 1
                swi.NO30 = NO30;
            else
                for i=1:ncl
                    swi.NO30(i) = NO30*(i-1)/(ncl-1);
                end
            end
            
            npl = min(10,ncl);
            
            res.swi = swi;
            
            % calculate
            res.zTOC = benthic_zTOC(res.bsd);
            res = res.zTOC.calc(res.bsd,res.swi, res);
            
            res.zNO3 = benthic_zNO3(res.bsd, res.swi);
            
            figure;
            zgrid = 0:0.1:res.bsd.zinf;
            
            % invent some intf locations
            res.zox = res.bsd.zbio/2;
            res.zxf = 1;
            res.zno3 = res.zox*1.5;
            
            
            [flxzno3, conczno3, flxswi,res] = res.zNO3.calcbc(res.zno3, res.bsd, res.swi, res, 1);
            for i=1:length(zgrid)
                [NO3(i,:), flxNO3(i,:)] = res.zNO3.calcNO3(zgrid(i), res.bsd, res.swi, res);
            end
            subplot(2,2,1)
            for i=1:npl
                plot(NO3(:,i),-zgrid);
                hold all;
            end
            xlabel('NO3');
            subplot(2,2,2);
            for i=1:npl
                plot(flxNO3(:,i),-zgrid);
                hold all;
            end
            xlabel('flx NO3');
            
            res.zno3 = res.bsd.zbio*1.5;
            
            [flxzno3, conczno3, flxswi, res] = res.zNO3.calcbc(res.zno3, res.bsd, res.swi, res, 1);
            for i=1:length(zgrid)
                [NO3(i,:), flxNO3(i,:)] = res.zNO3.calcNO3(zgrid(i), res.bsd, res.swi, res);
            end
            subplot(2,2,3)
            for i=1:npl
                plot(NO3(:,i),-zgrid);
                hold all;
            end
            xlabel('NO3');
            subplot(2,2,4);
            for i=1:npl
                plot(flxNO3(:,i),-zgrid);
                hold all;
            end
            xlabel('flxNO3');
            
            if res.bsd.ncl == 1
                figure;
                
                flxzno3(1)=NaN;
                flxswi(1)=NaN;
                for i=2:length(zgrid)
                    [flxzno3z, conczno3, flxswiz, res] = res.zNO3.calcbc(zgrid(i), res.bsd, res.swi, res, 1);
                    flxzno3(i)=flxzno3z;
                    flxswi(i)=flxswiz;
                end
                subplot(2,2,1)
                plot(flxzno3,-zgrid);
                title('flx zno3');
                subplot(2,2,2);
                plot(flxswi,-zgrid);
                xlabel('flxswi');
            end
        end
        
        function res = test_zSO4(swi)
            res.bsd = benthic_main();
            
            if nargin < 1  % some defaults
                %bottom water concentrations
                swi.T=20.0;                                                     %temperature (degree C)
                swi.C01=0.1/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm^3)
                swi.C02=0.1/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm^3)
                swi.SO40=8000.0e-9;                                             %SO4 concentration at SWI (mol/cm^3)
            end
            
            res.swi = swi;
            
            % calculate
            res.zTOC = benthic_zTOC(res.bsd);
            res = res.zTOC.calc(res.bsd,res.swi, res);
            
            res.zSO4 = benthic_zSO4(res.bsd, res.swi);
            
            figure;
            zgrid = 0:0.1:res.bsd.zinf;
            
            % invent some intf locations
            res.zox = res.bsd.zbio/2;
            res.zxf = 1;
            %res.zno3 = res.zox*1.5;
            res.zno3 = res.zox;
            
            res.zso4 = res.bsd.zbio*2;
            
            [flxzso4, conczso4, flxswi,res] = res.zSO4.calcbc(res.zso4, res.bsd, res.swi, res, 1);
            for i=1:length(zgrid)
                [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), res.bsd, res.swi, res);
            end
            subplot(2,2,1)
            plot(SO4,-zgrid);
            xlabel('SO4');
            subplot(2,2,2);
            plot(flxSO4,-zgrid);
            xlabel('flx SO4');
            
            res.zso4 = res.bsd.zbio*4;
            [flxzso4, conczso4, flxswi, res] = res.zSO4.calcbc(res.zso4, res.bsd, res.swi, res, 1);
            for i=1:length(zgrid)
                [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), res.bsd, res.swi, res);
            end
            subplot(2,2,3)
            plot(SO4,-zgrid);
            xlabel('SO4');
            subplot(2,2,4);
            plot(flxSO4,-zgrid);
            xlabel('flxSO4');
            
            
            figure;
            
            flxzso4(1)=NaN;
            flxswi(1)=NaN;
            for i=2:length(zgrid)
                [flxzso4z, conczso4z, flxswiz, res] = res.zSO4.calcbc(zgrid(i), res.bsd, res.swi, res, 1);
                flxzso4(i)=flxzso4z;
                flxswi(i)=flxswiz;
                FSO4(i) = res.zSO4.calcFSO4(zgrid(i), res.bsd, res.swi, res);
            end
            subplot(2,2,1)
            plot(flxzso4,-zgrid,-FSO4,-zgrid);
            title('flx zso4, FSO4');
            subplot(2,2,2);
            plot(flxswi,-zgrid);
            xlabel('flxswi');
            
            % test zero flux at zso4 bc
            res.zso4 = res.bsd.zinf;
            
            [flxzso4,conczso4, flxswi,res] = res.zSO4.calcbc(res.zso4, res.bsd, res.swi, res, 2);
            for i=1:length(zgrid)
                [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), res.bsd, res.swi, res);
            end
            subplot(2,2,3)
            plot(SO4,-zgrid);
            xlabel('SO4');
            subplot(2,2,4);
            plot(flxSO4,-zgrid);
            xlabel('flxSO4');
            
        end
        
        
        function test_matchsoln
            E_l     = 1;    F_l     = 2;   G_l      = 3;
            dEdx_l  = 4;   dFdx_l   = 5;   dGdx_l   = 6;
            E_r     = 7;    F_r     = 8;   G_r      = 9;
            dEdx_r  = 10;  dFdx_r   = 11;  dGdx_r   = 12;
            Vb = 20;
            Db = 30;
            
            [a, b, c, d, e, f] = benthic_utils.matchsoln(E_l, F_l, G_l, dEdx_l, dFdx_l, dGdx_l, ...
                E_r, F_r, G_r, dEdx_r, dFdx_r, dGdx_r, ...
                Vb, Db);
            
            
            fprintf('a %g b %g c %g d %g e %g f %g\n',a,b,c,d,e,f);
            
            
            A_r = 1.2;
            B_r = 2.2;
            A_l = a*A_r + b*B_r + e;
            B_l = c*A_r + d*B_r + f;
            
            y_r     = A_r*E_r       + B_r*F_r       + G_r;
            dydx_r  = A_r*dEdx_r    + B_r*dFdx_r    + dGdx_r;
            
            y_l     = A_l*E_l       + B_l*F_l       + G_l;
            dydx_l  = A_l*dEdx_l    + B_l*dFdx_l    + dGdx_l;
            
            fprintf('r soln y %g  dydx %g\n', y_r, dydx_r);
            
            fprintf('l soln y %g  dydx %g\n', y_l, dydx_l);
            
            y_err       = y_l + Vb - y_r;
            dydx_err    = dydx_l + Db - dydx_r;
            
            fprintf('err   y %g  dydx %g\n', y_err, dydx_err);
            
        end
        
        function test_fzero_vec
            
            fzerooptions = optimset('Display','iter');
            
            for i =1:2
                fun = @(x) scalfun(x,i);
                xfz = fzero(fun,[-100 100], fzerooptions);
            end
            
            [xfzv,fval,cmplt,output]  = fzero_vec(@vecfun, -100*[1 1],100*[1 1],fzerooptions);
            
            fprintf('xfzv %g\n',xfzv);
            output
            
            function vf = vecfun(x)
                vf(1) = 1-x(1);
                vf(2) = 2-x(2)^3;
            end
            
            function f = scalfun(x,i)
                function test_k_SWIflux()
                    clear
                    swi=benthic_test.default_swi()
                    %            % set date-time
                    %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
                    res=benthic_test.test_benthic(1,swi);
                    benthic_test.plot_column(res, false, swi, 'k_0.1')
                end
                
                
                function run_OMEN_BRNS()
                    % compare OMEN results with BRNS
                    clear
                    swi=benthic_test.default_swi()
                    %             % set date-time
                    %             str_date = datestr(now,'ddmmyy_HH_MM_SS');
                    res=benthic_test.test_benthic(1,swi);
                    benthic_test.plot_OMEN_BRNS(res, swi, '0107')
                end
                
                xv =zeros(1,2);
                xv(i) = x;
                vf = vecfun(xv);
                f=vf(i);
            end
        end
        
        
    end
    
end

