classdef benthic_test
    % test cases for benthic layer model
    
    properties
    end
    
    methods(Static)
        function swi = default_swi()
            bsd = benthic_main();
            %bottom water concentrations
            swi.T=8.1; %20.0;                         %temperature (degree C)
            % see caption for Fig 1.2 - two equal TOC fractions 0.02 0.2 2
            swi.C01=0.05*1e-2/12*bsd.rho_sed; %1.45*0.5*400.0e-6*1/(10.0.^(-0.87478367-0.00043512*500.0)*3.3); % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm3 bulk phase)
            swi.C02=0.05*1e-2/12*bsd.rho_sed; %6.5*0.5*400.0e-6*1/(10.0.^(-0.87478367-0.00043512*500.0)*3.3); % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm3 bulk phase)
            %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm3 bulk phase)
            %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm3 bulk phase)
            swi.O20=300.e-9;   %was    300.0e-9                            %O2  concentration at SWI (mol/cm3)
            swi.NO30=20.0e-9;             % was 20.0e-9                                  %NO3 concentration at SWI (mol/cm3)
            swi.Nitrogen=true;
            swi.NH40=0.0e-9;                                                %NH4 concentration at SWI (mol/cm3)
            swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm3)
            %swi.SO40 = 100e-9;
            swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm3)
            swi.PO40=0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm3)
            swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996        

            swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm3)
            swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm3)
            swi.S0=35;                                                      %Salinity at SWI
        end
        
        function swi = sensitivity_swi()
            
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
            a=1000;
            for n = 1:a
            swi.C01=n*0.0005*1e-2*2.5;                                %TOC concentration at SWI (wt%) -> (mol/cm3 bulk phase)
            swi.C02=n*0.0005*1e-2*2.5;                                %TOC concentration at SWI (wt%) -> (mol/cm3 bulk phase)

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
            plot((1:a)/100, -res.O2depth, 'b',(1:a)/100, -res.NO3depth, 'g', (1:a)/100, -res.SO4depth, 'r')            
            xlabel('TOC (wt%)');
            ylabel('TEA penetration depth (cm)');
            hleg=legend('Oxygen', 'Nitrate', 'Sulfate');
            set(hleg,'Location','SouthWest')
            print('-dpsc2', ['TOC_LOAD.ps']);

            
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
            res.zH2S = benthic_zH2S(res.bsd, res.swi);
            res.zPO4_M = benthic_zPO4_M(res.bsd, res.swi);
            res.zDIC = benthic_zDIC(res.bsd, res.swi);
            res.zALK = benthic_zALK(res.bsd, res.swi);
   
            tic;
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
            res = res.zDIC.calc(res.bsd, res.swi, res);
            res = res.zALK.calc(res.bsd, res.swi, res);
            toc;
            
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
            %%% WRITE EXACT FLUX
            FO2_exact=res.zO2.calcFO2_exact(res.zox,res.bsd, res.swi, res);
            fprintf('exact F_O2 flux (mol cm^{-2} yr^{-1}) %g \n',  FO2_exact);
            
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
        
	function plot_OMEN_BRNS(res, swi)
            % plot single sediment column vs depth and compare with BRNS
            
            g1=load('../FortranFilesDominik_2Test/g1.dat','ascii');
            g2=load('../FortranFilesDominik_2Test/g2.dat','ascii');
            zzo2=load('../FortranFilesDominik_2Test/zzo2.dat','ascii');
            zno3=load('../FortranFilesDominik_2Test/zno3.dat','ascii');
            zso4=load('../FortranFilesDominik_2Test/zso4.dat','ascii');
            zpo4=load('../FortranFilesDominik_2Test/zpo4.dat','ascii');
            znh4=load('../FortranFilesDominik_2Test/znh4.dat','ascii');
            zh2s=load('../FortranFilesDominik_2Test/zh2s.dat','ascii');
   
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
            t=xlim;         % to draw penetration depths the correct lengths
            plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
            plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
            plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
            plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
            axis ([-Inf Inf -ylim 0.0])
  %          axis([0 1.5*10^(-9) -100 0])
            xlabel ('PO_4 (mol/cm3)')
            ylabel('Depth (cm)')
            title ('PO_4 (mol/cm3)')
           
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
            xlabel ('Fe-bound P (mol/cm3)')
            ylabel('Depth (cm)')
            title ('Fe-bound P (mol/cm3)')
            
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
            xlabel ('DIC (mol/cm3)')
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
            xlabel ('ALK (mol/cm3)')
            ylabel('Depth (cm)')
            

%                ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%                text(0.5, 1,'\bf Test 5: 5000m oxic (adjusted) - all PO_4 processes','HorizontalAlignment','center','VerticalAlignment', 'top')
            print('-dpsc2', ['PO4_PROFILES.ps']);
        end
           
            
	if(true)      
       % CONCENTRATIONS WITHOUT PO4
	set(0,'defaultLineLineWidth', 1)
	set(0,'DefaultAxesFontSize',12)

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
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
                axis ([-Inf Inf -5.0 0.0])
                xlabel ('O_2 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('O2 (mol/cm3)')

                % NO3

                for i=1:length(zgrid)
                    [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,2,5)
                plot(NO3, -zgrid, 'b')
                hold on
                % BRNS
                plot(zno3(mm*l+(1:l),1), -zno3(mm*l+(1:l),2), 'r')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')             
                axis ([-Inf Inf -5.0 0.0])
                xlabel ('NO_3 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('NO3 (mol/cm3)')


                for i=1:length(zgrid)
                    [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,2,4)
                plot(NH4, -zgrid, 'b')
                hold on
                % BRNS
                plot(znh4(mm*l+(1:l),1), -znh4(mm*l+(1:l),2), 'r')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
                hold off
                axis ([-Inf Inf -ylim 0.0])
                xlabel ('NH_4 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('NH4 (mol/cm3)')

                subplot(3,2,2)
                for i=1:length(zgrid)
                    [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                end
                plot(SO4, -zgrid, 'b')
                hold on
                % BRNS
                plot(zso4(mm*l+(1:l),1), -zso4(mm*l+(1:l),2), 'r--')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
                hold off
                axis ([-Inf Inf -ylim 0.0])
                %xlim([0 SO40])
                xlabel ('SO_4 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('SO4 (mol/cm3)')

                subplot(3,2,6)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                end
                plot(H2S, -zgrid, 'b')
                hold on
                % BRNS
                plot(zh2s(mm*l+(1:l),1), -zh2s(mm*l+(1:l),2), 'r--')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')          
                axis ([-Inf Inf -ylim 0.0])
                xlabel ('H_2S (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('H2S (mol/cm3)')

%                ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%                text(0.5, 1,'\bf OMEN (blue) vs BRNS (red) - Test 5: 5000m oxic (adjusted)','HorizontalAlignment','center','VerticalAlignment', 'top')
                print('-dpsc2', ['ALL_PROFILES.ps']);

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
                xlabel ('O_2 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('O2 (mol/cm3)')

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
                xlabel ('SO_4 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('SO4 (mol/cm3)')

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
                xlabel ('H_2S (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('H2S (mol/cm3)')


                print('-dpsc2', ['ALL_PROFILES.ps']);
            end

    end

  


        end
        
        function plot_column(res, debug, swi)
            % plot single sediment column vs depth
            
            set(0,'defaultLineLineWidth', 2)
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
            t=xlim;         % to draw penetration depths the correct lengths
            plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
            plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
            plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
            plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
  %          axis([0 1.5*10^(-9) -100 0])
            xlabel ('PO_4 (mol/cm3)')
            ylabel('Depth (cm)')
            title ('PO_4 (mol/cm3)')
           
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
            xlabel ('Fe-bound P (mol/cm3)')
            ylabel('Depth (cm)')
            title ('Fe-bound P (mol/cm3)')
            
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
            xlabel ('DIC (mol/cm3)')
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
            xlabel ('ALK (mol/cm3)')
            ylabel('Depth (cm)')
            

            print('-dpsc2', ['PO4_PROFILES.ps']);
        end
           
            
	if(true)      
       % CONCENTRATIONS WITHOUT PO4
	set(0,'defaultLineLineWidth', 2)
	set(0,'DefaultAxesFontSize',12)

             if(swi.Nitrogen)
                figure;
                % TOC
                subplot(3,2,1)
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
                subplot(3,2,3)
                plot(O2, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')               
                xlabel ('O_2 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('O2 (mol/cm3)')

                % NO3

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
                xlabel ('NO_3 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('NO3 (mol/cm3)')


                for i=1:length(zgrid)
                    [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,2,4)
                plot(NH4, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
                hold off
                xlabel ('NH_4 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('NH4 (mol/cm3)')

                subplot(3,2,2)
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
                xlabel ('SO_4 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('SO4 (mol/cm3)')

                subplot(3,2,6)
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
                xlabel ('H_2S (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('H2S (mol/cm3)')


                print('-dpsc2', ['ALL_PROFILES.ps']);

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
                xlabel ('O_2 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('O2 (mol/cm3)')

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
                xlabel ('SO_4 (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('SO4 (mol/cm3)')

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
                xlabel ('H_2S (mol/cm3)')
                ylabel('Depth (cm)')
    %            title ('H2S (mol/cm3)')


                print('-dpsc2', ['ALL_PROFILES.ps']);
            end

    end

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
            xlabel ('H_2S (mol/cm3)')
            ylabel('Depth (cm)')
%            title ('H2S (mol/cm3)')

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


        end
 
    
   if(false)      
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
            xlabel ('O2 (mol/cm3)')
            ylabel('Depth (cm)')
            title ('O2 (mol/cm3)')
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
%             xlabel ('NO3 (mol/cm3)')
%             ylabel('Depth (cm)')
%             title ('NO3 (mol/cm3)')
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
%             xlabel ('NH4 (mol/cm3)')
%             ylabel('Depth (cm)')
%             title ('NH4 (mol/cm3)')
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
            xlabel ('SO4 (mol/cm3)')
            ylabel('Depth (cm)')
            title ('SO4 (mol/cm3)')
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
            xlabel ('H2S (mol/cm3)')
            ylabel('Depth (cm)')
            title ('H2S (mol/cm3)')
            subplot(3,4,12)
            plot(flxH2S, -zgrid, 'b');          
            xlabel ('H2S trsp(mol cm^{-2}yr^{-1})')
            ylabel('Depth (cm)')
            title ('H2S vert transport')
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
            swi.C01=0.1/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm3)
            swi.C02=0.1/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm3)
            
            O20=4*900.0e-5*1e-3;                                                 %O2  concentration at SWI (mol/cm3)
            
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
            swi.C01=0.1/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm3)
            swi.C02=0.1/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm3)
            NO30=50.0e-9;                                               %NO3 concentration at SWI (mol/cm3)
            
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
                swi.C01=0.1/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm3)
                swi.C02=0.1/12*res.bsd.rho_sed;                                         %TOC concentration at SWI (mol/cm3)
                swi.SO40=8000.0e-9;                                             %SO4 concentration at SWI (mol/cm3)
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
                 xv =zeros(1,2);
                 xv(i) = x;
                 vf = vecfun(xv);
                 f=vf(i);
             end
         end
         
         
    end
    
end

