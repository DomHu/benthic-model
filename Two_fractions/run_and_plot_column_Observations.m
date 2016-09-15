         function run_and_plot_column_Observations()
            % run OMEN + plot single sediment column vs Observations
            
            % TODO: 
            % 1) give observations as argument -> change boundary
            % conditions and load appropriate observations!
             
             
            clear
            
            % was here: swi=benthic_test.default_swi()
            
            % initialise main model parameters with standard values
            bsd = benthic_main();
           
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            %sediment characteristics
            bsd.rho_sed=2.6; %was 2.5                           % sediment density (g/cm3)
            bsd.wdepth=585.0;     % Dom was 600.0                       % water depth (m)
            bsd.zbio=0.2;                              % bioturbation depth (cm)
            bsd.zinf=100;                               %Inifinity (cm)
            bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
            bsd.w = 10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper
            
            %bottom water concentrations
            swi.T=5.85; %20.0;                         %temperature (degree C)
            swi.C01= 2.0*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            swi.C02= 3.5*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            swi.O20=15.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
            swi.NO30=50.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
            swi.Nitrogen=true;
            swi.NH40=0.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
            swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
            %swi.SO40 = 100e-9;
            swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
            swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
            swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
            swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
            swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
            swi.S0=35;                                                      %Salinity at SWI


            %            % set date-time
            %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
            res=benthic_test.test_benthic(1,swi);


            str_date = 'Reimers_1409_nozbio_no3_80_Corg21';
            % LOAD Observations
            %Data Reimer et al. 1996
            data.TOC=load('../Observations/Reimers_SantaBarbara/BC21_Corg.dat','ascii');
%            TOC2=load('../Observations/Reimers_SantaBarbara/BC21_Corg.dat','ascii');
            data.O2=load('../Observations/Reimers_SantaBarbara/IMP_O2.dat','ascii');
            PW_data=load('../Observations/Reimers_SantaBarbara/BC_PW.dat','ascii');            
            data.NO3=PW_data(:,[1 9]);
            data.NH4=PW_data(:,[1 3]);
            data.SO4=PW_data(:,[1 10]);
            data.H2S=PW_data(:,[1 11]);
            data.PO4=PW_data(:,[1 2]);
% %             % Dale 
% %             data.TOC=load('../Observations/AndyDale/M92_Corg_250m_17MUC5_198MUC34.dat','ascii');
% %             data.SO4=load('../Observations/AndyDale/M92_SO4_250m_17MUC5_198MUC34.dat','ascii');
% %             data.H2S=load('../Observations/AndyDale/M92_H2S_250m_17MUC5_198MUC34.dat','ascii');
% %             data.NH4=load('../Observations/AndyDale/M92_NH4_250m_17MUC5_198MUC34.dat','ascii');
% %             data.PO4=load('../Observations/AndyDale/M92_PO4_250m_17MUC5.dat','ascii');
% %            % data.PO4=load('../Observations/AndyDale/M92_PO4_250m_17MUC5_198MUC34.dat','ascii');
    
            plot_column_Observations(res, str_date, data)
            

       
         end
         
         function plot_column_Observations(res, str_date, data)
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12) % plots 18
            
            bsd = res.bsd;
            zgrid = 0:0.1:bsd.zinf;
            
                figure;
                % TOC
                subplot(2,4,1)
                for i=1:length(zgrid)
                    [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                % TOC wt %
                scatter(data.TOC(:,2), -data.TOC(:,1),'k','filled')
                hold on
                plot(100*C1*12/bsd.rho_sed, -zgrid, 'r--')
                plot(100*C2*12/bsd.rho_sed, -zgrid, 'g--')
                plot(100*C*12/bsd.rho_sed, -zgrid, 'b')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')  
                box on;
    %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                ylim([-60.0 0.0])
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
    %            title('Total TOC (wt%)')
    
                % O2
                for i=1:length(zgrid)
                    [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                end
                subplot(2,4,2)
                scatter(data.O2(:,2).*1e-9, -data.O2(:,1),'k','filled')
                hold on
                plot(O2, -zgrid, 'b')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')               
                ylim([-5 0.0])
                box on;
                xlabel ('O_2 (mol/cm^3)')
                ylabel('Depth (cm)')
    %            title ('O2 (mol/cm^3)')

                % NO3

                for i=1:length(zgrid)
                    [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                end
                subplot(2,4,3)
                scatter(data.NO3(:,2).*1e-9, -data.NO3(:,1),'k','filled')
                hold on
                plot(NO3, -zgrid, 'b')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')             
                ylim([-10.0 0.0])
                xlabel ('NO_3 (mol/cm^3)')
                ylabel('Depth (cm)')
                box on;
    %            title ('NO3 (mol/cm^3)')


                for i=1:length(zgrid)
                    [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                end
                subplot(2,4,5)
                scatter(data.NH4(:,2).*1e-9, -data.NH4(:,1),'k','filled')
                hold on
                plot(NH4, -zgrid, 'b')
                xlim([0 2e-6])
                box on;
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
                hold off        
                ylim([-bsd.zinf 0.0])                
                xlabel ('NH_4 (mol/cm^3)')
%                ylabel('Depth (cm)')
    %            title ('NH4 (mol/cm^3)')

                subplot(2,4,4)
                for i=1:length(zgrid)
                    [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                end
                scatter(data.SO4(:,2).*1e-9, -data.SO4(:,1),'k','filled')
                hold on
                plot(SO4, -zgrid, 'b')
                box on;
%                xlim([2.7e-5 swi.SO40])     
          %      xlim([2.7e-5 swi.SO40])     
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
                hold off
                ylim([-bsd.zinf 0.0]) 
                xlabel ('SO_4 (mol/cm^3)')
%                ylabel('Depth (cm)')
    %            title ('SO4 (mol/cm^3)')

                subplot(2,4,6)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                end
                scatter(data.H2S(:,2).*1e-9, -data.H2S(:,1),'k','filled')
                hold on
                plot(H2S, -zgrid, 'b')
%                xlim([0 4e-7])
                box on;
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')  
                ylim([-bsd.zinf 0.0]) 
                xlabel ('H_2S (mol/cm^3)')
 %               ylabel('Depth (cm)')
    %            title ('H2S (mol/cm^3)')


            % PO4
            subplot(2,4,7)
            for i=1:length(zgrid)                
                [PO4(i), flxPO4(i), M(i), flxM(i), e_M(i), f_M(i), p_M(i), q_M(i), g_M(i), dedz_M(i), dfdz_M(i), dpdz_M(i), dqdz_M(i), dgdz_M(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
            end
         	scatter(data.PO4(:,2).*1e-9, -data.PO4(:,1),'k','filled')
            hold on            
            plot(PO4, -zgrid, 'b')
            box on;
            t=xlim;         % to draw penetration depths the correct lengths
            plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
            plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
            plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
            plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
  %          axis([0 1.5*10^(-9) -100 0])
            ylim([-bsd.zinf 0.0])
            xlabel ('PO_4 (mol/cm^3)')
            ylabel('Depth (cm)')
%            title ('PO_4 (mol/cm^3)')
           
            % Fe-bound P (M)
            subplot(2,4,8)
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
            ylim([-bsd.zinf 0.0])
            xlabel ('Fe-bound P (mol/cm^3)')
%            ylabel('Depth (cm)')
%            title ('Fe-bound P (mol/cm^3)')
    
           print('-depsc2', ['eps_output/0_ALL_PROFILES_' str_date '.eps']);

           
            figure
            % PO4
            subplot(1,2,1)
            for i=1:length(zgrid)                
                [PO4(i), flxPO4(i), M(i), flxM(i), e_M(i), f_M(i), p_M(i), q_M(i), g_M(i), dedz_M(i), dfdz_M(i), dpdz_M(i), dqdz_M(i), dgdz_M(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
            end
            plot(PO4, -zgrid, 'b')
            hold on            
         	scatter(data.PO4(:,2).*1e-9, -data.PO4(:,1),'k','filled')
            t=xlim;         % to draw penetration depths the correct lengths
            plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
            plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
            plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
            plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
  %          axis([0 1.5*10^(-9) -100 0])
            ylim([-bsd.zinf 0.0])
            xlabel ('PO_4 (mol/cm^3)')
            ylabel('Depth (cm)')
%            title ('PO_4 (mol/cm^3)')
           
            % Fe-bound P (M)
            subplot(1,2,2)
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
            ylim([-bsd.zinf 0.0])
            xlabel ('Fe-bound P (mol/cm^3)')
%            ylabel('Depth (cm)')
%            title ('Fe-bound P (mol/cm^3)')
            

 %           print('-depsc2', ['eps_output/0_PO4_PROFILES_' str_date '.eps']);
          
            
                
            % CONCENTRATIONS WITHOUT PO4

                figure;
                % TOC
%                subplot(3,2,1)
                for i=1:length(zgrid)
                    [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                % TOC wt %
                scatter(data.TOC(:,2), -data.TOC(:,1),'k','filled')
                hold on
                plot(100*C1*12/bsd.rho_sed, -zgrid, 'r--')
                plot(100*C2*12/bsd.rho_sed, -zgrid, 'g--')
                plot(100*C*12/bsd.rho_sed, -zgrid, 'b')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')  
                box on;
    %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                ylim([-60.0 0.0])
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
    %            title('Total TOC (wt%)')
    
           

         end

         
% % % %  Observations Dale        
% % % %     if(false)
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 
% % % %             
% % % %             %%%%    Station 1: Middle shelf water depth 74m
% % % %             
% % % %             %sediment characteristics
% % % %             bsd.rho_sed=2.0; %was 2.5                           % sediment density (g/cm3)
% % % %             bsd.wdepth=74.0;     % Dom was 600.0                       % water depth (m)
% % % %             bsd.zbio=4.0;                              % bioturbation depth (cm)
% % % %             bsd.zinf=100;                               %Inifinity (cm)
% % % %             bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %10.0;
% % % %             bsd.w = 0.45;
% % % % 
% % % %             %bottom water concentrations
% % % %             swi.T=14.0; %20.0;                         %temperature (degree C)
% % % %             swi.C01= 4*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             swi.C02= 4*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             swi.O20=0.0;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
% % % %             swi.NO30=0.1e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
% % % %             swi.Nitrogen=true;
% % % %             swi.NH40=2.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
% % % %             swi.SO40=29000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
% % % %             %swi.SO40 = 100e-9;
% % % %             swi.H2S0=30.0e-9;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
% % % %             swi.PO40=40.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
% % % %             swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
% % % %             swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
% % % %             swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
% % % %             swi.S0=35;                                                      %Salinity at SWI
% % % % 
% % % % 
% % % %             %            % set date-time
% % % %             %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
% % % %             res=benthic_test.test_benthic(1,swi);
% % % % 
% % % % 
% % % %             str_date = 'MiddleShelf_1808_zbio4_Dbiowdepth_k2004';
% % % %             % LOAD Observations
% % % %             data.TOC=load('../Observations/AndyDale/M92_Corg_74m_220MUC39.dat','ascii');
% % % %             data.SO4=load('../Observations/AndyDale/M92_SO4_74m_220MUC39.dat','ascii');
% % % %             data.H2S=load('../Observations/AndyDale/M92_H2S_74m_220MUC39.dat','ascii');
% % % %             data.NH4=load('../Observations/AndyDale/M92_NH4_74m_220MUC39.dat','ascii');
% % % %             data.PO4=load('../Observations/AndyDale/M92_PO4_74m_220MUC39.dat','ascii');
% % % % 
% % % %             plot_column_Observations(res, str_date, data)
% % % %             
% % % % 
% % % %             %%%%    Station 5: Outer shelf water depth 195m
% % % %             
% % % %             %sediment characteristics
% % % %             bsd.rho_sed=2.0; %was 2.5                           % sediment density (g/cm3)
% % % %             bsd.wdepth=195.0;     % Dom was 600.0                       % water depth (m)
% % % %             bsd.zbio=4.0;                              % bioturbation depth (cm)
% % % %             bsd.zinf=50;                               %Inifinity (cm)
% % % %             bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %1.0;
% % % %             bsd.w = 0.1;
% % % % 
% % % %             %bottom water concentrations
% % % %             swi.T=13.0; %20.0;                         %temperature (degree C)
% % % %             swi.C01= 3*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             swi.C02= 11.5*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             swi.O20=0.0;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
% % % %             swi.NO30=7.8e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
% % % %             swi.Nitrogen=true;
% % % %             swi.NH40=2.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
% % % %             swi.SO40=29000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
% % % %             %swi.SO40 = 100e-9;
% % % %             swi.H2S0=0.0e-9;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
% % % %             swi.PO40=40.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
% % % %             swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
% % % %             swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
% % % %             swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
% % % %             swi.S0=35;                                                      %Salinity at SWI
% % % % 
% % % % 
% % % %             %            % set date-time
% % % %             %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
% % % %             res=benthic_test.test_benthic(1,swi);
% % % % 
% % % % 
% % % %             str_date = 'OuterShelf_1808_zbio4_050008';
% % % %             % LOAD Observations
% % % %             data.TOC=load('../Observations/AndyDale/M92_Corg_195m_247MUC45.dat','ascii');
% % % %             data.SO4=load('../Observations/AndyDale/M92_SO4_195m_247MUC45.dat','ascii');
% % % %             data.H2S=load('../Observations/AndyDale/M92_H2S_195m_247MUC45.dat','ascii');
% % % %             data.NH4=load('../Observations/AndyDale/M92_NH4_195m_247MUC45.dat','ascii');
% % % %             data.PO4=load('../Observations/AndyDale/M92_PO4_195m_247MUC45.dat','ascii');
% % % % 
% % % %             plot_column_Observations(res, str_date, data)         
% % % % 
% % % %             
% % % %             %%%%    Station 10: Below OMZ water depth 1015m
% % % %             
% % % %             %sediment characteristics
% % % %             bsd.rho_sed=2.0; %was 2.5                           % sediment density (g/cm3)
% % % %             bsd.wdepth=1015.0;     % Dom was 600.0                       % water depth (m)
% % % %             bsd.zbio=2.0;                              % bioturbation depth (cm)
% % % %             bsd.zinf=30;                               %Inifinity (cm)
% % % %             bsd.Dbio=0.01;
% % % %             bsd.w = 0.06;
% % % %             bsd.por=0.76;
% % % %             
% % % %             %bottom water concentrations
% % % %             swi.T=4.4; %20.0;                         %temperature (degree C)
% % % %             swi.C01= 2*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             swi.C02= 2*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
% % % %             swi.O20=50.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
% % % %             swi.NO30=47.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
% % % %             swi.Nitrogen=true;
% % % %             swi.NH40=0.697e-9;                                                %NH4 concentration at SWI (mol/cm^3)
% % % %             swi.SO40=29000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
% % % %             %swi.SO40 = 100e-9;
% % % %             swi.H2S0=0.0e-9;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
% % % %             swi.PO40=40.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
% % % %             swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
% % % %             swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
% % % %             swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
% % % %             swi.S0=35;                                                      %Salinity at SWI
% % % % 
% % % % 
% % % %             %            % set date-time
% % % %             %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
% % % %             res=benthic_test.test_benthic(1,swi);
% % % % 
% % % % 
% % % %             str_date = 'BelowOMZ_1708';
% % % %             % LOAD Observations
% % % %             data.TOC=load('../Observations/AndyDale/M92_Corg_1015m_163MUC30.dat','ascii');
% % % %             data.SO4=load('../Observations/AndyDale/M92_SO4_1015m_163MUC30.dat','ascii');
% % % %             data.H2S=load('../Observations/AndyDale/M92_H2S_1015m_163MUC30.dat','ascii');
% % % %             data.NH4=load('../Observations/AndyDale/M92_NH4_1015m_163MUC30.dat','ascii');
% % % %             data.PO4=load('../Observations/AndyDale/M92_PO4_1015m_163MUC30.dat','ascii');
% % % % 
% % % %             plot_column_Observations(res, str_date, data)
% % % %             
% % % %             
% % % % 
% % % %             
% % % %     end            
% % % %             
% % % %          
