         function run_and_plot_column_Observations()
            % run OMEN + plot single sediment column vs Observations
            
            % TODO: 
            % 1) give observations as argument -> change boundary
            % conditions and load appropriate observations!
             
             
            clear
%         gamma=0.95   except 4298m: 0.9                             %fraction of NH4 that is oxidised in oxic layer
%         gammaH2S=0.95;                           %fraction of H2S that is oxidised in oxic layer
%         gammaCH4=0.99;                           %fraction of CH4 that is oxidised at SO4
           
            
            % was here: swi=benthic_test.default_swi()
%% Step 1:  initialise main model parameters with standard values & run model
            bsd = benthic_main();
            Obs = 10;       
            %sediment characteristics
            switch Obs
                    case 1  % OMEXDIA_2809_108m all solutes in Micromoles/litre
% %         k1= 0.65;                                                %TOC degradation rate constnat (1/yr)
% %         k2=0.00001;  
% %         KPO41=200.0;  %Adsorption coefficient in oxic layer (-) 
% %         KPO42=1.3;    %Adsorption coefficient in anoxic layer (-)
% %         ksPO4=1.0;    %Rate constant for kinetic P sorption (1/yr)   0.12 fits 1.CASE; 2.2 fits 2. CASE DOM: was 0.5*365 from Nicolas; Slomp ea 1996 0.26
% %         kmPO4=2.2e-6*24*365;          %Rate constant for Fe-bound P release upon Fe oxide reduction   DOM: was 1.8e-6 Slomp ea 1996 0.00053*365 
% %         kaPO4=10.0; %Rate constant for authigenic P formation (1/yr)    DOM: was 0.004*365 from Nicolas; Slomp ea 1996 0.001
% %         PO4s=1.0e-9;        %Equilibrium concentration for P sorption (mol/cm3)       was 1.5e-9; ; Slomp ea 1996
% %         PO4a= 1.5e-8; % %Equilibrium concentration for authigenic P formation (mol/cm3) was 0.7e-9
% %         Minf=1.0e-10;       % asymptotic concentration for Fe-bound P (mol/cm3)      TODO/CHECK: good value? is from Slomp et al. 1996 Dom was 1.99e-6

                        bsd.rho_sed=2.6; %was 2.5                           % sediment density (g/cm3)
                        bsd.wdepth=108.0;     % Dom was 600.0                       % water depth (m)
                        bsd.zbio=1.0;                              % bioturbation depth (cm)
                        bsd.zinf=50;                               %Inifinity (cm)
                        bsd.Dbio=0.02; %5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        bsd.w = 10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=12.5; %20.0;                         %temperature (degree C)
                        swi.C01= 2.64*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.C02= 1.8*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.O20=210.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=9.6e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.4e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
                        
                    case 2 % OMEXDIA_2809_2213m
%        k1= 0.15;                                                %TOC degradation rate constnat (1/yr)
%        k2=0.0009;     
%         ksPO4=1.0; %0.26*365;      %Rate constant for kinetic P sorption (1/yr)   0.12 fits 1.CASE; 2.2 fits 2. CASE DOM: was 0.5*365 from Nicolas; Slomp ea 1996 0.26
%        % ksPO4=1e-15;
%         %kmPO4= 1e-15 ;
%         kmPO4=2.2e-6*24*365;          % Dom was from Slomp 0.00053*365;	%Rate constant for Fe-bound P release upon Fe oxide reduction   DOM: was 1.8e-6 Slomp ea 1996 0.00053*365 
%         %kaPO4 = 0.0;
%         kaPO4=10.0; % Dom was 0.001*365;	%Rate constant for authigenic P formation (1/yr)    DOM: was 0.004*365 from Nicolas; Slomp ea 1996 0.001
%         PO4s=1.0e-9;        %Equilibrium concentration for P sorption (mol/cm3)       was 1.5e-9; ; Slomp ea 1996
%         PO4a= 0.5e-8; %47e-9;  %was 3.7e-9      %Equilibrium concentration for authigenic P formation (mol/cm3) was 0.7e-9

                        bsd.rho_sed=2.6; %was 2.5                           % sediment density (g/cm3)
                        bsd.wdepth=2213.0;     % Dom was 600.0                       % water depth (m)
                        bsd.zbio=10.0;                              % bioturbation depth (cm)
                        bsd.zinf=50;                               %Inifinity (cm)
                        bsd.Dbio=0.17; %5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        bsd.w = 10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=3.2; %20.0;                         %temperature (degree C)
                        swi.C01= 0.45*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.C02= 0.6*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.O20=250.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=25.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.6e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
           
                    case 3  % OMEXDIA_2809_3097m
%        k1= 0.13;                                                %TOC degradation rate constnat (1/yr)
%        k2=0.0008;     
                        bsd.rho_sed=2.6; %was 2.5                           % sediment density (g/cm3)
                        bsd.wdepth=3097.0;     % Dom was 600.0                       % water depth (m)
                        bsd.zbio=10.0;                              % bioturbation depth (cm)
                        bsd.zinf=50;                               %Inifinity (cm)
                        bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        bsd.w = 10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T= 2.6; %20.0;                         %temperature (degree C)
                        swi.C01= 2.0*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.C02= 1.8*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.O20=243.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=21.3e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.55e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
                       
                    case 4  % OMEXDIA_3371m
%        k1= 0.13;                                                %TOC degradation rate constnat (1/yr)
%        k2=0.0008;     
                        bsd.rho_sed=2.6; %was 2.5                           % sediment density (g/cm3)
                        bsd.wdepth=3371.0;     % Dom was 600.0                       % water depth (m)
                        bsd.zbio=5.0;                              % bioturbation depth (cm)
                        bsd.zinf=50;                               %Inifinity (cm)
                        bsd.Dbio=0.11; %5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        bsd.w = 10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T= 3.1; %20.0;                         %temperature (degree C)
                        swi.C01= 0.55*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.C02= 0.6*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.O20=245.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=22.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.27e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=1.4e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
                        
                    case 5  % OMEXDIA_4941m
%        k1= 0.13;                                                %TOC degradation rate constnat (1/yr)
%        k2=0.0008;     
                        bsd.rho_sed=2.6; %was 2.5                           % sediment density (g/cm3)
                        bsd.wdepth=4941.0;     % Dom was 600.0                       % water depth (m)
                        bsd.zbio=10.0;                              % bioturbation depth (cm)
                        bsd.zinf=50;                               %Inifinity (cm)
                        bsd.Dbio=0.44; %5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        bsd.w = 10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T= 2.5; %20.0;                         %temperature (degree C)
                        swi.C01= 0.25*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.C02= 0.3*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.O20=230.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=28.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=4.5e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=1.4e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
                        
                    case 6  % Reimers et al. 1996
% %         k1= 0.1;                                                %TOC degradation rate constnat (1/yr)
% %         k2=0.001;  
% %         KPO41=200.0;  %Adsorption coefficient in oxic layer (-) 
% %         KPO42=1.3;    %Adsorption coefficient in anoxic layer (-)
% %         ksPO4=1.0;    %Rate constant for kinetic P sorption (1/yr)   0.12 fits 1.CASE; 2.2 fits 2. CASE DOM: was 0.5*365 from Nicolas; Slomp ea 1996 0.26
% %         kmPO4=2.2e-6*24*365;          %Rate constant for Fe-bound P release upon Fe oxide reduction   DOM: was 1.8e-6 Slomp ea 1996 0.00053*365 
% %         kaPO4=10.0; %Rate constant for authigenic P formation (1/yr)    DOM: was 0.004*365 from Nicolas; Slomp ea 1996 0.001
% %         PO4s=1.0e-9;        %Equilibrium concentration for P sorption (mol/cm3)       was 1.5e-9; ; Slomp ea 1996
% %         PO4a= 9.0e-8; % %Equilibrium concentration for authigenic P formation (mol/cm3) was 0.7e-9
% %         Minf=1.0e-10;       % asymptotic concentration for Fe-bound P (mol/cm3)      TODO/CHECK: good value? is from Slomp et al. 1996 Dom was 1.99e-6

                        bsd.rho_sed=2.6; %was 2.5                           % sediment density (g/cm3)
                        bsd.wdepth=585.0;     % Dom was 600.0                       % water depth (m)
                        bsd.zbio=0.01;                              % bioturbation depth (cm)
                        bsd.zinf=50;                               %Inifinity (cm)
                        bsd.Dbio=0.02; %5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        bsd.w = 0.42; % 10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=5.85; %20.0;                         %temperature (degree C)
                        swi.C01= 2.0*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.C02= 3.5*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.O20=10.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=25e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=50.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
                        
                    case 7  % OMEXDIA_4908m - Obs not so nice
%        k1= 0.13;                                                %TOC degradation rate constnat (1/yr)
%        k2=0.0008;     
                        bsd.rho_sed=2.6; %was 2.5                           % sediment density (g/cm3)
                        bsd.wdepth=4908.0;     % Dom was 600.0                       % water depth (m)
                        bsd.zbio=10.0;                              % bioturbation depth (cm)
                        bsd.zinf=50;                               %Inifinity (cm)
                        bsd.Dbio=0.15; %5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        bsd.w = 10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T= 3.1; %20.0;                         %temperature (degree C)
                        swi.C01= 0.72*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.C02= 0.2*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.O20=245.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=22.7e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.38e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
                    
                    case 8  % OMEXDIA_2809_343m all solutes in Micromoles/litre
%        k1= 0.15;                                                %TOC degradation rate constnat (1/yr)
%        k2=0.0009;     
                        bsd.rho_sed=2.6; %was 2.5                           % sediment density (g/cm3)
                        bsd.wdepth=343.0;     % Dom was 600.0                       % water depth (m)
                        bsd.zbio=5.0;                              % bioturbation depth (cm)
                        bsd.zinf=50;                               %Inifinity (cm)
                        bsd.Dbio=0.08; % 5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        bsd.w = 10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=2.3; %20.0;                         %temperature (degree C)
                        swi.C01= 1.7*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.C02= 1.9*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.O20=206.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=11.4e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.22e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
                        
                      case 9 % Dale
                          
                      case 10  % OMEXDIA_4298m does not look nice
%        k1= 0.055;                                                %TOC degradation rate constnat (1/yr)
%        k2=0.00001;  
%         ksPO4=1.0; %0.26*365;      %Rate constant for kinetic P sorption (1/yr)   0.12 fits 1.CASE; 2.2 fits 2. CASE DOM: was 0.5*365 from Nicolas; Slomp ea 1996 0.26
%        % ksPO4=1e-15;
%         %kmPO4= 1e-15 ;
%         kmPO4=2.2e-6*24*365;          % Dom was from Slomp 0.00053*365;	%Rate constant for Fe-bound P release upon Fe oxide reduction   DOM: was 1.8e-6 Slomp ea 1996 0.00053*365 
%         %kaPO4 = 0.0;
%         kaPO4=10.0; % Dom was 0.001*365;	%Rate constant for authigenic P formation (1/yr)    DOM: was 0.004*365 from Nicolas; Slomp ea 1996 0.001
%         PO4s=1.0e-9;        %Equilibrium concentration for P sorption (mol/cm3)       was 1.5e-9; ; Slomp ea 1996
%         PO4a= 0.5e-8; %47e-9;  %was 3.7e-9      %Equilibrium concentration for authigenic P formation (mol/cm3) was 0.7e-9
% gamma=0.9;                                %fraction of NH4 that is oxidised in oxic layer
                        bsd.rho_sed=2.6; %was 2.5                           % sediment density (g/cm3)
                        bsd.wdepth=4298.0;     % Dom was 600.0                       % water depth (m)
                        bsd.zbio=4.2;                              % bioturbation depth (cm)
                        bsd.zinf=50;                               %Inifinity (cm)
                        bsd.Dbio=0.18; %5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        bsd.w = 10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T= 2.5; %20.0;                         %temperature (degree C)
                        swi.C01= 1.0*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.C02= 1.2*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C01=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02=0.0005*1e-2*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        swi.O20=243.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=30.1e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.22e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2000.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
       
                    otherwise
                        error('unrecognized Obs  %g\n',Obs);
                end

            res=benthic_test.test_benthic(1,swi);

 %% Step 2: LOAD Observations   

           switch Obs
                    case 1  % OMEXDIA_2809_108m all solutes in Micromoles/litre
                        str_date = '108m_OMEXDIA_2110_';
                        data.TOC=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','Corg','C2:D24');     % in wt%
                        data.O2=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','O2','C2:D126'); 
                        data.NO3=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','NO3','C2:D25'); 
                        data.NH4=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','NH4','C2:D25');
                        data.PO4=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','PO4','C2:D25');
            %            data.SO4=PW_data(:,[1 10]);
            %            data.H2S=PW_data(:,[1 11]);
                        
                    case 2 % OMEXDIA_2809_2213m
                        str_date = '2213m_OMEXDIA_2110_';
                        data.TOC=xlsread('../Observations/OMEXDIA/2_PE121_98-4_2213m.xlsx','Corg','C2:D39');     % in wt%
                        data.O2=xlsread('../Observations/OMEXDIA/2_PE121_98-4_2213m.xlsx','O2','C113:D1352'); % C2:D112
                        data.NO3=xlsread('../Observations/OMEXDIA/2_PE121_98-4_2213m.xlsx','NO3','C2:D18'); 
                        data.NH4=xlsread('../Observations/OMEXDIA/2_PE121_98-4_2213m.xlsx','NH4','C2:D18');
                        data.PO4=xlsread('../Observations/OMEXDIA/2_PE121_98-4_2213m.xlsx','PO4','C2:D18');
            %            data.SO4=PW_data(:,[1 10]);
            %            data.H2S=PW_data(:,[1 11]);
            
                    case 3  % OMEXDIA_2809_3097m
                         str_date = '3097m_OMEXDIA_2110_';
                        data.TOC=xlsread('../Observations/OMEXDIA/3_PE138_99-14.xlsx','Corg','C2:D31');     % in wt%
                        data.O2=xlsread('../Observations/OMEXDIA/3_PE138_99-14.xlsx','O2','C2:D391'); % C2:D112
                        data.NO3=xlsread('../Observations/OMEXDIA/3_PE138_99-14.xlsx','NO3','C2:D26'); 
                        data.NH4=xlsread('../Observations/OMEXDIA/3_PE138_99-14.xlsx','NH4','C2:D26');
                        data.PO4=xlsread('../Observations/OMEXDIA/3_PE138_99-14.xlsx','PO4','C2:D25');
            %            data.SO4=PW_data(:,[1 10]);
            %            data.H2S=PW_data(:,[1 11]);
            
                     case 4  % OMEXDIA_3371m
                        str_date = '3371m_OMEXDIA_2110_';
                        data.TOC=xlsread('../Observations/OMEXDIA/9_PE121_98-02_3371m.xlsx','Corg','C2:D29');     % in wt%
                        data.O2=xlsread('../Observations/OMEXDIA/9_PE121_98-02_3371m.xlsx','O2','C2:D130'); % C2:D112
                        data.NO3=xlsread('../Observations/OMEXDIA/9_PE121_98-02_3371m.xlsx','NO3','C2:D18'); 
                        data.NH4=xlsread('../Observations/OMEXDIA/9_PE121_98-02_3371m.xlsx','NH4','C2:D18');
                        data.PO4=xlsread('../Observations/OMEXDIA/9_PE121_98-02_3371m.xlsx','PO4','C2:D18');
            %            data.SO4=PW_data(:,[1 10]);
            %            data.H2S=PW_data(:,[1 11]);                        
                      
                    case 5  % OMEXDIA_4941m
                        str_date = '4941m_OMEXDIA_2110_';
                        data.TOC=xlsread('../Observations/OMEXDIA/8_PE121_98-01_4941m.xlsx','Corg','C2:D21');     % in wt%
                        data.O2=xlsread('../Observations/OMEXDIA/8_PE121_98-01_4941m.xlsx','O2','C2:D87'); % C2:D112
                        data.NO3=xlsread('../Observations/OMEXDIA/8_PE121_98-01_4941m.xlsx','NO3','C2:D24'); 
                        data.NH4=xlsread('../Observations/OMEXDIA/8_PE121_98-01_4941m.xlsx','NH4','C2:D24');
                        data.PO4=xlsread('../Observations/OMEXDIA/8_PE121_98-01_4941m.xlsx','PO4','C2:D24');
            %            data.SO4=PW_data(:,[1 10]);
            %            data.H2S=PW_data(:,[1 11]);                        
                        
                    case 6  % Reimers et al. 1996   NOT really clear which measurements fit together -> leave out!
                        str_date = '585m_Reimers_BC68_2110_';
                        data.TOC=load('../Observations/Reimers_SantaBarbara/BC68_Corg.dat','ascii');
            %            TOC2=load('../Observations/Reimers_SantaBarbara/BC21_Corg.dat','ascii');
                        data.O2=load('../Observations/Reimers_SantaBarbara/IMP_O2.dat','ascii');
                        data.O2(:,1)=data.O2(:,1)/10;
                        PW_data=load('../Observations/Reimers_SantaBarbara/BC_PW.dat','ascii');   
                        data.NO3=PW_data(:,[1 9]);
                        data.NH4=load('../Observations/Reimers_SantaBarbara/BC68_NH4.dat','ascii');
                        %data.NH4=PW_data(:,[1 3]);
                        data.SO4=PW_data(:,[1 10]);
                        data.H2S=load('../Observations/Reimers_SantaBarbara/BC68_H2S.dat','ascii');
                        %data.H2S=PW_data(:,[1 11]);
                        data.PO4=load('../Observations/Reimers_SantaBarbara/BC68_PO4.dat','ascii');
                       % data.PO4=PW_data(:,[1 2]);
                        
                    case 7  % 
                        str_date = '4908m_OMEXDIA_2110_';
                        data.TOC=xlsread('../Observations/OMEXDIA/7_PE121_98-06_4908m.xlsx','Corg','C2:D38');     % in wt%
                        data.O2=xlsread('../Observations/OMEXDIA/7_PE121_98-06_4908m.xlsx','O2','C2:D87'); % C2:D112
                        data.NO3=xlsread('../Observations/OMEXDIA/7_PE121_98-06_4908m.xlsx','NO3','C2:D20'); 
                        data.NH4=xlsread('../Observations/OMEXDIA/7_PE121_98-06_4908m.xlsx','NH4','C2:D20');
                        data.PO4=xlsread('../Observations/OMEXDIA/7_PE121_98-06_4908m.xlsx','PO4','C2:D20');
            %            data.SO4=PW_data(:,[1 10]);
            %            data.H2S=PW_data(:,[1 11]);

                       
                    case 8  % OMEXDIA_2809_343m all solutes in Micromoles/litre
                        str_date = '343m_OMEXDIA_2110_';
                        data.TOC=xlsread('../Observations/OMEXDIA/1_PE138_99-12.xlsx','Corg','C2:D25');     % in wt%
                        data.O2=xlsread('../Observations/OMEXDIA/1_PE138_99-12.xlsx','O2','C2:D249'); 
                        data.NO3=xlsread('../Observations/OMEXDIA/1_PE138_99-12.xlsx','NO3','C2:D27'); 
                        data.NH4=xlsread('../Observations/OMEXDIA/1_PE138_99-12.xlsx','NH4','C2:D27');
                        data.PO4=xlsread('../Observations/OMEXDIA/1_PE138_99-12.xlsx','PO4','C2:D26');
            %            data.SO4=PW_data(:,[1 10]);
            %            data.H2S=PW_data(:,[1 11]);
                        
                    case 9 % Dale 
                        data.TOC=load('../Observations/AndyDale/M92_Corg_250m_17MUC5_198MUC34.dat','ascii');
                        data.SO4=load('../Observations/AndyDale/M92_SO4_250m_17MUC5_198MUC34.dat','ascii');
                        data.H2S=load('../Observations/AndyDale/M92_H2S_250m_17MUC5_198MUC34.dat','ascii');
                        data.NH4=load('../Observations/AndyDale/M92_NH4_250m_17MUC5_198MUC34.dat','ascii');
                        data.PO4=load('../Observations/AndyDale/M92_PO4_250m_17MUC5.dat','ascii');
                       % data.PO4=load('../Observations/AndyDale/M92_PO4_250m_17MUC5_198MUC34.dat','ascii');                   

                    case 10  % OMEXDIA_4298m
                        str_date = '4298m_OMEXDIA_2110_';
                        data.TOC=xlsread('../Observations/OMEXDIA/4_PE138_99-17_4298m.xlsx','Corg','C2:D38');     % in wt%
                        data.O2=xlsread('../Observations/OMEXDIA/4_PE138_99-17_4298m.xlsx','O2','C2:D815'); % C2:D112
                        data.NO3=xlsread('../Observations/OMEXDIA/4_PE138_99-17_4298m.xlsx','NO3','C2:D20'); 
                        data.NH4=xlsread('../Observations/OMEXDIA/4_PE138_99-17_4298m.xlsx','NH4','C2:D20');
                        data.PO4=xlsread('../Observations/OMEXDIA/4_PE138_99-17_4298m.xlsx','PO4','C2:D20');
            %            data.SO4=PW_data(:,[1 10]);
            %            data.H2S=PW_data(:,[1 11]);                        
  
                    otherwise
                        error('unrecognized Obs  %g\n',Obs);
                end


%% Step 3: PLOT MODEL + Observations           

            plot_column_Observations(res, str_date, data, Obs)
            

       
         end
         
         function plot_column_Observations(res, str_date, data, Obs)
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
                ylim([-bsd.zinf 0.0])
                xlabel ('TOC (wt%)')
  %              ylabel('Depth (cm)')
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
                ylim([-10 0.0])
                box on;
                xlabel ('O_2 (mol/cm^3)')
    %            ylabel('Depth (cm)')
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
                ylim([-50.0 0.0])
                xlabel ('NO_3 (mol/cm^3)')
    %            ylabel('Depth (cm)')
                box on;
    %            title ('NO3 (mol/cm^3)')


                for i=1:length(zgrid)
                    [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                end
                subplot(2,4,5)
                scatter(data.NH4(:,2).*1e-9, -data.NH4(:,1),'k','filled')
                hold on
                plot(NH4, -zgrid, 'b')
                xlim([0.0 100e-9])     
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
                if(Obs == 6 || Obs == 9)
                    scatter(data.SO4(:,2).*1e-9, -data.SO4(:,1),'k','filled')
                end
                hold on
                plot(SO4, -zgrid, 'b')
                xlim([2.0e-5 3.0e-5])     
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
                if(Obs == 6 || Obs == 9)
                    scatter(data.H2S(:,2).*1e-9, -data.H2S(:,1),'k','filled')
                end
                hold on
                plot(H2S, -zgrid, 'b')
                xlim([0 10e-9])
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
            xlim([0.0 10e-9])
            t=xlim;         % to draw penetration depths the correct lengths
            plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
            plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
            plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
            plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
  %          axis([0 1.5*10^(-9) -100 0])
            ylim([-bsd.zinf 0.0])
            xlabel ('PO_4 (mol/cm^3)')
 %           ylabel('Depth (cm)')
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
    
           print('-depsc2', ['eps_output/' str_date '_PROFILES.eps']);

% %            
% %             figure
% %             % PO4
% %             subplot(1,2,1)
% %             for i=1:length(zgrid)                
% %                 [PO4(i), flxPO4(i), M(i), flxM(i), e_M(i), f_M(i), p_M(i), q_M(i), g_M(i), dedz_M(i), dfdz_M(i), dpdz_M(i), dqdz_M(i), dgdz_M(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
% %             end
% %             plot(PO4, -zgrid, 'b')
% %             hold on            
% %          	scatter(data.PO4(:,2).*1e-9, -data.PO4(:,1),'k','filled')
% %             t=xlim;         % to draw penetration depths the correct lengths
% %             plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
% %             plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
% %             plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
% %             plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
% %   %          axis([0 1.5*10^(-9) -100 0])
% %             ylim([-bsd.zinf 0.0])
% %             xlabel ('PO_4 (mol/cm^3)')
% %             ylabel('Depth (cm)')
% % %            title ('PO_4 (mol/cm^3)')
% %            
% %             % Fe-bound P (M)
% %             subplot(1,2,2)
% %             %for i=1:length(zgrid)                
% %             %    [PO4(i), flxPO4(i), M(i), flxM(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
% %             %end
% %             plot(M, -zgrid, 'b')
% %             hold on
% % %            plot([0,max(M)], [-bsd.zbio,-bsd.zbio], 'k--')        
% %             t=xlim;         % to draw penetration depths the correct lengths
% %             plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
% %             plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
% %             plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
% %             plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')  
% %             ylim([-bsd.zinf 0.0])
% %             xlabel ('Fe-bound P (mol/cm^3)')
% % %            ylabel('Depth (cm)')
% % %            title ('Fe-bound P (mol/cm^3)')
% %             
% % 
% % %           print('-depsc2', ['eps_output/0_PO4_PROFILES_' str_date '.eps']);
% %           
            
                
            % CONCENTRATIONS WITHOUT PO4
    if(false)
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
