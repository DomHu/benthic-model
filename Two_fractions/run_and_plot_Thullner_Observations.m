         function run_and_plot_Thullner_Observations(pObs, pplot)
            % run OMEN + plot global hypsometry
            
 % --> pObs in [1, 7]
 % --> pplot true, false
             
%            clear
%         gamma=0.95   except 4298m: 0.9                             %fraction of NH4 that is oxidised in oxic layer
%         gammaH2S=0.95;                           %fraction of H2S that is oxidised in oxic layer
%         gammaCH4=0.99;                           %fraction of CH4 that is oxidised at SO4
           
            swi=benthic_test.default_swi();            

%            Obs = 1;       
            %sediment characteristics

            switch pObs
                    case 1  % 100m
                        res.bsd.wdepth=100.0;     % Dom was 600.0                       % water depth (m)
                        res.bsd = benthic_main(1, res.bsd.wdepth);

                        res.zTOC = benthic_zTOC(res.bsd);
                        res.zTOC.k1= 0.00001*0.221;                                                %TOC degradation rate constnat (1/yr)
                        res.zTOC.k2=0.221;  

                        res.bsd.rho_sed=2.5; %was 2.5                           % sediment density (g/cm3)
                        res.bsd.zbio=10.0;                              % bioturbation depth (cm)
                        res.bsd.zinf=100;                               %Inifinity (cm)
                        res.bsd.por=0.85; 
                        %bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        %bsd.w = %10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=10.3; %20.0;                         %temperature (degree C)
                        %swi.C01_nonbio= 2.64*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02_nonbio= 1.8*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                     	swi.FPOM = 510e-6;      % Total flux after Thullner et al. 2009 [mol/(cm2 yr)]
                        swi.Fnonbio1 = 0.0000001*swi.FPOM;    % [mol/(cm2 yr)] according non-bioturbated flux
                        swi.Fnonbio2 = swi.FPOM;
                        swi.C01 = 0.0;  % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
                        swi.C02 = 0.0;
                        swi.O20=132.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=17.3e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        % for Pacific
%                         swi.O20=0.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
%                         swi.NO30=40.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2400.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
                        
                    case 2 %  200m

                        res.bsd.wdepth=200.0;     % Dom was 600.0                       % water depth (m)
                        res.bsd = benthic_main(1, res.bsd.wdepth);
                        res.zTOC = benthic_zTOC(res.bsd);
                       	res.zTOC.k1= 0.00001*0.2085;                                                %TOC degradation rate constnat (1/yr)
                       	res.zTOC.k2=0.2085;  

                        res.bsd.rho_sed=2.5; %was 2.5                           % sediment density (g/cm3)
                        res.bsd.zbio=10.0;                              % bioturbation depth (cm)
                        res.bsd.zinf=100;                               %Inifinity (cm)
                        res.bsd.por=0.85; 
                        %bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        %bsd.w = %10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=9.7; %20.0;                         %temperature (degree C)
                        %swi.C01_nonbio= 2.64*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02_nonbio= 1.8*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                     	swi.FPOM = 467e-6;      % Total flux after Thullner et al. 2009 [mol/(cm2 yr)]
                        swi.Fnonbio1 = 0.00001*swi.FPOM;    % [mol/(cm2 yr)] according non-bioturbated flux
                        swi.Fnonbio2 = swi.FPOM;
                        swi.C01 = 0.0;  % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
                        swi.C02 = 0.0;
                        swi.O20=129.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=18.6e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        % for Pacific
%                         swi.O20=10.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
%                         swi.NO30=80.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2400.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI

                    case 3 %  500m
                        % %         k1= 0.00001*0.174;                                                %TOC degradation rate constnat (1/yr)
                        % %         k2=0.174;  

                        res.bsd.wdepth=500.0;     % Dom was 600.0                       % water depth (m)
                        res.bsd = benthic_main(1, res.bsd.wdepth);
                        res.zTOC = benthic_zTOC(res.bsd);
                       	res.zTOC.k1= 0.00001*0.174;                                                %TOC degradation rate constnat (1/yr)
                       	res.zTOC.k2=0.174;  

                        res.bsd.rho_sed=2.5; %was 2.5                           % sediment density (g/cm3)
                        res.bsd.zbio=10.0;                              % bioturbation depth (cm)
                        res.bsd.zinf=100;                               %Inifinity (cm)
                        res.bsd.por=0.80; 
                        %bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        %bsd.w = %10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=8.1; %20.0;                         %temperature (degree C)
                        %swi.C01_nonbio= 2.64*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02_nonbio= 1.8*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                     	swi.FPOM = 357e-6;      % Total flux after Thullner et al. 2009 [mol/(cm2 yr)]
                        swi.Fnonbio1 = 0.00001*swi.FPOM;    % [mol/(cm2 yr)] according non-bioturbated flux
                        swi.Fnonbio2 = swi.FPOM;
                        swi.C01 = 0.0;  % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
                        swi.C02 = 0.0;
                        swi.O20=121.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=22.1e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        % for Pacific
%                         swi.O20=10.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
%                         swi.NO30=80.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2400.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
  
                    case 4 %  1000m
                        % %         k1= 0.00001*0.130;                                                %TOC degradation rate constnat (1/yr)
                        % %         k2=0.130;  

                        res.bsd.wdepth=1000.0;     % Dom was 600.0                       % water depth (m)
                        res.bsd = benthic_main(1, res.bsd.wdepth);
                        res.zTOC = benthic_zTOC(res.bsd);
                       	res.zTOC.k1= 0.00001*0.130;                                                %TOC degradation rate constnat (1/yr)
                       	res.zTOC.k2=0.130;  

                        res.bsd.rho_sed=2.5; %was 2.5                           % sediment density (g/cm3)
                        res.bsd.zbio=10.0;                              % bioturbation depth (cm)
                        res.bsd.zinf=100;                               %Inifinity (cm)
                        res.bsd.por=0.80; 
                        %bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        %bsd.w = %10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=5.8; %20.0;                         %temperature (degree C)
                        %swi.C01_nonbio= 2.64*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02_nonbio= 1.8*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                     	swi.FPOM = 228e-6;      % Total flux after Thullner et al. 2009 [mol/(cm2 yr)]
                        swi.Fnonbio1 = 0.00001*swi.FPOM;    % [mol/(cm2 yr)] according non-bioturbated flux
                        swi.Fnonbio2 = swi.FPOM;
                        swi.C01 = 0.0;  % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
                        swi.C02 = 0.0;
                        swi.O20=114.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=26.5e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        % for Pacific
%                         swi.O20=10.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
%                         swi.NO30=80.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2400.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
                        
                    case 5 %  2000m
                        % %         k1= 0.00001*0.0718;                                                %TOC degradation rate constnat (1/yr)
                        % %         k2=0.0718;  

                        res.bsd.wdepth=2000.0;     % Dom was 600.0                       % water depth (m)
                        res.bsd = benthic_main(1, res.bsd.wdepth);
                        res.zTOC = benthic_zTOC(res.bsd);
                       	res.zTOC.k1= 0.00001*0.0718;                                                %TOC degradation rate constnat (1/yr)
                       	res.zTOC.k2=0.0718;  

                        res.bsd.rho_sed=2.5; %was 2.5                           % sediment density (g/cm3)
                        res.bsd.zbio=10.0;                              % bioturbation depth (cm)
                        res.bsd.zinf=100;                               %Inifinity (cm)
                        res.bsd.por=0.80; 
                        %bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        %bsd.w = %10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=3.0; %20.0;                         %temperature (degree C)
                        %swi.C01_nonbio= 2.64*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02_nonbio= 1.8*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                     	swi.FPOM = 93e-6;      % Total flux after Thullner et al. 2009 [mol/(cm2 yr)]
                        swi.Fnonbio1 = 0.00001*swi.FPOM;    % [mol/(cm2 yr)] according non-bioturbated flux
                        swi.Fnonbio2 = swi.FPOM;
                        swi.C01 = 0.0;  % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
                        swi.C02 = 0.0;
                        swi.O20=116.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=31.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        % for Pacific
%                         swi.O20=10.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
%                         swi.NO30=80.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2400.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
 
                
                case 6 %  3500m
                        % %         k1= 0.00001*0.0296;                                                %TOC degradation rate constnat (1/yr)
                        % %         k2=1.0*0.0296;  

                        res.bsd.wdepth=3500.0;     % Dom was 600.0                       % water depth (m)
                        res.bsd = benthic_main(1, res.bsd.wdepth);
                        res.zTOC = benthic_zTOC(res.bsd);
                       	res.zTOC.k1= 0.00001*0.0296;                                                %TOC degradation rate constnat (1/yr)
                       	res.zTOC.k2=0.0296;  

                        res.bsd.rho_sed=2.5; %was 2.5                           % sediment density (g/cm3)
                        res.bsd.zbio=10.0;                              % bioturbation depth (cm)
                        res.bsd.zinf=100;                               %Inifinity (cm)
                        res.bsd.por=0.80; 
                        %bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        %bsd.w = %10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=1.5; %20.0;                         %temperature (degree C)
                        %swi.C01_nonbio= 2.64*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02_nonbio= 1.8*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                     	swi.FPOM = 24.3e-6;      % Total flux after Thullner et al. 2009 [mol/(cm2 yr)]
                        swi.Fnonbio1 = 0.0001*swi.FPOM;    % [mol/(cm2 yr)] according non-bioturbated flux
                        swi.Fnonbio2 = 1.0*swi.FPOM;
                        swi.C01 = 0.0;  % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
                        swi.C02 = 0.0;
                        swi.O20=135.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=31.6e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        % for Pacific
%                         swi.O20=10.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
%                         swi.NO30=80.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                        swi.Nitrogen=true;
                        swi.NH40=0.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2400.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI

                    case 7 %  5000m
                        % %         k1= 0.00001*0.0122;                                                %TOC degradation rate constnat (1/yr)
                        % %         k2=1.0*0.0122;  

                        res.bsd.wdepth=5000.0;     % Dom was 600.0                       % water depth (m)
                        res.bsd = benthic_main(1, res.bsd.wdepth);
                        res.zTOC = benthic_zTOC(res.bsd);
                       	res.zTOC.k1= 0.00001*0.0122;                                                %TOC degradation rate constnat (1/yr)
                       	res.zTOC.k2=0.0122;  

                        res.bsd.rho_sed=2.5; %was 2.5                           % sediment density (g/cm3)
                        res.bsd.zbio=10.0;                              % bioturbation depth (cm)
                        res.bsd.zinf=100;                               %Inifinity (cm)
                        res.bsd.por=0.8; 
                        %bsd.Dbio=5.2*(10.0^(0.7624-0.0003972*bsd.wdepth)); %0.5;
                        %bsd.w = %10.0.^(-0.87478367-0.00043512*bsd.wdepth)*3.3; % or check 0.42 for Reimers et al. 1996 as stated in the paper

                        %bottom water concentrations
                        swi.T=1.4; %20.0;                         %temperature (degree C)
                        %swi.C01_nonbio= 2.64*1e-2/12*bsd.rho_sed; % adjusted Test 2+4: 1.45* Test5: 35* Dom was 0.06*1e-2/12*bsd.rho_sed; %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                        %swi.C02_nonbio= 1.8*1e-2/12*bsd.rho_sed; % adjusted Test2+4: 6.5* Test5: 190* Dom was 0.06*1e-2/12*bsd.rho_sed;                                %TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
                     	swi.FPOM = 6.33e-6;      % Total flux after Thullner et al. 2009 [mol/(cm2 yr)]
                        swi.Fnonbio1 = 0.000001*swi.FPOM;    % [mol/(cm2 yr)] according non-bioturbated flux
                        swi.Fnonbio2 = 1.0*swi.FPOM;
                        swi.C01 = 0.0;  % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
                        swi.C02 = 0.0;
                        swi.O20=141.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
                        swi.NO30=31.6e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                         % for Pacific
%                         swi.O20=10.0e-9;   %was    300.0e-9  20              %O2  concentration at SWI (mol/cm^3)
%                         swi.NO30=80.0e-9;             % was 20.0e-9      %NO3 concentration at SWI (mol/cm^3)
                       swi.Nitrogen=true;
                        swi.NH40=0.0e-9;                                                %NH4 concentration at SWI (mol/cm^3)
                        swi.SO40=28000.0e-9;                                            %SO4 concentration at SWI (mol/cm^3)
                        %swi.SO40 = 100e-9;
                        swi.H2S0=0.0e-13;         %was 0.0e-9                            %H2S concentration at SWI (mol/cm^3)
                        swi.PO40=0.0e-9 ;%0.06e-8; % Dom was 1e-9;    % Sandra played with 3e-9                                              %PO4 concentration at SWI (mol/cm^3)
                        swi.Mflux0=365*0.2e-10; % Sandra played with 10e-9; ;   % = 7.3e-9    %flux of M to the sediment (mol/(cm2*yr))   TODO/CHECK: good value+right conversion? is from Slomp et al. 1996
                        swi.DIC0=2400.0e-9;                                             %DIC concentration at SWI (mol/cm^3)
                        swi.ALK0=2400.0e-9;                                             %ALK concentration at SWI (mol/cm^3)
                        swi.S0=35;                                                      %Salinity at SWI
                        
                    otherwise
                        error('unrecognized Obs  %g\n',pObs);
                end

            res.swi = swi;
            
            % calculate 
            res.zO2 = benthic_zO2(res.bsd, res.swi);           
            res.zNO3 = benthic_zNO3(res.bsd, res.swi);
            res.zSO4 = benthic_zSO4(res.bsd, res.swi);
            res.zNH4 = benthic_zNH4(res.bsd, res.swi);
            res.zH2S = benthic_zH2S(res.bsd, res.swi);
            res.zPO4_M = benthic_zPO4_M(res.bsd, res.swi);
            res.zDIC = benthic_zDIC(res.bsd, res.swi);
            res.zALK = benthic_zALK(res.bsd, res.swi);
   
            tic;
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
% %                res.zso4=res.zox;   % for test-case with just TOC & O2
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
            fprintf('POC1 wtpc SWI %g \n',  res.swi.C01*100*12/res.bsd.rho_sed);
            fprintf('POC2 wtpc SWI %g \n',  res.swi.C02*100*12/res.bsd.rho_sed);
            % calculate depth integrated OM degradation rates
            Cox_rate.Cox_total = res.zTOC.calcReac(0.0, res.bsd.zinf, 1, 1, res.bsd, swi, res);
            Cox_rate.Cox_aerobic = res.zTOC.calcReac(0.0, res.zox, 1, 1, res.bsd, swi, res);
            if(swi.Nitrogen)
                Cox_rate.Cox_denitr = res.zTOC.calcReac(res.zox, res.zno3, 1, 1, res.bsd, swi, res);
            end
                Cox_rate.Cox_sulfred = res.zTOC.calcReac(res.zno3, res.bsd.zinf, 1, 1, res.bsd, swi, res)

 %% Step 2: LOAD Observations   

%            switch Obs
%                     case 1  % OMEXDIA_2809_108m all solutes in Micromoles/litre
%                         str_date = '108m_OMEXDIA_1503_';
%                         data.TOC=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','Corg','C2:D24');     % in wt%
%                         data.O2=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','O2','C2:D126'); 
%                         data.NO3=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','NO3','C2:D25'); 
%                         data.NH4=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','NH4','C2:D25');
%                         data.PO4=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','PO4','C2:D25');
%                         data.DIC=xlsread('../Observations/OMEXDIA/5_PE138_99-06_108m.xlsx','PO4','C2:D25');
%             %            data.H2S=PW_data(:,[1 11]);
%                         
%                     case 2 % OMEXDIA_2809_2213m
%                         str_date = '2213m_OMEXDIA_2110_';
%                         data.TOC=xlsread('../Observations/OMEXDIA/2_PE121_98-4_2213m.xlsx','Corg','C2:D39');     % in wt%
%                         data.O2=xlsread('../Observations/OMEXDIA/2_PE121_98-4_2213m.xlsx','O2','C113:D1352'); % C2:D112
%                         data.NO3=xlsread('../Observations/OMEXDIA/2_PE121_98-4_2213m.xlsx','NO3','C2:D18'); 
%                         data.NH4=xlsread('../Observations/OMEXDIA/2_PE121_98-4_2213m.xlsx','NH4','C2:D18');
%                         data.PO4=xlsread('../Observations/OMEXDIA/2_PE121_98-4_2213m.xlsx','PO4','C2:D18');
%             %            data.SO4=PW_data(:,[1 10]);
%             %            data.H2S=PW_data(:,[1 11]);
%             
%                     case 3  % OMEXDIA_2809_3097m
%                          str_date = '3097m_OMEXDIA_2110_';
%                         data.TOC=xlsread('../Observations/OMEXDIA/3_PE138_99-14.xlsx','Corg','C2:D31');     % in wt%
%                         data.O2=xlsread('../Observations/OMEXDIA/3_PE138_99-14.xlsx','O2','C2:D391'); % C2:D112
%                         data.NO3=xlsread('../Observations/OMEXDIA/3_PE138_99-14.xlsx','NO3','C2:D26'); 
%                         data.NH4=xlsread('../Observations/OMEXDIA/3_PE138_99-14.xlsx','NH4','C2:D26');
%                         data.PO4=xlsread('../Observations/OMEXDIA/3_PE138_99-14.xlsx','PO4','C2:D25');
%             %            data.SO4=PW_data(:,[1 10]);
%             %            data.H2S=PW_data(:,[1 11]);
%             
%                     otherwise
%                         error('unrecognized Obs  %g\n',Obs);
%                 end
% 

%% Step 3: PLOT MODEL + Observations           
            if(pplot)
                benthic_test.plot_column(res, false, swi, 'Thullner')
            end
            

       
         end
         
         function plot_column_Observations(res, str_date, data, Obs)
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12) % plots 18
            
            bsd = res.bsd;
            zgrid = 0:0.1:bsd.zinf;
            
                figure;
                % TOC
                subplot(2,5,1)
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
                subplot(2,5,2)
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
    %            ylabel('Depth (cm)')
    %            title ('O2 (mol/cm^3)')

                % NO3

                for i=1:length(zgrid)
                    [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                end
                subplot(2,5,3)
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
                subplot(2,5,5)
                scatter(data.NH4(:,2).*1e-9, -data.NH4(:,1),'k','filled')
                hold on
                plot(NH4, -zgrid, 'b')
                xlim([0.0 2000e-9])     
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

                subplot(2,5,4)
                for i=1:length(zgrid)
                    [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                end
                if(Obs == 6 || Obs == 9)
                    scatter(data.SO4(:,2).*1e-9, -data.SO4(:,1),'k','filled')
                end
                hold on
                plot(SO4, -zgrid, 'b')
                xlim([0.0e-5 3.0e-5])     
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

                subplot(2,5,6)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                end
                if(Obs == 6 || Obs == 9)
                    scatter(data.H2S(:,2).*1e-9, -data.H2S(:,1),'k','filled')
                end
                hold on
                plot(H2S, -zgrid, 'b')
                xlim([0 4000e-9])
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
                subplot(2,5,7)
                for i=1:length(zgrid)                
                    [PO4(i), flxPO4(i), M(i), flxM(i), e_M(i), f_M(i), p_M(i), q_M(i), g_M(i), dedz_M(i), dfdz_M(i), dpdz_M(i), dqdz_M(i), dgdz_M(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
                end
                scatter(data.PO4(:,2).*1e-9, -data.PO4(:,1),'k','filled')
                hold on            
                plot(PO4, -zgrid, 'b')
                box on;
                xlim([0.0 100e-9])
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
                subplot(2,5,8)
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


                subplot(2,5,9)
                for i=1:length(zgrid)
                    [DIC(i), flxDIC(i)] = res.zDIC.calcDIC(zgrid(i), bsd, res.swi, res);
                end
                if(Obs == 1)
                    scatter(data.DIC(:,2).*1e-9, -data.DIC(:,1),'k','filled')
                end
                hold on
                plot(DIC, -zgrid, 'b')
    %            xlim([2.0e-5 3.0e-5])     
                box on;
    %                xlim([2.7e-5 swi.DIC0])     
          %      xlim([2.7e-5 swi.DIC0])     
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')     
                hold off
                ylim([-bsd.zinf 0.0]) 
                xlabel ('DIC (mol/cm^3)')
    %                ylabel('Depth (cm)')
    %            title ('SO4 (mol/cm^3)')

                subplot(2,5,10)
                for i=1:length(zgrid)
                    [ALK(i), flxALK(i)] = res.zALK.calcALK(zgrid(i), bsd, res.swi, res);
                end
                if(Obs == 6 || Obs == 9)
                    scatter(data.ALK(:,2).*1e-9, -data.ALK(:,1),'k','filled')
                end
                hold on
                plot(ALK, -zgrid, 'b')
%                xlim([0 10e-9])
                box on;
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')     
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')     
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')     
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')  
                ylim([-bsd.zinf 0.0]) 
                xlabel ('ALK (mol/cm^3)')
     %               ylabel('Depth (cm)')
        %            title ('ALK (mol/cm^3)')

    
           print('-depsc2', ['eps_output/' str_date 'PROFILES.eps']);

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
