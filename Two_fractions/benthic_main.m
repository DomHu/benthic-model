classdef benthic_main < handle
    % Global properties for benthic model
    %
    % 
    
    properties
        
        ncl;                                    % number of sediment columns
        usescalarcode = true;                   % use scalar code 
        
        tol_const = 1e-18;                      % non-zero constant to avoid numerical errors (e.g. division by zero)
        
        %sediment characteristics
        rho_sed=2.5; %was 2.5                           % sediment density (g/cm3)
        wdepth=100.0; %3575.5999999999999;     % Dom was 600.0                       % water depth (m)
        w;                                      % burial velocity  (cm/yr)
        z0  = 0;                                % surface
        zbio=10.0; %10.0;                              % bioturbation depth (cm)       
        
        zinf=100;                               %Inifinity (cm)
        %zinf = 1000;
        %zlow=100;                              
        Dbio; % was 29.06.2016 =5.2*(10.0^(0.7624-0.0003972*obj.wdepth)); % Dom was 3;                                 %bioturbation coefficient (cm2/yr)
        Dunbio=0.01;
        por=0.85;      %was 0.85                          %porosity (-)
        tort=3.0;                               %tortuosity (-)
        irrigationFactor=1.0;                   %irrigation factor (-)
        dispFactor;                             %dispersion factor (-)
        
        %stoichiometric factors
        X_C;                                    % Carbon Redfield stoichiometry
        Y_N;                                    % Nitrogen Redfield stoichiometry
        Z_P;                                    % Phosphorous Redfield stoichiometry
        SD;                                     %volume factor solid->dissolved phase
        OC;                                     %O2/C (mol/mol)
        NC1;                                    %N/C first TOC fraction (mol/mol)
        NC2;                                    %N/C second TOC fraction (mol/mol)
        PC1;                                    %P/C first TOC fraction (mol/mol)
        PC2;                                    %P/C second TOC fraction (mol/mol)
        SO4C;                                   %SO4/C (mol/mol)
        O2H2S;                                  % O2/H2S ratio for oxidation of H2S (mol/mol)
        DICC1;                                  %DIC/C until zSO4 (mol/mol)
        DICC2;                                  %DIC/C below zSO4 (mol/mol)
        MC;                                     %CH4/C (mol/mol)
        gamma=0.05;                                %fraction of NH4 that is oxidised in oxic layer
        gammaH2S=0.05;                           %fraction of H2S that is oxidised in oxic layer
        gammaCH4=0.99;                           %fraction of CH4 that is oxidised at SO4
        satSO4=0.0;                               % SO4 saturation
        NO3CR;                                  % NO3 consumed by Denitrification
        ALKROX;                                              % Aerobic degradation                     
        ALKRNIT;                                                 % Nitrification    
        ALKRDEN;                                           % Denitrification
        ALKRSUL;                                             % Sulfato reduction
        ALKRH2S;                                                 % H2S oxydation (CHECK THIS VALUE!!!)
        ALKRMET;                                             % Methanogenesis
        ALKRAOM;                                                  % AOM

        
        zoxgf = 0.0;              % was 0.1              % cm, rolloff NH4, H2S oxidation for small zox depth
        
        % Diagnostic output from root finder
        %fzerooptions;
        %fzerooptions = optimset('Display','iter');
        %fzerooptions = optimset('Display','final');
        %fzerooptions = optimset('TolX',0.001);
        fzerooptions = optimset('TolX',100*eps);
    end
    
    methods
        function obj = benthic_main(ncl, wdepth)
            if nargin > 0
                obj.ncl = ncl;                 
                obj.wdepth = obj.wdepth*ones(1,obj.ncl);
                obj.zbio = obj.zbio*ones(1,obj.ncl);
                obj.Dbio = obj.Dbio*ones(1,obj.ncl);
                obj.zinf = obj.zinf*ones(1,obj.ncl);
                obj.z0 = obj.z0*ones(1,obj.ncl);
            else
                obj.ncl = 1;
            end
            
            if nargin > 1
                    obj.wdepth = wdepth;
            end
            
            obj.usescalarcode = (obj.ncl == 1);
            
            obj.w=benthic_main.sedrate(obj.wdepth); 
            obj.Dbio=benthic_main.biorate(obj.wdepth);
            obj.dispFactor=obj.por.^(obj.tort-1.0).*obj.irrigationFactor;                 %dispersion factor (-)
            obj.SD=(1-obj.por)./obj.por;   % Sandra played with 1.0
            
            obj.X_C=106;                                                    % Carbon Redfield stoichiometry
            obj.Y_N=16;                                                     % Nitrogen Redfield stoichiometry
            obj.Z_P=1;                                                      % Phosphorous Redfield stoichiometry
            obj.OC=(138/106)*obj.SD; % 1.0*obj.SD; %                      	%O2/C (mol/mol)
            obj.NC1= 0.067*obj.SD; %16.0/106.0*obj.SD; % Bohlen 0.067*obj.SD;    %   % 0.0; % 0.1509                                     %N/C first TOC fraction: 16/106 (mol/mol)
            obj.NC2= 0.067*obj.SD; %16.0/106.0*obj.SD; % Bohlen 0.067*obj.SD;    %0.1509*obj.SD; %0.13333*obj.SD;                                         %N/C second TOC fraction (mol/mol)            
            obj.PC1=1/106*obj.SD;  % was 0.0094 Sandra played with 1e-20;             %P/C first TOC fraction  1/106 (mol/mol)
            obj.PC2=1/106*obj.SD;  % was 0.0094 Sandra played with  1e-20;            %P/C second TOC fraction 1/106 (mol/mol)
            obj.SO4C=(138.0/212.0)*obj.SD; % 0.5*obj.SD;                                                %SO4/C (mol/mol)
            obj.O2H2S=2.0;
            obj.DICC1=1.0*obj.SD;                                               %DIC/C until zSO4 (mol/mol)
            obj.DICC2=0.5*obj.SD;                                               %DIC/C below zSO4 (mol/mol)
            obj.MC=0.5*obj.SD;                                                  %CH4/C (mol/mol)
            obj.NO3CR=(94.4/106)*obj.SD;                                        % NO3 consumed by Denitrification
            
%           For implicit Nitrogen (e.g. coupled to GENIE):  
            obj.ALKROX=-(obj.Y_N)/obj.X_C*obj.SD;                        % -(obj.Y_N+2*obj.Z_P)/obj.X_C*obj.SD was +15      % Aerobic degradation                     
            obj.ALKRNIT=0.0;  % no ALK                                             % Nitrification    
            obj.ALKRDEN=0.0;       %  462/530 was 93.4;    % Denitrification
%       For explicit Nitrogen:
%            obj.ALKROX=(obj.Y_N-2*obj.Z_P)/obj.X_C*obj.SD;              % Aerobic degradation
%            obj.ALKRNIT=-2.0;  % no ALK                                             % Nitrification    
%            obj.ALKRDEN=(4*obj.X_C+3*obj.Y_N-10*obj.Z_P)/(5*obj.X_C)*obj.SD;       %  462/530 was 93.4;    % Denitrification

%           For implicit Nitrogen (e.g. coupled to GENIE):
            obj.ALKRSUL=(obj.X_C+obj.Y_N)/obj.X_C*obj.SD;                % 122/106  was 15      	% Sulfato reduction
%           For explicit Nitrogen:     
%            obj.ALKRSUL=(obj.X_C+obj.Y_N-2*obj.Z_P)/obj.X_C*obj.SD; %(obj.X_C+obj.Y_N)/obj.X_C*obj.SD;                % 120/106  was 15      	% Sulfato reduction
            
            obj.ALKRH2S=-2.0;      % no ALK                                         % H2S oxydation (CHECK THIS VALUE!!!)
            obj.ALKRMET= (obj.Y_N-2*obj.Z_P)/obj.X_C;      % no ALK                	% 14/106 was 14       	% Methanogenesis
            obj.ALKRAOM= 2.0;        % no ALK                                       	% AOM

        end
    end
    
    methods(Static)
        
        function w = sedrate(wdepth)
            % sedimentation rate, cm/yr
            w = 10.0.^(-0.87478367-0.00043512*wdepth)*3.3; % 0.2668 is at 500m this is of  10.0.^(-0.87478367-0.00043512*wdepth)*3.3; % 0.03; 
            % w = 5.0e-4; %12.0E-003; 
        end
        
         function Dbio = biorate(wdepth)
            % bioturbation coeff, cm^2/yr
            Dbio= 5.2*(10.0^(0.7624-0.0003972*wdepth)); %5.2*(10.0^(0.7624-0.0003972*wdepth)); % Dom was 3.0;                                 %bioturbation coefficient (cm2/yr)
            %Dbio=0.1; %0.08; %10.0;
        end       
        
        
        
    end
    
end

