classdef benthic_main < handle
    % Global properties for benthic model
    %
    % 
    
    properties
        
        ncl;                                    % number of sediment columns
        usescalarcode = true;                   % use scalar code 
        
        %sediment characteristics
        rho_sed=2.5;                            % sediment density (g/cm3)
        wdepth=600;                            % water depth (m)
        w;                                      % burial velocity  (cm/yr)
        z0  = 0;                                % surface
        zbio=5;            % Sandra played with 0.0                    % bioturbation depth (cm)       
        
        zinf=100;                               %Inifinity (cm)
        %zinf = 1000;
        %zlow=100;                              
        Dbio=3;                                 %bioturbation coefficient (cm2/yr)
        por=0.8;                                %porosity (-)
        tort=3.0;                               %tortuosity (-)
        irrigationFactor=1.0;                   %irrigation factor (-)
        dispFactor;                             %dispersion factor (-)
        
        %stoichiometric factors
        SD;                                     %volume factor solid->dissolved phase
        OC;                                     %O2/C (mol/mol)
        NC1;                                    %N/C first TOC fraction (mol/mol)
        NC2;                                    %N/C second TOC fraction (mol/mol)
        PC1;                                    %P/C first TOC fraction (mol/mol)
        PC2;                                    %P/C second TOC fraction (mol/mol)
        SO4C;                                   %SO4/C (mol/mol)
        DICC1;                                  %DIC/C until zSO4 (mol/mol)
        DICC2;                                  %DIC/C below zSO$ (mol/mol)
        MC;                                     %CH4/C (mol/mol)
        gamma=1;                                %fraction of NH4 that is oxidised in oxic layer
        gammaH2S=0.8;                           %fraction of H2S that is oxidised in oxic layer
        gammaCH4=0.9;                           %fraction of CH4 that is oxidised at SO4
        satSO4=0;                               % SO4 saturation
        NO3CR;                                  % NO3 consumed by Denitrification
        
        zoxgf = 0.1;                            % cm, rolloff NH4, H2S oxidation for small zox depth
        
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
            obj.dispFactor=obj.por.^(obj.tort-1.0).*obj.irrigationFactor;                 %dispersion factor (-)
            obj.SD=(1-obj.por)./obj.por;   % Sandra played with 1.0
            
            obj.OC=1.0*obj.SD;                                                  %O2/C (mol/mol)
            obj.NC1=0.1509*obj.SD;                                              %N/C first TOC fraction (mol/mol)
            obj.NC2=0.13333*obj.SD;                                             %N/C second TOC fraction (mol/mol)            
            obj.PC1=0.0094*obj.SD;  % Sandra played with 1e-20;                                            %P/C first TOC fraction (mol/mol)
            obj.PC2=0.0094*obj.SD;  % Sandra played with  1e-20;                                         %P/C second TOC fraction (mol/mol)
            obj.SO4C=0.5*obj.SD;                                                %SO4/C (mol/mol)
            obj.DICC1=1.0*obj.SD;                                               %DIC/C until zSO4 (mol/mol)
            obj.DICC2=0.5*obj.SD;                                               %DIC/C below zSO$ (mol/mol)
            obj.MC=0.5*obj.SD;                                                  %CH4/C (mol/mol)
            obj.NO3CR=(94.4/106)*obj.SD;                                        % NO3 consumed by Denitrification

        end
    end
    
    methods(Static)
        
        function w = sedrate(wdepth)
            % sedimentation rate, cm/yr
            w=10.0.^(-0.87478367-0.00043512*wdepth)*3.3; % 0.03; 
        end
        
        
        
        
        
    end
    
end

