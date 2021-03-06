%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   HISTORY

% 24/09/2015    Dominik: started implementing 3rd TOC fraction - vectorized form not maintained

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef benthic_zPO4_M
    % Solve PO4
    
    properties                      
        qdispPO4=309.0528;          %PO4 diffusion coefficient in water (cm2/yr)
        adispPO4=12.2640;           %PO4 linear coefficient for temperature dependence (cm2/yr/oC)
        DPO41;                      %PO4 diffusion coefficient in bioturbated layer (cm2/yr)
        DPO42;                      %PO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
        
        KPO41=10.0;         %Adsorption coefficient in oxic layer (-)
        KPO42=1.3;          %Adsorption coefficient in anoxic layer (-)
        ksPO4=2.2*365;      %Rate constant for kinetic PO4 sorption (1/yr) was  0.12 fits 1.CASE; 2.2 fits 2. CASE DOM; from Nicolas was 0.5*365
        PO4s=1.0e-9;        %Equilibrium concentration for P sorption (mol/cm3)       was 1.5e-9
        kaPO4=0.004*365;	%Rate constant for authigenic P formation (1/yr)    DOM: was 0.004*365 from Nicolas
        PO4a=3.7e-9;        %Equilibrium concentration for authigenic P formation (mol/cm3) was 0.7e-9
        kmPO4=1.8e-6*365;	%Rate constant for Fe-bound P release upon Fe oxide reduction   DOM: was 0.0005*365 from Nicolas
        Minf=1.99e-6;       % asymptotic concentration for Fe-bound P (mol/cm3)      TODO/CHECK: good value? is from Slomp et al. 1996 Dom was 1.99e-6

   % OLD FROM NICOLAS     
%        FePFlux=0.01/1000;                                          %mol/m2/hr     Flux to the surface
%        FePdeepConc=1000e-9;                                        %Asymptotic concentration of FeP
%        FeSCt=FePdeepConc*(1.-por);                                 %convert FePdeepConc from /solid to /bulk

        % OM reactive terms in oxic layer
        reac1_ox;      
        reac2_ox;
        reac3_ox;
        % OM reactive terms in anoxic layer
        reac1_anox;
        reac2_anox;        
        reac3_anox;        
    end
    
    methods
        function obj = benthic_zPO4_M(bsd, swi)
            obj.DPO41=(obj.qdispPO4+obj.adispPO4*swi.T).*bsd.dispFactor+bsd.Dbio;            %PO4 diffusion coefficient in bioturbated layer (cm2/yr)
            obj.DPO42=(obj.qdispPO4+obj.adispPO4*swi.T).*bsd.dispFactor;                     %PO4 diffusion coefficient in non-bioturbated layer (cm2/yr)

  
             %reactive terms: OM degradation
             obj.reac1_ox=1/(1+obj.KPO41)*bsd.PC1;
             obj.reac2_ox=1/(1+obj.KPO41)*bsd.PC2;
             obj.reac3_ox=1/(1+obj.KPO41)*bsd.PC3;
             obj.reac1_anox=1/(1+obj.KPO42)*bsd.PC1;
             obj.reac2_anox=1/(1+obj.KPO42)*bsd.PC2;
             obj.reac3_anox=1/(1+obj.KPO42)*bsd.PC3;
         
        end
        
        function r = calc(obj, bsd, swi, r)

            
            % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
            
            % layer 1: 0 < z < zox, OM degradation (-) Sorption to sediment Fe-oxides (ktemp) 
            %           ls =  prepfg_l12_PO4M(   bsd, swi, r, reac1P,        reac2P,        reac3P,        kP                  ,     QP                          , zU, zL  ,   D1P    , D2P      , alphaP
            rPO4_M.ls1 = r.zTOC.prepfg_l12_PO4_M(bsd, swi, r, obj.reac1_ox, obj.reac2_ox, obj.reac3_ox, -obj.ksPO4/(1+obj.KPO41), obj.PO4s*obj.ksPO4/(1+obj.KPO41), 0, r.zox, obj.DPO41, obj.DPO42, 0, ...
                                                0, 0, bsd.Dbio, 0, (1/bsd.SD)*obj.ksPO4);
                                      % for M  kM,QM,     D1M ,D2M,   alphaM
                                
            % layer 2: zox < z < zinf, 
            % OM degradation (-) authigenic P formation (ktemp) (+) P desorption due to Fe-bound P release upon Fe oxide reduction

            %             ls =  prepfg_l12_PO4M( bsd, swi, r,     reac1P,        reac2P,        reac3P,        kP                ,     QP,                                zU,   zL,           D1P,       D2P,    alphaP
            rPO4_M.ls2 = r.zTOC.prepfg_l12_PO4_M(bsd, swi, r, obj.reac1_anox, obj.reac2_anox, obj.reac3_anox, -obj.kaPO4/(1+obj.KPO42), obj.PO4a*obj.kaPO4/(1+obj.KPO42), r.zox, bsd.zinf, obj.DPO41, obj.DPO42, bsd.SD*obj.kmPO4/(1+obj.KPO42), ...
                                                -obj.kmPO4, obj.kmPO4.*obj.Minf, bsd.Dbio, 0, 0);
                                  % for M           kM,       QM,              D1M,  D2M, alphaM)
                                  
            % Work up from the bottom, matching solutions at boundaries
            % Basis functions at bottom of layer 2 zinf
        %           calcfg_l12_PO4_M(obj, z, bsd, swi, res,     reac1P,        reac2P,              kP,                   QtempP,                          alphaP,                 ls , ....
        %    arguments for M              kM,        QM,           alphaM)               
            [ e2_zinf_P, dedz2_zinf_P, f2_zinf_P, dfdz2_zinf_P, g2_zinf_P, dgdz2_zinf_P, p2_zinf_P, dpdz2_zinf_P, q2_zinf_P, dqdz2_zinf_P, ...
                e2_zinf_M, dedz2_zinf_M, f2_zinf_M, dfdz2_zinf_M, g2_zinf_M, dgdz2_zinf_M, p2_zinf_M, dpdz2_zinf_M, q2_zinf_M, dqdz2_zinf_M] ...
                = r.zTOC.calcfg_l12_PO4_M(bsd.zinf, bsd, swi, r, obj.reac1_anox, obj.reac2_anox, obj.reac3_anox, -obj.kaPO4/(1+obj.KPO42), obj.PO4a*obj.kaPO4/(1+obj.KPO42), bsd.SD*obj.kmPO4/(1+obj.KPO42), rPO4_M.ls2, -obj.kmPO4, obj.kmPO4.*obj.Minf, 0);
            % calcfg_l12_PO4_M(obj,            z,   bsd, swi, res,     reac1P,         reac2P,         reac3P,            kP,                      QtempP,                       alphaP,                     ls,           kM,           QtempM,        alphaM)
                        
            % Match at zox, layer 1 - layer 2 (continuity and flux)                        
            % basis functions at bottom of layer 1
            [ e1_zox_P, dedz1_zox_P, f1_zox_P, dfdz1_zox_P, g1_zox_P, dgdz1_zox_P, p1_zox_P, dpdz1_zox_P, q1_zox_P, dqdz1_zox_P, ...
                e1_zox_M, dedz1_zox_M, f1_zox_M, dfdz1_zox_M, g1_zox_M, dgdz1_zox_M, p1_zox_M, dpdz1_zox_M, q1_zox_M, dqdz1_zox_M]...
                = r.zTOC.calcfg_l12_PO4_M(r.zox, bsd, swi, r,  obj.reac1_ox, obj.reac2_ox, obj.reac3_ox, -obj.ksPO4/(1+obj.KPO41), obj.PO4s*obj.ksPO4/(1+obj.KPO41), 0, rPO4_M.ls1, 0, 0, (1/bsd.SD)*obj.ksPO4);
                % calcfg_l12_PO4_M(obj, z, bsd, swi, res,   reac1P,       reac2P,               reac3P,  kP,                       QtempP,                        alphaP, ls,      kM, QtempM, alphaM)
            
            % ... and top of layer 2
            [ e2_zox_P, dedz2_zox_P, f2_zox_P, dfdz2_zox_P, g2_zox_P, dgdz2_zox_P, p2_zox_P, dpdz2_zox_P, q2_zox_P, dqdz2_zox_P, ...
                e2_zox_M, dedz2_zox_M, f2_zox_M, dfdz2_zox_M, g2_zox_M, dgdz2_zox_M, p2_zox_M, dpdz2_zox_M, q2_zox_M, dqdz2_zox_M] ...
                = r.zTOC.calcfg_l12_PO4_M(r.zox, bsd, swi, r, obj.reac1_anox, obj.reac2_anox, obj.reac3_anox, -obj.kaPO4/(1+obj.KPO42), obj.PO4a*obj.kaPO4/(1+obj.KPO42), bsd.SD*obj.kmPO4/(1+obj.KPO42), rPO4_M.ls2, -obj.kmPO4, obj.kmPO4.*obj.Minf, 0);
              % calcfg_l12_PO4_M(obj, z, bsd, swi, res,     reac1P,         reac2P,                 reac3P,          kP,                        QtempP,                       alphaP,                      ls,         kM,           QtempM,    alphaM)
            
            % match solutions at zox - continuous concentration and flux
            % organize the data in matrices and let matlab do the calculation
            %  |x1        |   | A_l |      | y1        | | A_r|    |z1|    always PO4 continuity  
            %  |    .     |   | B_l |      |    .      | | B_r|    |z2|    always PO4 flux  
            %  |      .   |   | C_l |   =  |      .    | | C_r|  + |z3|    always M continuity  
            %  |       x16|   | D_l |      |        y16| | D_r|    |z4|    always M flux  

            % discontinuity constants
            Vb = 0;
            Fb = 0;
            
            if(r.zox <= bsd.zbio)   % 1. CASE: 4 int const. in each layer
                X = [e1_zox_P, f1_zox_P, p1_zox_P, q1_zox_P; ...
                     dedz1_zox_P, dfdz1_zox_P, dpdz1_zox_P, dqdz1_zox_P; ...
                     e1_zox_M, f1_zox_M, p1_zox_M, q1_zox_M; ...
                     dedz1_zox_M, dfdz1_zox_M, dpdz1_zox_M, dqdz1_zox_M];
                Y = [e2_zox_P, f2_zox_P, p2_zox_P, q2_zox_P; ...
                     dedz2_zox_P, dfdz2_zox_P, dpdz2_zox_P, dqdz2_zox_P; ...
                     e2_zox_M, f2_zox_M, p2_zox_M, q2_zox_M; ...
                     dedz2_zox_M, dfdz2_zox_M, dpdz2_zox_M, dqdz2_zox_M];
                Z = [g2_zox_P-g1_zox_P + Vb; ... 
                     dgdz2_zox_P - dgdz1_zox_P + Fb - bsd.w.*Vb; ...
                     g2_zox_M-g1_zox_M + Vb; ... 
                     dgdz2_zox_M - dgdz1_zox_M + Fb - bsd.w.*Vb];      
                 case_flag=1; 
            else    % 2. CASE: 3 int const. in each layer
                X = [e1_zox_P, f1_zox_P, p1_zox_P; ...
                     dedz1_zox_P, dfdz1_zox_P, dpdz1_zox_P; ...
                     e1_zox_M, f1_zox_M, p1_zox_M; ...
                     dedz1_zox_M, dfdz1_zox_M, dpdz1_zox_M];
                Y = [e2_zox_P, f2_zox_P, p2_zox_P; ...
                     dedz2_zox_P, dfdz2_zox_P, dpdz2_zox_P; ...
                     e2_zox_M, f2_zox_M, p2_zox_M; ...
                     dedz2_zox_M, dfdz2_zox_M, dpdz2_zox_M];
                Z = [g2_zox_P-g1_zox_P + Vb; ... 
                     dgdz2_zox_P - dgdz1_zox_P + Fb - bsd.w.*Vb; ...
                     g2_zox_M-g1_zox_M + Vb; ... 
                     dgdz2_zox_M - dgdz1_zox_M + Fb - bsd.w.*Vb];  
                 case_flag=2; 
            end
                    
            [zox.C, zox.D] = benthic_utils.matchsoln_PO4_M(X, Y, Z,case_flag);

                       
            % Solution at swi, top of layer 1
            [ e1_0_P, dedz1_0_P, f1_0_P, dfdz1_0_P, g1_0_P, dgdz1_0_P, p1_0_P, dpdz1_0_P, q1_0_P, dqdz1_0_P, ...
                e1_0_M, dedz1_0_M, f1_0_M, dfdz1_0_M, g1_0_M, dgdz1_0_M, p1_0_M, dpdz1_0_M, q1_0_M, dqdz1_0_M]...
                = r.zTOC.calcfg_l12_PO4_M(0, bsd, swi, r, obj.reac1_ox, obj.reac2_ox, obj.reac3_ox, -obj.ksPO4/(1+obj.KPO41), obj.PO4s*obj.ksPO4/(1+obj.KPO41), 0,  rPO4_M.ls1, 0, 0, (1/bsd.SD)*obj.ksPO4);
            % calcfg_l12_PO4_M(obj, z, bsd, swi, res, reac1P,              reac2P,          reac3P,          kP,                       QtempP,             alphaP,   ls,    kM, QtempM, alphaM)
          
            % transform to use coeffs from l2
            % Now find 'transformed' basis functions such that in layer 1, O2 = A_2*et + B_2*ft + gt  (ie layer 1 soln written in terms of layer 2 coeffs A_2, B_2)
                            
            EFPQ_P = [e1_0_P; f1_0_P; p1_0_P; q1_0_P];   
            dEFPQdz_P = [dedz1_0_P; dfdz1_0_P; dpdz1_0_P; dqdz1_0_P];
            EFPQ_M = [e1_0_M; f1_0_M; p1_0_M; q1_0_M];
            dEFPQdz_M = [dedz1_0_M; dfdz1_0_M; dpdz1_0_M; dqdz1_0_M];
            [EFPQ_P, g1_0_P, dEFPQdz_P, dgdz1_0_P, EFPQ_M, g1_0_M, dEFPQdz_M, dgdz1_0_M] ...
            = benthic_utils.xformsoln_PO4_M(EFPQ_P, EFPQ_M, dEFPQdz_P, dEFPQdz_M, g1_0_P, g1_0_M,dgdz1_0_P, dgdz1_0_M, zox.C, zox.D);
        
% % % %     VERSION WITH SINGLE VALUES (NO MATRICES)            
% % % %             [e1_0_P, f1_0_P, p1_0_P, q1_0_P, g1_0_P, dedz1_0_P, dfdz1_0_P, dpdz1_0_P, dqdz1_0_P, dgdz1_0_P, e1_0_M, f1_0_M, p1_0_M, q1_0_M, g1_0_M, dedz1_0_M, dfdz1_0_M, dpdz1_0_M, dqdz1_0_M, dgdz1_0_M] ...
% % % %             = benthic_utils.xformsoln_PO4_M(e1_0_P, f1_0_P, p1_0_P, q1_0_P, e1_0_M, f1_0_M, p1_0_M, q1_0_M, dedz1_0_P, dfdz1_0_P, dpdz1_0_P, dqdz1_0_P, dedz1_0_M, dfdz1_0_M, dpdz1_0_M, dqdz1_0_M, g1_0_P, g1_0_M,dgdz1_0_P, dgdz1_0_M, zox.C, zox.D);
            
            % Solve for APO4M, BPO4M, CPO4M, DPO4M given boundary conditions (expressed in terms of transformed basis fns, layer 2 A, B, C, D)
            % BC for PO4: zero flux at z=oo AND concentration equal at z=0
            % BC for M: Known flux at z=0 AND concentration = Minf at z=oo
            
            % APO4M*dedz2_zinf_P   +  BPO4M*dfdz2_zinf_P  + CPO4M*dpdz2_zinf_P   +  DPO4M*dqdz2_zinf_P + dgdz2_zinf_P = 0;
            % APO4M*e1_0_P         +  BPO4M*f1_0_P        + CPO4M*e1_0_P         +  DPO4M*f1_0_P       + g1_0_P       = swi.PO40;
            % APO4M*e2_zinf_M      +  BPO4M*f2_zinf_M     + CPO4M*p2_zinf_M      +  DPO4M*q2_zinf_M    + g2_zinf_M    = obj.Minf;
            % APO4M*dedz1_0_M      +  BPO4M*dfdz1_0_M     + CPO4M*dpdz1_0_M      +  DPO4M*dqdz1_0_M    + dgdz1_0_M    = swi.Mflux0;          
                        
            % | dedz2_zinf_P dfdz2_zinf_P dpdz2_zinf_P dqdz2_zinf_P|  |APO4M|     | - dgdz2_zinf_P       |
            % |    e1_0_P      f1_0_P        p1_0_P      q1_0_P    |  |BPO4M|     | swi.PO40 - g1_0_P    |
            % | e2_zinf_M    f2_zinf_M     p2_zinf_M    q2_zinf_M  |  |CPO4M|   = | obj.Minf - g2_zinf_M |
            % | dedz1_0_M    dfdz1_0_M     dpdz1_0_M    dqdz1_0_M  |  |DPO4M|     |swi.Mflux0 - dgdz1_0_M|
            
            if(case_flag==2)    % DEAL WITH VARIABLES SHORT FROM LAYER BELOW IN 2. CASE 
                EFPQ_P(3)=p1_0_P;
                dEFPQdz_P(3)=dpdz1_0_P;
                EFPQ_P(4)=q1_0_P;
                dEFPQdz_P(4)=dqdz1_0_P;

                EFPQ_M(4)=q1_0_M;
                dEFPQdz_M(4)=dqdz1_0_M;
            end           
        
            X = [dedz2_zinf_P, dfdz2_zinf_P, dpdz2_zinf_P, dqdz2_zinf_P; ...
                 EFPQ_P(1), EFPQ_P(2), EFPQ_P(3), EFPQ_P(4); ...
                 e2_zinf_M, f2_zinf_M, p2_zinf_M, q2_zinf_M; ...
                 dEFPQdz_M(1), dEFPQdz_M(2), dEFPQdz_M(3), dEFPQdz_M(4)];
            Y = [-dgdz2_zinf_P; ...
                 swi.PO40 - g1_0_P; ...
                 obj.Minf - g2_zinf_M; ...
                 swi.Mflux0 - dgdz1_0_M];
             
            [ rPO4_M.A2, rPO4_M.B2, rPO4_M.C2, rPO4_M.D2]  = benthic_utils.solve2eqn_PO4_M(X,Y);
            
            
% % % %      % calculate PO4 conc and flux at zinf (flux is ZERO anyway)
% % % %             % CHECK/TODO: Why not calculate concentration at zinf?
% % % %             conczPO4 = rPO4_M.A2.*e2_zinf_P + rPO4_M.B2.*f2_zinf_P + rPO4_M.C2.*p2_zinf_P + rPO4_M.D2.*q2_zinf_P + g2_zinf_P;
% % % %             D = (bsd.zinf <= bsd.zbio).*obj.DPO41 + (zPO4 > bsd.zbio).*obj.DPO42;
% % % %             flxzPO4 = D.*(rPO4_M.A3.*dedz3_zPO4+rPO4_M.B3.*dfdz3_zPO4 + dgdz3_zPO4);        % includes 1/por ie flux per (cm^2 pore area)
            

            % flux PO4 at swi - DO include por so this is per cm^2 water column area
            r.flxswi_P = bsd.por.*obj.DPO41.*(rPO4_M.A2.*dEFPQdz_P(1)+rPO4_M.B2.*dEFPQdz_P(2) + rPO4_M.C2.*dEFPQdz_P(3)+rPO4_M.D2.*dEFPQdz_P(4) + dgdz1_0_P);   % NB: use A2, B2, C2, D2 as these are _xformed_ layer 1 basis functions
            
            % flux M at swi - DO include por so this is per cm^2 water column area
            r.flxswi_M = bsd.por.*obj.DPO41.*(rPO4_M.A2.*dEFPQdz_M(1)+rPO4_M.B2.*dEFPQdz_M(2) + rPO4_M.C2.*dEFPQdz_M(3)+rPO4_M.D2.*dEFPQdz_M(4) + dgdz1_0_M);   % NB: use A2, B2, C2, D2 as these are _xformed_ layer 1 basis functions
            
            
            % save coeffs for layer 1          
            L1_P = zox.C*[ rPO4_M.A2; rPO4_M.B2; rPO4_M.C2; rPO4_M.D2]+zox.D;
            
            
            rPO4_M.A1 = L1_P(1);
            rPO4_M.B1 = L1_P(2);
            rPO4_M.C1 = L1_P(3);
            rPO4_M.D1 = L1_P(4);
                      
            
            r.rPO4_M = rPO4_M;
            
        end
                     
        
               
        function FPO4 = calcFPO4(obj, zPO4, bsd, swi, r)
           % Calculate PO4 consumption below zPO4, by organic matter and indirectly via methane oxidation 
           
            tmpreac1    = bsd.PO4C.*bsd.gammaCH4;
            tmpreac2    = bsd.PO4C.*bsd.gammaCH4;
            tmpreac3    = bsd.PO4C.*bsd.gammaCH4;
       
            FPO4 = r.zTOC.calcReac(zPO4, bsd.zinf, tmpreac1, tmpreac2, tmpreac3, bsd, swi, r);
            % TODO confirm (1-bsd.por)*  has been added (to k1 & k2 ?)
        end
        
        
                                          
        function [PO4, flxPO4, M, flxM] = calcPO4_M(obj, z, bsd, swi, r)
            % Calculate PO4 concentration and flux at depth z from solution
            
                rPO4_M = r.rPO4_M;
                if z <= bsd.zbio
                    D_P = obj.DPO41;
                    D_M = bsd.Dbio;
                else
                    D_P = obj.DPO42;
                    D_M = 0;
                end
                
                if (z <= r.zox)   % layer 1                    
                    [ e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, ...
                    e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M]...
                    = r.zTOC.calcfg_l12_PO4_M(z, bsd, swi, r,  obj.reac1_ox, obj.reac2_ox, obj.reac3_ox, obj.ksPO4/(1+obj.KPO41), obj.PO4s*obj.ksPO4/(1+obj.KPO41), 0, rPO4_M.ls1, 0, 0, (1/bsd.SD)*obj.ksPO4);
                    % calcfg_l12_PO4_M(obj, z, bsd, swi, res,   reac1P,       reac2P,          reac3P,          kP,                     QtempP,                alphaP,      ls,   kM, QtempM,  alphaM)
                    
                    
                    PO4     = r.rPO4_M.A1.*e_P + r.rPO4_M.B1.*f_P + r.rPO4_M.C1.*p_P + r.rPO4_M.D1.*q_P + g_P;
                    flxPO4  = D_P.*(r.rPO4_M.A1.*dedz_P+r.rPO4_M.B1.*dfdz_P + r.rPO4_M.C1.*dpdz_P+r.rPO4_M.D1.*dqdz_P + dgdz_P);
                    M     = r.rPO4_M.A1.*e_M + r.rPO4_M.B1.*f_M + r.rPO4_M.C1.*p_M + r.rPO4_M.D1.*q_M + g_M;
                    flxM  = D_M.*(r.rPO4_M.A1.*dedz_M+r.rPO4_M.B1.*dfdz_M + r.rPO4_M.C1.*dpdz_M+r.rPO4_M.D1.*dqdz_M + dgdz_M);

                else        % layer 2
                    [ e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P, ...
                    e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M]...
                     = r.zTOC.calcfg_l12_PO4_M(z, bsd, swi, r, obj.reac1_anox, obj.reac2_anox, obj.reac3_anox, obj.kaPO4/(1+obj.KPO42), obj.PO4a*obj.kaPO4/(1+obj.KPO42), bsd.SD*obj.kmPO4/(1+obj.KPO42), rPO4_M.ls2, -obj.kmPO4, obj.kmPO4.*obj.Minf, 0);
                    % calcfg_l12_PO4_M(obj, z, bsd, swi, res,     reac1P,         reac2P,          reac3P,          kP,                        QtempP,                     alphaP,                    ls,      kM,           QtempM,    alphaM)
                   
                    
                    PO4     = r.rPO4_M.A2.*e_P + r.rPO4_M.B2.*f_P + r.rPO4_M.C2.*p_P + r.rPO4_M.D2.*q_P + g_P;
                    flxPO4  = D_P.*(r.rPO4_M.A2.*dedz_P+r.rPO4_M.B2.*dfdz_P + r.rPO4_M.C2.*dpdz_P+r.rPO4_M.D2.*dqdz_P + dgdz_P);
                    M     = r.rPO4_M.A2.*e_M + r.rPO4_M.B2.*f_M + r.rPO4_M.C2.*p_M + r.rPO4_M.D2.*q_M + g_M;
                    flxM  = D_M.*(r.rPO4_M.A2.*dedz_M+r.rPO4_M.B2.*dfdz_M + r.rPO4_M.C2.*dpdz_M+r.rPO4_M.D2.*dqdz_M + dgdz_M);
                end
                
        end
        
    end
    
end

