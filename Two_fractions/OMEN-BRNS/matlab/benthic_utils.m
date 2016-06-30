classdef benthic_utils
    % Utility functions for solution matching etc
    
    
    
    methods(Static)
       
        
        
        function [a, b, c, d, e, f] = matchsoln(E_l, F_l, G_l, dEdx_l, dFdx_l, dGdx_l, ...
                                                E_r, F_r, G_r, dEdx_r, dFdx_r, dGdx_r, ...
                                                Vb, Db)
            % Match two solutions at a boundary:  
            % 'left' solution   y_l(x) = A_l*E_l(x) + B_l*F_l(x) + G_l(x)
            % 'right' solution  y_r(x) = A_r*E_r(x) + B_r*F_r(x) + G_r(x)
            %
            % (Dis)continuity conditions at boundary:
            %                   y_r(xb)    = y_l(xb)     + Vb
            %                   dydx_r(xb) = dydx_l(xb)  + Db
            %
            % Find a,b,c,d,e,f such that:
            %         | A_l |   =  | a  b | | A_r|  + |e|  
            %         | B_l |      | c  d | | B_r|    |f|
            
            alden = dFdx_l.*E_l - F_l.*dEdx_l;
            a     = (dFdx_l.*E_r     - F_l.*dEdx_r)./alden;
            b     = (dFdx_l.*F_r     - F_l.*dFdx_r)./alden;
            e     = (F_l.*(dGdx_l - dGdx_r + Db) + dFdx_l.*(-G_l + G_r - Vb))./alden; 
            blden = dEdx_l.*F_l      -   E_l.*dFdx_l;
            c     = (dEdx_l.*E_r     -   E_l.*dEdx_r)./blden;
            d     = (dEdx_l.*F_r     -   E_l.*dFdx_r)./blden;
            f     = (E_l.*(dGdx_l-dGdx_r+Db) + dEdx_l.*(-G_l+G_r - Vb))./blden;
            
            %
        end

        function [Et, Ft, Gt, dEtdx, dFtdx, dGtdx] = xformsoln(E, F, G, dEdx, dFdx, dGdx, ...
                                                               a , b , c , d , e ,f)
            % Find 'transformed' soln such that in layer l, 
            %    y_l = A_r*et + B_r*ft + gt  
            % (ie l soln written in terms of r solution coefficents A_r, B_r)                                               
                                                           
            Et      = a.*E    + c.*F;
            dEtdx   = a.*dEdx + c.*dFdx;
            Ft      = b.*E    + d.*F;
            dFtdx   = b.*dEdx + d.*dFdx;
            Gt      = G      + e.*E      + f.*F;
            dGtdx   = dGdx   + e.*dEdx   + f.*dFdx;
        end
        
        
        function [ x, y]      = solve2eqn(a, b, c, d, e, f)
            % Find soln of
            % | a    b |  |x|   = | e |
            % | c    d |  |y|     | f |
            
            det = a.*d-b.*c;
            x    =  (e.*d-b.*f)./det;
            y    =  (a.*f-e.*c)./det;
                        
        end

        
                
        function [C,D] = matchsoln_PO4_M(X, Y, Z, case_flag)
            % Match four solutions at a boundary:  
            %  for PO4
            % 'left' solution   y_l(z) = A_l*E_l(z) + B_l*F_l(z) + C_l*P_l(z) + D_l*Q_l(z) + G_l(z)
            % 'right' solution  y_r(z) = A_r*E_r(z) + B_r*F_r(z) + C_r*P_r(z) + D_r*Q_r(z) + G_r(z)
            %
            %  and the same for M
            %
            % (Dis)continuity conditions at boundary:
            %                   y_r(xb)    = y_l(xb)     + Vb
            %                   dydx_r(xb) = dydx_l(xb)  + Db
            
            %         | A_l |         | A_r|      
            %         | B_l |         | B_r|    
            %     X   | C_l |   =  Y  | C_r|  + Z  
            %         | D_l |         | D_r|               
            
            
            
            %
            % Find C and D such that:
            %         | A_l |         | A_r|      
            %         | B_l |         | B_r|    
            %         | C_l |   =  C  | C_r|  +  D  
            %         | D_l |         | D_r|               
            %    
            C = X\Y;
            D = X\Z;
            
            if(case_flag==2)            % add line and column of zero for missing Q
                C = [C;0,0,0];          % add line
                %C = [C [0;0;0;1]];     % add column OR JUST ZEROS ?????
                C = [C [0;0;0;0]];      % SD - just zeros
                D = [D;0];
            end
            
            %
        end
   
        function [EFPQ_P_t, G_P_t, dEFPQ_P_t, dG_P_t, EFPQ_M_t, G_M_t, dEFPQ_M_t, dG_M_t] = xformsoln_PO4_M(EFPQ_P, EFPQ_M, dEFPQdz_P, dEFPQdz_M, g_P, g_M,dgdz_P, dgdz_M, C, D)
            % Find 'transformed' soln such that in layer l, 
            %    y_l = A_r*et + B_r*ft + gt  
            % (ie l soln written in terms of r solution coefficents A_r, B_r)   
            %
            % NOW: save values in matrices! 
            
            
            EFPQ_P_t = C'*EFPQ_P;
            dEFPQ_P_t = C'*dEFPQdz_P;
            G_P_t = D'*EFPQ_P+g_P;
            dG_P_t = D'*dEFPQdz_P+dgdz_P;
            
            EFPQ_M_t = C'*EFPQ_M;
            dEFPQ_M_t = C'*dEFPQdz_M;
            G_M_t = D'*EFPQ_M+g_M;
            dG_M_t = D'*dEFPQdz_M+dgdz_M;

            
        end
        
% % % %     VERSION WITH SINGLE VALUES (NOT MATRICES)        
% % % %         function [e1_0_P_t, f1_0_P_t, p1_0_P_t, q1_0_P_t, g1_0_P_t, dedz1_0_P_t, dfdz1_0_P_t, dpdz1_0_P_t, ...
% % % %                     dqdz1_0_P_t, dgdz1_0_P_t, e1_0_M_t, f1_0_M_t, p1_0_M_t, q1_0_M_t, g1_0_M_t, dedz1_0_M_t, dfdz1_0_M_t, dpdz1_0_M_t, dqdz1_0_M_t, dgdz1_0_M_t] ...
% % % %                 = xformsoln_P_tO4_M(e1_0_P, f1_0_P, p1_0_P, q1_0_P, e1_0_M, f1_0_M, p1_0_M, q1_0_M, dedz1_0_P, dfdz1_0_P, dpdz1_0_P, dqdz1_0_P, ...
% % % %                                     dedz1_0_M, dfdz1_0_M, dpdz1_0_M, dqdz1_0_M, g1_0_P, g1_0_M,dgdz1_0_P, dgdz1_0_M, C, D)
% % % %             % Find 'transformed' soln such that in layer l, 
% % % %             %    y_l = A_r*et + B_r*ft + gt  
% % % %             % (ie l soln written in terms of r solution coefficents A_r, B_r)   
% % % %             %
% % % %             % NOW: save values in matrices! 
% % % %             
% % % % 
% % % %                       
% % % %         end
        

        function [ A, B, C, D]      = solve2eqn_PO4_M(X, Y)
            % Find soln of
            % | x1  .  . x4|  |A|     | y1 |
            % |     .      |  |B|     | y2 |
            % |       .    |  |C|   = | y3 |
            % | .       x16|  |D|     | y4 |
            
            Z = X\Y;
            
            A = Z(1);
            B = Z(2);
            C = Z(3);
     %       D = Z(4);                        
        end
        
       
    end
    
end

