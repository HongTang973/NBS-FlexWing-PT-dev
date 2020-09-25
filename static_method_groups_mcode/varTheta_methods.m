classdef varTheta_methods < handle
    %A collection of methods responsible for the calculation of varTheta and it's derivatives
    
    methods (Static = true)
        
        %method to return group 1 varTheta variables
        function [dvarTheta_dqa_G_3x1xnsxnqa,dvarTheta_dt_G,d2varTheta_dt2_G_star] = get_varTheta_group1(...
                ...
                ss, cs, ey_WE, dE_dt_W,...
                ...
                nszrs,ns,...
                R_W_WE,R_G_W,...
                dvarTheta_dt_root_G,dvarTheta_dt_root_G_skew,...
                ...
                dth_dt,dsi_dt,dph_dt,...
                B_th_tr,B_si_tr,B_ph_tr,...
                ...
                d2varTheta_dt2_root_G_star)
            
            %--------------------------------------------------------------
            dthVec_WE = [cs ;-ss ; nszrs];
            dthVec_G = mult_Anm1_Bmpz(R_G_W*R_W_WE,dthVec_WE);
            %dthVec_W = mult_Anm1_Bmpz(R_W_WE,dthVec_WE);
            %dthVec_G = ct_methods.map_WtoG(dthVec_W,R_G_W);
            %--------------------------------------------------------------
            dsiVec_W =-R_W_WE(:,3);                                        % = R_W_WE*[0;0;-1]
            dsiVec_G = mult_Anmz_Bmp1(R_G_W,dsiVec_W);
            if size(dsiVec_G,3)==1, dsiVec_G = repmat(dsiVec_G,[1 1 ns]); end
            %--------------------------------------------------------------
            dphVec_WE = ey_WE;
            dphVec_G = mult_Anm1_Bmpz(R_G_W*R_W_WE,dphVec_WE);
            %--------------------------------------------------------------
            
            dzeta_a_dt_tr = [dth_dt,dsi_dt,dph_dt];
            [dvarTheta_dqa_G,dvarTheta_dt_G] = varTheta_methods.get_dvarTheta(dthVec_G,dsiVec_G,dphVec_G,B_th_tr,B_si_tr,B_ph_tr,dzeta_a_dt_tr,dvarTheta_dt_root_G);
            dvarTheta_dqa_G_3x1xnsxnqa = permute(dvarTheta_dqa_G,[1 4 3 2]);
%             if ~isempty(dvarTheta_dqe_G)
%                 dvarTheta_dqr123_G = zeros(3);
%                 dvarTheta_dqr456_G = eye(3);
%                 dvarTheta_dqr_G = [dvarTheta_dqr123_G,dvarTheta_dqr456_G];
%             else
%                 dvarTheta_dqr_G = zeros(3,0);
%             end
            
            dthdsiVec_WE = [-ss ;-cs ; nszrs];
            dthdsiVec_W = mult_Anm1_Bmpz(R_W_WE,dthdsiVec_WE);
           %dthdsiVec_G = ct_methods.map_WtoG(dthdsiVec_W,R_G_W);
            
            dey_dt_W = dE_dt_W(:,2,:);
            
            d2varTheta_dt2_W_star = ...
                bsxfun(@times,dth_dt.*dsi_dt,dthdsiVec_W) +...
                bsxfun(@times,dph_dt,dey_dt_W);
            
            d2varTheta_dt2_G_star_temp = ...
                mult_Anm1_Bmpz(dvarTheta_dt_root_G_skew,dvarTheta_dt_G)...
              + ct_methods.map_WtoG(d2varTheta_dt2_W_star,R_G_W);
            d2varTheta_dt2_G_star = bsxfun(@plus,d2varTheta_dt2_root_G_star,d2varTheta_dt2_G_star_temp);
        
        end
        
    end
    
    
    methods (Static = true, Access = private) %Group 1 Static Methods
        
        %dvarTheta
        %----------------------------------------------------------------------
        function [dvarTheta_dq_G,dvarTheta_dt_G] = get_dvarTheta(dthVec_G,dsiVec_G,dphVec_G,B_th_tr,B_si_tr,B_ph_tr,dzeta_a_dt_tr,dvarTheta_dt_root_G)
            dvarTheta_dq_G = [bsxfun(@times,B_th_tr,dthVec_G) bsxfun(@times,B_si_tr,dsiVec_G) bsxfun(@times,B_ph_tr,dphVec_G)];
            %---------------------------
            dvarTheta_dt_G = ...
                bsxfun(@plus , sum([dthVec_G,dsiVec_G,dphVec_G].*[dzeta_a_dt_tr;dzeta_a_dt_tr;dzeta_a_dt_tr],2) , dvarTheta_dt_root_G);
        end
        %----------------------------------------------------------------------
        
    end
    
end

