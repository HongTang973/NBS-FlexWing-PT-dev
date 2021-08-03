classdef kappa_methods < handle
    %A collection of methods responsible for the calculation of varTheta and it's derivatives
    
    methods (Static = true)
        
        %method to return group 1 varTheta variables
        [KAPPA_I,dKAPPA_dt_I,dKAPPA_dqa_I_tr_Dim1x3xnsxnqa] = get_kappa_group1(...
            ...
            st,ct,sp,cp,dth_ds,dsi_ds,dph_ds,nszrs,R_W_WE,B_th,B_si,B_ph,dB_th,dB_si,dB_ph,st_sp,st_cp,ct_sp,ct_cp,yidx2,TD,...
            dzeta_a_dt_tr,d2th_dsdt,d2si_dsdt,d2ph_dsdt)
        
    end
    
    methods (Static = true, Access = private) %Group 1 Static Methods
        
        %dKAPPA_I
        %----------------------------------------------------------------------
        [dKAPPA_dt_I,dKAPPA_dq_I_tr] = get_dKAPPA_I(st,ct,sp,cp,dth_ds,dsi_ds,nszrs,R_W_WE_tr,B_th,B_si,B_ph,dB_th,dB_si,dB_ph,dzeta_a_dt_tr,d2th_dsdt,d2si_dsdt,d2ph_dsdt,st_sp,st_cp,ct_sp,ct_cp,yidx2,TD)
        
        %%%%%%%%%
        
        %         function val = spinCross(vec11,vec12,vec21,vec22,vec31,vec32)
        %             val = 1/2*(MultiProd_(getSkewMat(vec11),vec12)+MultiProd_(getSkewMat(vec21),vec22)+MultiProd_(getSkewMat(vec31),vec32));
        %         end
        
    end
    
end

