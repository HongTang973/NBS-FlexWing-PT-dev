classdef ct_methods
    
    methods (Static = true)
        
        M_W = map_WEtoW(M_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD)
        
        M_G = map_WtoG(M_W,R_G_W)

        M_A = map_WtoA(M_W,R_A_W)

        M_G = map_AtoG(M_A,R_G_A)

        M_W_flat = map_WEtoW_stackDim2(M_WE_flat,Rs_W_WE_tr,Rv_W_WE_tr,TD)

    end
end

