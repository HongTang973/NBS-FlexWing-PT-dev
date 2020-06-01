function O = rigidPart_ExampleHTP(master_object,parent_object,name)

O = NBS_rigidPart(master_object,name,'Parent',parent_object);

%--------------------------------------------------------------------------
L = 6;

O.s = reshape([-L/2 0 L/2],1,1,[]); ns = numel(O.s);
O.s_aero = O.s; nAnode = numel(O.s_aero);
O.isAero = true;
    
O.Gamma_W = O.s.*[0;1;0];

O.w = cat(3,1,1.5,1);
O.c = O.w;
O.h = ones(1,1,ns)*0.2;
O.beam_cntr = ones(1,1,ns)*0.5;
O.AICs = ones(1,1,nAnode-1);
%--------------------------------------------------------------------------
O.root_idx = 1;
O.Mu = 10;
O.rMu_W = [0;L/2;0];
%--------------------------------------------------------------------------
I_tau = 1;
O.I_varTheta = [
    I_tau 0     0;
    0 I_tau 0;
    0 0     I_tau];
%--------------------------------------------------------------------------
O.CrossSectionProfiles = 'NACA0012';
O.root_idx = 1;

end

