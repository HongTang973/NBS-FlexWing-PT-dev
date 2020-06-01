function O = rigidPart_ExampleFuselage(master_object,parent_object,name)

O = NBS_rigidPart(master_object,name,'Parent',parent_object);

%--------------------------------------------------------------------------
    L = 25;
O.s = reshape(linspace(-L*1/5,L*4/5,50),1,1,[]); ns = numel(O.s);
O.s_aero = reshape([0 L],1,1,[]);
O.isAero = false;
    
O.Gamma_W = O.s.*[0;1;0];

    temp = linspace(1,0,ns);
O.w = reshape(sin(temp.^2*pi)+0.3,1,1,[]);
O.h = reshape(sin(temp.^2*pi)+0.3,1,1,[]);
%--------------------------------------------------------------------------
O.root_idx = 1;
O.Mu = 10;
O.rMu_W = [0;L/2*0+1;0];
%--------------------------------------------------------------------------
I_tau = 1;
O.I_varTheta = [
    I_tau 0     0;
    0 I_tau 0;
    0 0     I_tau];
%--------------------------------------------------------------------------
O.CrossSectionProfiles = 'circ';
O.root_idx = 1;

end