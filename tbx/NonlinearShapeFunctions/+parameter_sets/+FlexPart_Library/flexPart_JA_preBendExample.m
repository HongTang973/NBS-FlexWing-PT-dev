function O = flexPart_JA_preBendExample(master_object,parent_object,name)

O = NBS_flexPart_nonlinear(master_object,name,'Parent',parent_object);

%% --------------------------------------------------------------------------
%Example using a pre-bend geometry

L = 10;
ns = 201;
s = permute(linspace(0,L,ns),[1 3 2]);
ds = s.*0; ds(1:end-1) = diff(s); ds(end) = ds(end-1);

%coordinate system
%X from hub in upstream direction
%Y in radial direction
%Z complete right hand set

%a projection of the blade profile, lets choose a 1-cos curve
YZ_projection = 1-cos(s.*3/4*(2*pi/L));
dYZ_projection = s.*0; 
dYZ_projection(1:end-1) = diff(YZ_projection); 
dYZ_projection(end) = dYZ_projection(end-1);

%a twist distribution, choose here some nonlinear function going from -45 to 45 degrees
twist_deg = ((s./L).^0.8-0.5)*90;

%information required

%initial Euler angle distributions
ph_0 = twist_deg*pi/180;
th_0 = atan(dYZ_projection./ds);
si_0 = s.*0;

%--------------------------------------------------------------------------

O.R_W_WE = eye(3);
get_kappa(O,s,th_0,si_0,ph_0);

%--------------------------------------------------------------------------
nAnodes = 17;
Anode_Skew_Factor = 1.0;
Anode_distr = linspace(1,0,nAnodes).^Anode_Skew_Factor;
s_aero_ = (1 - Anode_distr)*L;
O.s_aero = permute(s_aero_,[1 3 2]);
O.isAero = true;

O.w = ones(1,1,ns)*1;
O.h = ones(1,1,ns)*0.2;
%--------------------------------------------------------------------------
EIxx = 2e4; EIzz = 4e6; GJ  = 1e4;
O.StiffnessMatrix = [
    EIxx 0    0;
    0    GJ   0;
    0    0    EIzz];
O.StiffnessMatrix(4:6,4:6) = eye(3)*1e3;
%--------------------------------------------------------------------------
damping_factor = 0.04;
O.DampingMatrix = damping_factor*O.StiffnessMatrix;
%--------------------------------------------------------------------------
I_tau = 0.1;
O.I_varTheta_ps_I = [
    I_tau 0 0;
    0 I_tau 0;
    0 0 I_tau];
O.I_varTheta_discrete_I = 0;
%--------------------------------------------------------------------------
O.mps = 0.75;
O.msDiscrete = 0;
%--------------------------------------------------------------------------
O.c = ones(1,1,ns);
O.beam_cntr = ones(1,1,ns)*0.5;
%--------------------------------------------------------------------------
O.AICs = sqrt(1-(0.5*(O.s_aero(2:end)+O.s_aero(1:end-1))).^2./L^2)*0+1;
O.CrossSectionProfiles = 'NACA0012';
%--------------------------------------------------------------------------

end


function kappa = get_kappa(O,s,th,si,ph)

ds = diff(s);
dth_ds = th.*0; dsi_ds = si.*0; dph_ds = ph.*0;
dth_ds(1:end-1) = diff(th)./ds;
dsi_ds(1:end-1) = diff(si)./ds;
dph_ds(1:end-1) = diff(ph)./ds;

%-----------------------------------
kappa_x_IE = dsi_ds.*cos(th).*sin(ph) + dth_ds.*cos(ph);                   %(1)x(1)x(ns)
tau_IE     = dph_ds - dsi_ds.*sin(th);                                     %(1)x(1)x(ns)
kappa_z_IE =-dsi_ds.*cos(th).*cos(ph) + dth_ds.*sin(ph);                   %(1)x(1)x(ns)

kappa = [kappa_x_IE ; tau_IE ; kappa_z_IE];
%-----------------------------------

O.s = s;
O.th0 = th; O.si0 = si; O.ph0 = ph;
O.dth_ds0 = dth_ds; O.dsi_ds0 = dsi_ds; O.dph_ds0 = dph_ds;
O.KAPPA_0_I = kappa;

end
