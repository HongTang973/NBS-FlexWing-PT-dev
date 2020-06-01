function O = flexPart_windTurbineBlade(master_object,parent_object,name)

O = NBS_flexPart_nonlinear(master_object,name,'Parent',parent_object);

%% Taken from old code

        %wing turbine blade based on 5-MW reference offshore turbine
        load('TurbineProperties.mat'); %-> TurbineProperties_5MW
        
        L = TurbineProperties_5MW.s.vals(end)-TurbineProperties_5MW.s.vals(1);
        ns = 201;
        nAnodes = 201;
        s = permute(linspace(0,L,ns),[1 3 2]);
        del_s = s(2:end)-s(1:end-1);
        s_aero = permute(linspace(0,L,nAnodes),[1 3 2]);
        nszrs = s*0;
        
        span_samples = TurbineProperties_5MW.s.vals-TurbineProperties_5MW.s.vals(1);
        sample1 = sampleMat(span_samples,s);
        sampleMat_tr = sample1.';
        
        %properties per span
        x_ = 1-s/max(s);
       %EI1 = simple_smooth(permute(sampleMat_tr*TurbineProperties_5MW.EI1.vals,[2 3 1]));
        aEI1 = round([5.929391124711572;
           7.429623753729656;
           0.679124993171499;
           6.311480672605765;
          -2.638636441905949;
           2.599628620255143;
          -1.814730492485513;
          -1.709136041845677],4);
        aEI16 = aEI1(6) + aEI1(7) + aEI1(8);
        EI1 = exp( 12.04 + aEI1(2)*x_.^(1/aEI1(1)) + (aEI1(3)*x_ + aEI1(4)*x_.^2 + aEI1(5)*x_.^3 + aEI16*x_.^4) );
        
       %EI2 = simple_smooth(permute(sampleMat_tr*TurbineProperties_5MW.EI2.vals,[2 3 1]));
        aEI2 = round([4.035116264012688;
           7.269835194203548;
           0.155323950694095;
          -0.245882804627211;
           0.307151836152231;
           0.379466651682743;
          -0.122716886834057;
           0.209872194808633],4);
        aEI26 = aEI2(6) + aEI2(7) + aEI2(8);
        EI2 = exp( 15.43 + aEI2(2)*x_.^(1/aEI2(1)) + (aEI2(3)*x_ + aEI2(4)*x_.^2 + aEI2(5)*x_.^3 + aEI26*x_.^4) );
           
       %GJ  = simple_smooth(permute(sampleMat_tr*TurbineProperties_5MW.GJ.vals ,[2 3 1]));
        aGJ = round([3.546636560169235;
           6.886384671488536;
          -1.280936676900210;
           2.089223502020646;
           2.885835978006972;
          -0.948737779772829;
           0.726988215239715;
           0.253905252008049],4);
        aGJ6 = aGJ(6) + aGJ(7) + aGJ(8);
        GJ = exp( 12.15 + aGJ(2)*x_.^(1/aGJ(1)) + (aGJ(3)*x_ + aGJ(4)*x_.^2 + aGJ(5)*x_.^3 + aGJ6*x_.^4) );
        
        %Ixx = permute(sampleMat_tr*TurbineProperties_5MW.Ixx,[2 3 1]);
        %Iyy = permute(sampleMat_tr*TurbineProperties_5MW.Iyy,[2 3 1]);
        %Izz = permute(sampleMat_tr*TurbineProperties_5MW.Izz,[2 3 1]);
        
        mps_coeff(1) = -1.411872057825304;
        mps_coeff(2) =  17.861710218789344;
        mps_coeff(3) = -21.588031423815117;
        mps_coeff(4) = -9.732997830074641;
        mps_coeff(5) =  17.044945799643749;
        mps_coeff = round(mps_coeff,4);
        x_ = (61.7-s)/61.7;
        mps = 250*(x_).^(0.5).*(1+mps_coeff(1)*x_.^2 + mps_coeff(2)*x_.^3 + mps_coeff(3)*x_.^4 + mps_coeff(4)*x_.^5 + mps_coeff(5)*x_.^6);
        
        %mps = simple_smooth(permute(sampleMat_tr*TurbineProperties_5MW.mps.vals ,[2 3 1]));
        
        EIGJ = [EI1 nszrs nszrs ; nszrs GJ nszrs ; nszrs nszrs EI2];
        EIGJ_Damp = 0.05*EIGJ;
        I_varTheta_ps = [nszrs nszrs nszrs;nszrs nszrs+10 nszrs;nszrs nszrs nszrs];%temp torsional inertia value, not representative!
        
        sweep = 0; dihedral = 0;

%         rotor_diameter = 126;%m not taking into account pre-cone shortening
%         hub_diameter = 3;%m
%         pre_cone = 2.5;%deg %overestimate of pre-cone to account for neglecting pre-curvature
%         shaft_tilt = 5;%deg
%         rotor_mass = 110000;%kg
        
        h = nszrs + 0.1;
        w = simple_smooth(TurbineProperties_5MW.chord.vals(:,2).'*sampleMat(TurbineProperties_5MW.chord.vals(:,1).'-TurbineProperties_5MW.s.vals(1),squeeze(s).'));
        
        
        %pre curvature, curve fitted from data
        x_ = (s-30.75)/18.02;
        p(1) = -0.0021779;
        p(2) = 0.00047782;
        p(3) = 0.0080149;
        p(4) = -0.0074728;
        p(5) = -0.0031189;
        p(6) = 0.039748;
        p(7) = 0.029954;
        
        dy = polyval(p,x_).*(atan(1*(s-14))/pi+0.5);
        %dy = p1*x_.^6 + p2*x_.^5 + p3*x_.^4 + p4*x_.^3 + p5*x_.^2 + p6*x_ + p7;
        %[~,b] = find(dy<0); dy(1:max(b)) = 0;
        
        x_ = (s-30.75)/18.02;
        p(1) = -0.14539;
        p(2) = 0.30529;
        p(3) = 0.80577;
        p(4) = -1.7149;
        p(5) = -1.2031;
        p(6) = 2.6739;
        p(7) = 0.73322;
        p(8) = -0.71825;
        p(9) = -5.3307;
        p(10) = 6.5963;
        
        twist = (polyval(p,x_)-13.5).*(atan(1*(s-10))/pi+0.5)+13.5;
        %twist(1:11) = 13.52;
        
        th0 = -atan(dy);
        dth_ds0_ = (th0(2:end)-th0(1:end-1))./del_s;
        dth_ds0 = cat(3,dth_ds0_(1),(dth_ds0_(2:end)+dth_ds0_(1:end-1))/2,dth_ds0_(end));
        ph0 = -twist*pi/180;
        dph_ds0_ = (ph0(2:end)-ph0(1:end-1))./del_s;
        dph_ds0 = cat(3,dph_ds0_(1),(dph_ds0_(2:end)+dph_ds0_(1:end-1))/2,dph_ds0_(end));
        si0 = nszrs;
        dsi_ds0 = nszrs;
        
        kappa_x0 = dsi_ds0.*cos(th0).*sin(ph0) + dth_ds0.*cos(ph0);        %(1)x(1)x(ns)
        tau_0    = dph_ds0 - dsi_ds0.*sin(th0);                            %(1)x(1)x(ns)
        kappa_z0 =-dsi_ds0.*cos(th0).*cos(ph0) + dth_ds0.*sin(ph0);        %(1)x(1)x(ns)
        
        KAPPA_0_I = [kappa_x0;tau_0;kappa_z0];
        %th0 = th0*0; dth_ds0 = dth_ds0*0;
        %ph0 = ph0*0; dph_ds0 = dph_ds0*0;
        
        c = permute(w,[1 3 2]);
        aero_cntr = simple_smooth(permute(sampleMat_tr*TurbineProperties_5MW.aeroCentre.vals,[2 3 1]));
        %beam_cntr = simple_smooth(permute(0.5 - (sampleMat_tr*TurbineProperties_5MW.massCentre.vals)./w.',[2 3 1]));
        
        p(1) = -3.264e-3;
        p(2) = 5.961e-5;
        p(3) = -4.356e-7;
        
        beam_cntr = 0.5 + p(1)*s + p(2)*s.^2 + p(3)*s.^3;
        
        R_W_WE = eye(3);
        V = 24;
        rho = 0.0881;
        alpha_0 = -20;
        tip_force_global  = [0;0;0];
        tip_force_local   = [0;0;0];
        tip_moment_global = [0;0;0];
        tip_moment_local  = [0;0;0];
        grav_acc = 9.807*0;
        gravVec_G = [0;0;-1];
        aerodynamics = 0;
        AICs = sqrt(1-(s_aero/L).^2);
        Pvec_appliedLocal_I = bsxfun(@times,AICs,[0;0;1]*20000);%*1200 to take pre curvature out of the blade   ,   20000 large load test case
        q0 = zeros(1,10); dq0 = q0*0;
        p0 = zeros(1,10); dp0 = p0*0;
        r0 = zeros(1,10); dr0 = r0*0;
        IC = [q0,p0,r0,dq0,dp0,dr0];
        n = length(q0); m = length(p0); o = length(r0);
        nmo = [n m o];
        %//////////////////
        shape_format_bend = 'chebyshev_1st';
        shape_format_twist = 'chebyshev_1st';
        shape_BCs_bend  = [0 1;1 1;1 1];
        shape_BCs_twist = [0 1;1 1;1 1];
        
        plotBounds = [-L/2 L/2;0 L;-L/4 L*3/4];
        
        af1 = TurbineProperties_5MW.airfoil_sections.DU40_A17;
        af2 = TurbineProperties_5MW.airfoil_sections.NACA64_A17;
        transform_upper = sampleMat(-af1.vals_upper(:,1).',-af2.vals_upper(:,1).');
        transform_lower = sampleMat(af1.vals_lower(:,1).',af2.vals_lower(:,1).');
        AeroProfiles = {
            8.333-1.5,cos(linspace(0,pi,26))/2+0.5,sin(linspace(0,pi,26))/2,cos(linspace(pi,0,26))/2+0.5,-sin(linspace(pi,0,26))/2;%cylindrical section
            11.75-1.5,af1.vals_upper(:,1).'*transform_upper,af1.vals_upper(:,2).'*transform_upper,af1.vals_lower(:,1).'*transform_lower,af1.vals_lower(:,2).'*transform_lower;
            40.45-1.5,af1.vals_upper(:,1).'*transform_upper,af1.vals_upper(:,2).'*transform_upper,af1.vals_lower(:,1).'*transform_lower,af1.vals_lower(:,2).'*transform_lower;
            L,af2.vals_upper(:,1).',af2.vals_upper(:,2).',af2.vals_lower(:,1).',af2.vals_lower(:,2).'};

        
%% Build the flexPart object from the above variables

O.s = s;
O.s_aero = s_aero;

O.w = w;
O.h = w;
%--------------------------------------------------------------------------
O.StiffnessMatrix = EIGJ;
%--------------------------------------------------------------------------
damping_factor = 0.02*10;
O.DampingMatrix = damping_factor*O.StiffnessMatrix;
%--------------------------------------------------------------------------
O.I_varTheta_ps_I = I_varTheta_ps
O.I_varTheta_discrete_I = 0;
%--------------------------------------------------------------------------
O.mps = mps;
O.msDiscrete = 0;
%--------------------------------------------------------------------------
%O.sweep = sweep;
%O.dihedral = 0;
O.R_W_WE = R_W_WE;
O.KAPPA_0_I = KAPPA_0_I;
%--------------------------------------------------------------------------
O.c = c;
%O.aero_cntr = aero_cntr;
O.beam_cntr = beam_cntr;
%--------------------------------------------------------------------------
O.AICs = sqrt(1-(0.5*(O.s_aero(2:end)+O.s_aero(1:end-1))).^2./L^2)*0+1;
O.CrossSectionProfiles = AeroProfiles;
%--------------------------------------------------------------------------
O.isAero = true;

end