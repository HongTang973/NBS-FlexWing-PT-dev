
classdef NBS_flexPart_nonlinear < handle

    properties
        NBS_Master
        Parent
        root_idx = 1                                                       %spanwise index of the wing reference point (this point has zero velocity in the wing and aircraft coordinate systems)
        connection_idx_ParentObj = 1
        address
        partName
        QOI_Container                                                      %handle to flexPart QOI_Container
    end

    properties%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
     %...properties..... .........................
      %..................derived_properties.......
       %|---------------|-------------------------
        s                                                                  %[m] spanwise evaluations points along the wing
        %               del_s                                              %[m] distance between neighbouring s points
        %               ns                                                 %[-] number of s points
        %               L                                                  %[m] total length of the wing
        s_ms %currently unused                                             %[m] spanwise locations of discrete mass points
        s_aero                                                             %[m] spanwise locations of aero panel boundaries
        %               del_s_aero
        %               s_aeroMid                                          %[m] spanwise location of aero mid points
        %               nsAp                                               %[-] number of spanwise aero panels
                        Apn_idx                                            %[-] linear coefficients for sampling quantities at aero node points
                        Apm_idx                                            %[-] linear coefficients for sampling quantities at aero mid points
        sAp_idx_global
        StiffnessMatrix                                                    %[Nm/rad] structural stiffness matrix
        %               StiffnessMatrix_kappaHalf
        %               StiffnessMatrix_shearHalf
        DampingMatrix                                                      %[Nms/rad] structural damping matrix
        damping_factor
        %               DampingMatrix_kappaHalf
        %               DampingMatrix_shearHalf
        I_varTheta_ps_I                                                    %[kg m] rotational structural inertia matrix (per span)
                        I_varTheta_continuous_I                            %[kg m] rotational structural inertia matrix (per span)
        I_varTheta_discrete_I = 0                                          %[kg] optional ns point inertias to be used alongside I_varTheta_continuous_I in making up the inertial definition of the wing
                        I_varTheta                                         %[kg m^2] rotational inertia matrix for each del_s width element
        mps                                                                %[kg/m] structural mass per span
                        msContinuous                                       %[kg] ns point masses calculed from mass per span
        msDiscrete = 0                                                     %[kg] ns discrete point masses to be used alongside mps in making up the mass definition of the wing
                        ms                                                 %[kg] total ns point masses deriving from both discrete and continuous sources
        massOffsetContinuous_I                                             %[m] mass offsets for msContinuous in the intrinsic reference system
        massOffsetDiscrete_I                                               %[m] mass offsets for msDiscrete in the intrinsic reference system
                        massOffset_I                                       %[m] mass offsets for ms, taken from the beam reference line and expressed in the intrinsic coordinate system
                        massOffsetFlag                                     %[true/false] mass offset flag, true if non-zero offsets supplied
        c                                                                  %[m] wing chord lengths in direction of aerodynamic panels
                        c_pm                                               %[m] wing chord sampled at aero panel mid points
                        c_pn                                               %[m] wing chord sampled at aero panel nodal points
                        ApWidths_pm                                        %[m] aero panel widths in the eyAp direction
        w                                                                  %[m] wing width in the ex direction
        h                                                                  %[m] wing thicknesses
        aero_cntr                                                          %[-] assumed aero centres for each panel (percentage of chord)
        aero_cntr_pm

        aeroData = struct();
        beam_cntr                                                          %[-] beam centre locations (percentage of chord)
                        beam_cntr_pm                                       %[-] beam centre sampled at aero panel mid points
                        beam_cntr_pn                                       %[-] beam centre sampled at aero panel nodal points

        alpha_A                                                            %[deg] angle of attack of aero panels in the aircraft system
        alpha_I                                                            %[deg] angle of attack of aero panels in the intrinsic system (not used if alpha_A specified)
        ASD_CELL

        reflectedPart = false
        EAp_W
                        EAp_I                                              %[-] body fixed vector triad aligned with the aerodynamic panels and expressed in the intrinsic reference frame
                        EAp_I_pm                                           %[-] EAp_I sampled at the aero panel mid points
                        EAp_I_pn                                           %[-] EAp_I sampled at the aero panel nodal points

        varTheta0_W
        E0_W                                                               %initial E ortientation triad in the [W] system
        KAPPA_0_I = 0                                                      %[rad/m] intial undeformed wing curvature
                                                                           %--------Euler Angle Datum--------
        th0 = 0                                                            %[rad] initial theta distribution
        dth_ds0 = 0                                                        %[rad/m] initial spanwise theta derivative
        si0 = 0                                                            %[rad] initial psi distribution
        dsi_ds0 = 0                                                        %[rad/m] initial spanwise psi derivative
        ph0 = 0                                                            %[rad] initial phi distribution
        dph_ds0 = 0                                                        %[rad/m] initial spanwise phi derivative

        Gamma_approx_level = 0

        tip_force_global = [0;0;0]                                         %[N] globally applied tip force
        tip_force_local = [0;0;0]                                          %[N] locally applied tip force
        tip_moment_global = [0;0;0]                                        %[N] globally applied tip moment
        tip_moment_local  = [0;0;0]                                        %[N] locally applied tip moment

        harmonic_load_global = [];

        distributed_force_global = [0;0;0]                                 %[N] globally applied tip force
        distributed_force_local = [0;0;0]                                  %[N] locally applied tip force
        distributed_moment_global = [0;0;0]                                %[N] globally applied tip moment
        distributed_moment_local = [0;0;0]                                 %[N] locally applied tip moment


        reflect_aeroPanels = false                                         %[true/false] reflect aerodynamic panels in fuselage plane
        isAero
        AICs                                                               %[-] aerodynamic influence coefficients

        symmetry_plane_normal = [0;1;0];                                   %[-] normal vector of fuselage symmetry plane
        symmetry_point_in_plane = [0;0;0];                                 %[-] a point that lies in the fuselage symmetry plane

        %--                                                                %[-] useful rotational transforms between Euler, Wing and Aircraft reference systems
        R_W_WE = eye(3)
                        R_W_WE_tr = eye(3)
        R_A_W = eye(3)
                        R_A_W_tr
        TD = eye(3)
                        TD_tr

        idx_ribs                                                           %spanwise index locations of ribs, used only for plotting
        wingRoot_offset_A = [0;0;0]                                        %[m] wing root offset from rRoot_A to the [W] system origin

                        FLAG_shear                                         %[true/false] true if shear states included in problem

        CrossSectionProfiles = 'box'                                       %[-] cross sectional geometries of airfoil sections;
        zRot_Edraw_cell = {0}
        temp_properties

    end

    properties (SetAccess = private)

       %StiffnessMatrix
                        StiffnessMatrix_kappaHalf
                        StiffnessMatrix_shearHalf
       %DampingMatrix
                        DampingMatrix_kappaHalf
                        DampingMatrix_shearHalf

                        Pvec_appliedGlobal_G                               %[N] the set of all external global loads applied at the spanwise locations in 's'
                        Pvec_appliedLocal_I                                %[N] the set of all external local loads applied at the spanwise locations in 's'
                        Mvec_appliedGlobal_G                               %[N] the set of all external global moments applied at the spanwise locations in 's'
                        Mvec_appliedLocal_I                                %[N] the set of all external local moments applied at the spanwise locations in 's'

       %s                                                                  %[m] spanwise evaluations points along the wing
                        del_s                                              %[m] distance between neighbouring s points
                        ns                                                 %[-] number of s points
                        L                                                  %[m] total length of the wing

       %s_aero                                                             %[m] spanwise locations of aero panel boundaries
                        del_s_aero
                        s_aeroMid                                          %[m] spanwise location of aero mid points
                        nsAp                                               %[-] number of spanwise aero panels

    end

    properties
    %   ASD_CELL
                        alpha_root = 0                                     %[deg] structural angle of attack at the root (useful for wind tunnel modelling)
                        sweep_root = 0                                     %[deg] wing root sweep
                        dihedral_root = 0                                  %[deg] wing root dihedral
    end

    properties %Shape function property set
     %...properties..... .........................
      %..................derived_properties.......
       %|---------------|-------------------------
        shape_BCs_bend                                                     %[-] bending shape function boundary conditions
        shape_class_bend                                                   %[string] chosen bending shape function set e.g. 'chebyshev_1st'
                        shapeObject_bend
        shape_BCs_twist                                                    %[-] twist shape function boundary conditions
        shape_class_twist                                                  %[string] chosen twist shape function set
                        shapeObject_twist
        shape_BCs_shear = ones(3,2)                                        %[-] shear shape function boundary conditions
        shape_class_shear = 'chebyshev_1st'                                                 %[string] chosen shear shape function set
                        shapeObject_shear
        shape_BCs_exten = ones(3,2)                                        %[-] extension shape function boundary conditions
        shape_class_exten = 'chebyshev_1st'                                                 %[string] chosen extension shape function set
                        shapeObject_exten
        B_custom_bend                                                      %[-] optional custom bending shape functions to include in full shape set
        dB_custom_bend                                                     %[-] derivative of custom bending shape functions
        s_custom_bend                                                      %[-] s values at which custom bending functions are specified
        B_custom_twist                                                     %[-] optional custom twist shape functions to include in full shape set
        dB_custom_twist                                                    %[-] derivative of custom twist shape functions
        s_custom_twist                                                     %[-] s values at which custom twist functions are specified
        B_custom_shear                                                     %[-] optional custom shear shape functions to include in full shape set
        dB_custom_shear                                                    %[-] derivative of custom shear shape functions
        s_custom_shear                                                     %[-] s values at which custom shear functions are specified
        B_custom_exten                                                     %[-] optional custom extension shape functions to include in full shape set
        dB_custom_exten                                                    %[-] derivative of custom extension shape functions
        s_custom_exten                                                     %[-] s values at which custom extension functions are specified

        nmo                                                                %[-] [No._theta_functions, No._psi_functions, No._phi_functions, No._tau_z_functions, No._tau_x_functions, No._tau_y_functions]

                        %B, dB                                                  %[-] matrix defining each shape functions used
                        %----------temporary subdivisions of B
                        B_tr, dB_tr %temp
                        B_th, B_th_tr, dB_th,% B_th_tr%temp?
                        B_si, B_si_tr, dB_si,% B_si_tr%temp?
                        B_ph, B_ph_tr, dB_ph,% B_ph_tr%temp?
                        B_Sx, B_Sx_tr, dB_Sx
                        B_Sy, B_Sy_tr, dB_Sy
                        B_Sz, B_Sz_tr, dB_Sz

                        Ba, Ba_tr, dBa, dBa_tr

%                        BB, BB_tr %temp
%                        BdB
%                        dBB
                       %onsB, onsdB
        qth,qsi,qph,qSx,qSy,qSz,qAero,nQAero
        nqs, nqa, nq2nd
        th_idx, si_idx, ph_idx, Sx_idx, Sy_idx, Sz_idx
        dth_idx, dsi_idx, dph_idx, dSx_idx, dSy_idx, dSz_idx
    end

    properties (Dependent)
        t
    end

    %#ok<*MCSUP>

%==========================================================================
%//////////////////////////////////////////////////////////////////////////
%==========================================================================
    methods %constuctor and parameter methods

    function obj = NBS_flexPart_nonlinear(Master_Object,name,varargin)
        obj.Parent = utility_functions.get_option(varargin,'Parent',[]);
        obj.partName = name;
        obj.NBS_Master = Master_Object;
        obj.NBS_Master.flexParts_nonlinear.(obj.partName) = obj;
        obj.NBS_Master.flexParts_nonlinear_cell{end+1} = obj;
        obj.NBS_Master.allParts_cell{end+1} = obj;
        obj.NBS_Master.allParts_struct.(obj.partName) = obj;

        obj.qSx.n = 0; obj.qSy.n = 0; obj.qSz.n = 0;
        obj.qth.group = []; obj.qsi.group = []; obj.qph.group = [];
        obj.qSx.group = []; obj.qSy.group = []; obj.qSz.group = [];
    end

    function update_R_A_W(obj,varargin) %TODO check use of this method

        SDA_CELL = utility_functions.get_option(varargin,'sweep_dihedral_alpha_deg',{});
        if isempty(SDA_CELL)
            varargin = [varargin , 'R_A_W' , obj.R_A_W];
        end

        obj.R_A_W = obj.return_R_A_W(varargin{:},'reflect',obj.reflectedPart);

    end

    function populate_shape_set(obj,varargin)

        %what does function do?
        %primary function: generating shape function properties of
        %SimObject; e.g. B, dBB...
        %if shape objects are not already defined, generates these and
        %  writes them to the SimObject
        %if custom functions exist in object then add these to the front
        %  the the shape sets and orthogonalise the entire set


        if nargin == 0
            disp(' ')
            disp('method call of the form ''SimObject = SimObject.populate_shape_set(options)''')
            disp('...')
            disp(' ')
            disp('additional keyword/value pairs specified in ''options'' argument')
            disp(' ')
            disp('optional keywords (default values marked with ''*''):')
            disp('        halfShape: [false*/''left''/''right''] - if halfShape not false then return left or right half sets as specified')
            disp('        PLOT:      [false*/true] - will plot bending and twist shape function sets if true')
            disp('        setName:   [string] - assign a custom name to the generated shape objects')
            disp('        ...')
            return
        end

    %options:
        %halfShape
        %PLOT
        %setName
        %orthNorm
        %objectBend
        %objectTwist

    %required object fields:
        %O.nmo, O.s
        %O.shape_class_bend, O.shape_class_twist
        %O.shape_BCs_bend, O.shape_BCs_twist

    %special functionality if pre-existing properties
        %O.shapeObject_bend, O.shapeObject_twist
        %O.B_custom_bend, O.B_custom_twist

        halfShape = utility_functions.get_option(varargin,'halfShape',false); %options: 'false' / 'left' / 'right'
        PLOT = utility_functions.get_option(varargin,'PLOT',false);
        %setName = get_option(varargin,'setName',[]);
        setName = obj.partName;

        nqth = obj.qth.n; nqsi = obj.qsi.n; nqph = obj.qph.n;
        nqTx = obj.qSx.n; nqTy = obj.qSy.n; nqTz = obj.qSz.n;

        B_bend=[];B_twist=[];B_shear=[];B_exten=[];
        if isempty(obj.shapeObject_bend)
            obj.shapeObject_bend  = create_shape_object(obj.shape_class_bend ,obj.s,obj.B_custom_bend ,obj.dB_custom_bend ,'nShapes',max(nqth,nqsi)  ,'s',obj.s,'BCs',obj.shape_BCs_bend  ,'halfShape',halfShape,'setName',[setName ' (bend)']);
            y_bend  = obj.shapeObject_bend.y;  dy_ds_bend  = obj.shapeObject_bend.dy_ds;
            B_bend  = [reshape(y_bend,size(y_bend,1),1,[])   reshape(dy_ds_bend,size(dy_ds_bend,1),1,[])];
        end
        if isempty(obj.shapeObject_twist)
            obj.shapeObject_twist = create_shape_object(obj.shape_class_twist,obj.s,obj.B_custom_twist,obj.dB_custom_twist,'nShapes',nqph         ,'s',obj.s,'BCs',obj.shape_BCs_twist ,'halfShape',halfShape,'setName',[setName ' (twist)']);
            y_twist = obj.shapeObject_twist.y; dy_ds_twist = obj.shapeObject_twist.dy_ds;
            B_twist  = [reshape(y_twist,size(y_twist,1),1,[])   reshape(dy_ds_twist,size(dy_ds_twist,1),1,[])];
            if numel(B_twist)==0, B_twist = zeros(0,2,obj.ns); end
        end
        if isempty(obj.shapeObject_shear)% && (Sn~=0 || Sm~=0)
            obj.shapeObject_shear = create_shape_object(obj.shape_class_shear,obj.s,obj.B_custom_shear,obj.dB_custom_shear,'nShapes',max(nqTz,nqTx),'s',obj.s,'BCs',obj.shape_BCs_shear ,'halfShape',halfShape,'setName',[setName ' (shear)']);
            y_shear = obj.shapeObject_shear.y; dy_ds_shear = obj.shapeObject_shear.dy_ds;
            B_shear  = [reshape(y_shear,size(y_shear,1),1,[])   reshape(dy_ds_shear,size(dy_ds_shear,1),1,[])];
        end
        if isempty(obj.shapeObject_exten)% &&  So~=0
            obj.shapeObject_exten = create_shape_object(obj.shape_class_exten,obj.s,obj.B_custom_exten,obj.dB_custom_exten,'nShapes',nqTy        ,'s',obj.s,'BCs',obj.shape_BCs_exten ,'halfShape',halfShape,'setName',[setName ' (exten)']);
            y_exten = obj.shapeObject_exten.y; dy_ds_exten = obj.shapeObject_exten.dy_ds;
            B_exten  = [reshape(y_exten,size(y_exten,1),1,[])   reshape(dy_ds_exten,size(dy_ds_exten,1),1,[])];
        end

        obj.B_th = B_bend(1:nqth,1,:);  obj.B_si = B_bend(1:nqsi,1,:);   obj.B_ph = B_twist(1:nqph,1,:);
        obj.B_th_tr = permute(obj.B_th,[2 1 3]); obj.B_si_tr = permute(obj.B_si,[2 1 3]); obj.B_ph_tr = permute(obj.B_ph,[2 1 3]);
        obj.B_Sz = B_shear(1:nqTz,1,:); obj.B_Sx = B_shear(1:nqTx,1,:); obj.B_Sy = B_exten(1:nqTy,1,:);
        obj.B_Sz_tr = permute(obj.B_Sz,[2 1 3]); obj.B_Sx_tr = permute(obj.B_Sx,[2 1 3]); obj.B_Sy_tr = permute(obj.B_Sy,[2 1 3]);

        obj.dB_th = B_bend(1:nqth,2,:);  obj.dB_si = B_bend(1:nqsi,2,:);   obj.dB_ph = B_twist(1:nqph,2,:);
        obj.dB_Sz = B_shear(1:nqTz,2,:); obj.dB_Sx = B_shear(1:nqTx,2,:); obj.dB_Sy = B_exten(1:nqTy,2,:);

%        O.B_th_tr = permute(O.B_th,[2 1 3]); O.B_si_tr = permute(O.B_si,[2 1 3]); O.B_ph_tr = permute(O.B_ph,[2 1 3]);

        obj.Ba = [obj.B_th;obj.B_si;obj.B_ph]; obj.dBa = [obj.dB_th;obj.dB_si;obj.dB_ph];
        obj.Ba_tr = permute(obj.Ba,[2 1 3]); obj.dBa_tr = permute(obj.dBa,[2 1 3]);

%        O.temp.B = [O.B_th;O.B_si;O.B_ph;O.B_Sz;O.B_Sx;O.B_Sy]; O.dB = [O.dB_th;O.dB_si;O.dB_ph;O.dB_Sz;O.dB_Sx;O.dB_Sy];
%        O.temp.B_tr = permute(O.B,[2 1 3]); O.dB_tr = permute(O.dB,[2 1 3]);
        obj.temp_properties.BB = reshape(obj.Ba.*obj.Ba_tr,(nqth+nqsi+nqph)^2,1,obj.ns); obj.temp_properties.BB_tr = reshape(obj.temp_properties.BB,1,(nqth+nqsi+nqph)^2,obj.ns);
%        O.temp.onsB = reshape(bsxfun(@times,O.B*0+1,O.B_tr),(n+m+o)^2,1,O.ns);
        obj.temp_properties.BdB = reshape(obj.Ba.*obj.dBa_tr,(nqth+nqsi+nqph)^2,1,obj.ns);
% %        O.onsdB = reshape(bsxfun(@times,O.B*0+1,O.dB_tr),(n+m+o)^2,1,O.ns);
%         O.dBB = reshape(bsxfun(@times,O.dB,O.B_tr),(n+m+o)^2,1,O.ns);


        yidx2 = [zeros(1,nqth)+1 zeros(1,nqsi)+2 zeros(1,nqph)+3];
        yidx3 = [...
            repmat(yidx2,1,nqth)+3*0,...
            repmat(yidx2,1,nqsi)+3*1,...
            repmat(yidx2,1,nqph)+3*2];

        obj.temp_properties.yidx2 = yidx2;
        obj.temp_properties.yidx3 = yidx3;


        yidx22 = [zeros(1,nqth)+1 zeros(1,nqsi)+2 zeros(1,nqph)+3 zeros(1,nqTz)+4 zeros(1,nqTx)+5 zeros(1,nqTy)+6];
        yidx33 = [...
            repmat(yidx22,1,nqth)+6*0,...
            repmat(yidx22,1,nqsi)+6*1,...
            repmat(yidx22,1,nqph)+6*2,...
            repmat(yidx22,1,nqTz)+6*3,...
            repmat(yidx22,1,nqTx)+6*4,...
            repmat(yidx22,1,nqTy)+6*5];

        obj.temp_properties.yidx22 = yidx22;
        obj.temp_properties.yidx33 = yidx33;


        if isempty(obj.qth.group); obj.qth.group = ['qth_' obj.partName]; end
        if isempty(obj.qsi.group); obj.qsi.group = ['qsi_' obj.partName]; end 
        if isempty(obj.qph.group); obj.qph.group = ['qph_' obj.partName]; end
        if isempty(obj.qSx.group); obj.qSx.group = ['qSx_' obj.partName]; end 
        if isempty(obj.qSy.group); obj.qSy.group = ['qSy_' obj.partName]; end
        if isempty(obj.qSz.group); obj.qSz.group = ['qSz_' obj.partName]; end

        
        if PLOT
            obj.shapeObject_bend.plotShapes(obj.shapeObject_bend,'output_detail','final');
            obj.shapeObject_twist.plotShapes(obj.shapeObject_twist,'output_detail','final');
            if obj.shapeObject_shear.nShapes ~= 0
                obj.shapeObject_shear.plotShapes(obj.shapeObject_shear,'output_detail','final');
            end
            if obj.shapeObject_exten.nShapes ~= 0
                obj.shapeObject_exten.plotShapes(obj.shapeObject_exten,'output_detail','final');
            end
        end




        function shapeObject = create_shape_object(shapeTemplate,s_custom,B_custom,dB_custom,varargin)
            shapeObject  = ShapeFunctionObject(shapeTemplate,varargin{:});
            if ~isempty(B_custom)
                shapeObject = shapeObject.addCustomFunctions(s_custom,B_custom,dB_custom);
                if ~isempty(shapeObject.weightFunction)
                    shapeObject = shapeObject.orthNorm();
                    shapeObject = shapeObject.orth();
                end
            end
        end

    end

    function draw_part(obj,varargin)

        Tidx = utility_functions.get_option(varargin,'Tidx',obj.NBS_Master.nt);

        request_qoi_write = utility_functions.get_option(varargin,'qoiRequest',true);

        zRot_Edraw_cell_ = obj.zRot_Edraw_cell;
        assert(isa(zRot_Edraw_cell_,'cell'),'argument ''zRot_Edraw_cell'' must be a cell array');
        if numel(zRot_Edraw_cell_) == 2
            zRot_Edraw = cat(3,reshape(zRot_Edraw_cell_{1},1,1,[]),...
                zeros(1,1,obj.ns-numel(zRot_Edraw_cell_{1})-numel(zRot_Edraw_cell_{2})),...
                reshape(zRot_Edraw_cell_{2},1,1,[]));
        elseif numel(zRot_Edraw_cell_) == 1
            zRot_Edraw = zRot_Edraw_cell_{1};
        end

        zRmat_Edraw = utility_functions.r_matrix([0;0;1],reshape(zRot_Edraw,1,1,[]));

        if isempty(obj.c)
            width = obj.w./cos(reshape(zRot_Edraw,1,1,[]));
        else
            width = obj.c./cos(reshape(zRot_Edraw,1,1,[]));
        end
        height = obj.h+obj.s*0;

        if request_qoi_write

        %specify a qoi Request for the required plotting fields
        qoiRequest = {...
            true,'ex_G','1:ns',num2str(Tidx);...
            true,'ey_G','1:ns',num2str(Tidx);...
            true,'ez_G','1:ns',num2str(Tidx);...
            true,'Gamma_G','1:ns',num2str(Tidx)};

        %write qoi Request to the QOI_Container
        obj.QOI_Container.qoi_request_partLevel = qoiRequest;

        %populate any extra info required based on above request
        obj.Parent.QOI_master.write_QOI_values('partLevel');

        end

        %get E triad data
        E_draw = [...
            obj.Parent.get_qoiValue(obj,'ex_G','Tidx',Tidx,'generate_QOIs',false);...
            obj.Parent.get_qoiValue(obj,'ey_G','Tidx',Tidx,'generate_QOIs',false);...
            obj.Parent.get_qoiValue(obj,'ez_G','Tidx',Tidx,'generate_QOIs',false)];

        E_draw_3D = reshape(E_draw,3,3,[]);
        E_draw_3D_rotated = utility_functions.MultiProd_(E_draw_3D,zRmat_Edraw);

        E_draw = reshape(E_draw_3D_rotated,9,1,[]);

        %get Gamma data
        Gamma_draw = obj.Parent.get_qoiValue(obj,'Gamma_G','Tidx',Tidx,'generate_QOIs',false);

        %call the plotting function for this flexPart
        varargin = [varargin,{'CrossSectionProfiles'},{obj.CrossSectionProfiles}];
        obj.NBS_Master.draw_genericPart(obj.s,E_draw,Gamma_draw,width,height,obj.beam_cntr,varargin{:});

    end

    end



    methods (Access = {?NBS_Master})

        function set_dependent_properties(obj, tbf)
            %populate additional dependent wing parameters
            obj.R_W_WE_tr = obj.R_W_WE.';
            %obj.s_aeroMid = (obj.s_aero(2:end)+obj.s_aero(1:end-1))/2;
            %obj.nsAp = length(obj.s_aero)-1;
            %if numel(O.nmo)==3, O.nmo(6) = 0; end %pad zeros for shear if required
            obj.nqa = obj.qth.n + obj.qsi.n + obj.qph.n;
            obj.nqs = obj.qSx.n + obj.qSy.n + obj.qSz.n;
            obj.nq2nd = obj.nqa + obj.nqs;
            if obj.nqs ~= 0, obj.FLAG_shear = true; else, obj.FLAG_shear = false; end

            mass_resampling(obj);

            transfm_nodes2aero = utility_functions.sampleMat(obj.s,obj.s_aero);
            transfm_nodes2aeroMid = utility_functions.sampleMat(obj.s,obj.s_aeroMid);
            obj.Apn_idx = (1:obj.ns)*transfm_nodes2aero;
            obj.Apm_idx = (1:obj.ns)*transfm_nodes2aeroMid;

            assert(isa(obj.isAero,'logical'),'Property ''isAero'' must be set true or false');
            if obj.isAero, obj.NBS_Master.aeroPartNames{end+1} = obj.partName; end

            %--------------------
            if ~isempty(obj.varTheta0_W)
                obj.E0_W = rTransform_thetaVec_to_Rmat(obj.varTheta0_W);
            end

            %calculate E0_A
            %use E0_A to get th/si/ph 0 d/ds0
            %also calculate R_A_W from E0_A

            if ~isempty(obj.E0_W)
                E0_W_ = obj.E0_W;
            else
                nszrs(1,1,obj.ns) = 0;
                st0 = sin(obj.th0)+nszrs; ss0 = sin(obj.si0)+nszrs; sp0 = sin(obj.ph0)+nszrs;
                ct0 = cos(obj.th0)+nszrs; cs0 = cos(obj.si0)+nszrs; cp0 = cos(obj.ph0)+nszrs;

                ey_WE = [ct0.*ss0 ; ct0.*cs0 ; st0];
                ex_WE = [cs0.*cp0 + st0.*ss0.*sp0;- ss0.*cp0 + st0.*cs0.*sp0;-ct0.*sp0];
                ez_WE = [cs0.*sp0 - st0.*ss0.*cp0;- ss0.*sp0 - st0.*cs0.*cp0; ct0.*cp0];

                E0_WE = utility_functions.MultiProd_([ex_WE ey_WE ez_WE],obj.TD);
                E0_W_ = utility_functions.MultiProd_(obj.R_W_WE,utility_functions.MultiProd_(E0_WE,obj.R_W_WE_tr));
            end
            %------------------------------------------------------------------
            %calculate aerodynamic panel reference frame
            %if EAp_A is already provided, use this information
            %E0_A = MultiProd_(O.R_A_W,E0_W);
            E0_aero_ref_W = E0_W_;
            if ~isempty(obj.EAp_W)
                obj.EAp_I = utility_functions.MultiProd_(permute(E0_aero_ref_W,[2 1 3]),obj.EAp_W);
                %if EAp_A is not provided, calculate aerodynamic reference frame
                %via projection onto the wing surface.
            else
                aero_ref_direction_A = [1;0;0];
                aero_ref_direction_W = obj.R_A_W.'*aero_ref_direction_A;
                aero_ref_direction_I = utility_functions.mult_Anmz_Bmp1(permute(E0_aero_ref_W,[2 1 3]),aero_ref_direction_W); %calculate the [1;0;0] vector in the intrinsic system
                aero_ref_projection_I = [1;0;0].*aero_ref_direction_I; %projection of the intrinsic [1;0;0] vector onto the wing plane

                %note: EAp aligns with free-stream direction
                %regardless of which side of aircraft being considered
                exAp_pr_I = aero_ref_projection_I.*(1./sum(aero_ref_projection_I.^2).^0.5); %exAp_pr_I = normalised projection of the intrinsic [1;0;0] vector onto the wing plane
                eyAp_pr_I = utility_functions.MultiProd_([0 -1 0;1 0 0;0 0 0],exAp_pr_I); %eyAp_pr_I = cross( [0;0;1] , exAp_pr_I )
                ezAp_pr_I = repmat([0;0;1],[1 1 size(exAp_pr_I,3)]); %ezAp_pr_I = [0;0;1]

                EAp_pr_I = [exAp_pr_I eyAp_pr_I ezAp_pr_I];

                %now rotate EAp_pr_I projected frame to achieve required alpha_I or alpha_A distribution
                if isempty(obj.alpha_A) && isempty(obj.alpha_I)
                    obj.alpha_I = zeros(1,1,obj.ns);
                end
                if ~isempty(obj.alpha_A)
                    alpha_aero_ref_I = asin(aero_ref_direction_I(3,1,:)); %alpha value of the projected system EAp_pr before rotation
                    obj.EAp_I = utility_functions.MultiProd_(utility_functions.r_matrix(eyAp_pr_I,obj.alpha_A*pi/180-alpha_aero_ref_I),EAp_pr_I);
                    obj.alpha_I = [];
                elseif ~isempty(obj.alpha_I)
                    obj.EAp_I = utility_functions.MultiProd_(utility_functions.r_matrix(eyAp_pr_I,obj.alpha_I*pi/180),EAp_pr_I)*0 + eye(3);% [-1;-1;1].*eye(3);
                    obj.alpha_A = [];
                end

                obj.TD_tr = permute(obj.TD,[2 1 3]);

            end

            %------------------------------------------------------------------
            obj.c_pm = utility_functions.sample(obj.c,obj.Apm_idx,3);
            obj.c_pn = utility_functions.sample(obj.c,obj.Apn_idx,3);
            obj.beam_cntr_pm = utility_functions.sample(obj.beam_cntr,obj.Apm_idx,3);
            obj.aero_cntr_pm = utility_functions.sample(obj.aero_cntr,obj.Apm_idx,3);
            obj.beam_cntr_pn = utility_functions.sample(obj.beam_cntr,obj.Apn_idx,3);
            obj.del_s_aero = diff(obj.s_aero);
            obj.EAp_I_pm = utility_functions.sample(obj.EAp_I,obj.Apm_idx,3);
            obj.EAp_I_pn = utility_functions.sample(obj.EAp_I,obj.Apn_idx,3);
            obj.ApWidths_pm = abs(obj.del_s_aero.*obj.EAp_I_pm(2,2,:));

            %--------------------
            Rs_W_WE_tr = [1;1;1];
            %simplification of R_W_WE for 90 degree rotations; Rs_W_WE_tr = sign of each column entry, Rv_W_WE_tr = index of each column entry
            for i_ = 1:3
                [rR,~] = find(round(obj.R_W_WE_tr)); Rv_W_WE_tr = rR;
                [~,cRn] = find(min(obj.R_W_WE_tr,0)); Rs_W_WE_tr(cRn) = -1;
            end
            obj.temp_properties.Rs_W_WE_tr = Rs_W_WE_tr; obj.temp_properties.Rv_W_WE_tr = Rv_W_WE_tr;
            %--------------------
            obj.temp_properties.R_W_WE_iseye = isequal(obj.R_W_WE,eye(3));
            obj.temp_properties.R_A_WE = obj.R_A_W.*obj.R_W_WE;
            %         if obj.root_idx ~= 1
            %             obj.StiffnessMatrix(:,:,obj.root_idx) = zeros(3);
            %             obj.DampingMatrix(:,:,obj.root_idx) = zeros(3);
            %             obj.I_varTheta(:,:,obj.root_idx) = zeros(3);
            %             obj.ms(1,1,obj.root_idx) = 0;
            %         end

            %         szSM = size(obj.StiffnessMatrix);
            %         if szSM(1)==3 && szSM(2)==3
            %             obj.StiffnessMatrix_kappaHalf = [obj.StiffnessMatrix obj.StiffnessMatrix*0];
            %             obj.StiffnessMatrix_shearHalf = obj.StiffnessMatrix_kappaHalf*0;
            %         elseif szSM(1)==6 && szSM(2)==6
            %             obj.StiffnessMatrix_kappaHalf = obj.StiffnessMatrix(1:3,:,:);
            %             obj.StiffnessMatrix_shearHalf = obj.StiffnessMatrix(4:end,:,:);
            %             obj.DampingMatrix_kappaHalf = obj.DampingMatrix(1:3,:,:);
            %             obj.DampingMatrix_shearHalf = obj.DampingMatrix(4:end,:,:);
            %         else
            %             error('StiffnessMatrix property must have dimension 3x3 or 6x6');
            %         end


                if strcmp(obj.NBS_Master.aerodynamics,'strip_unsteady')
                    obj.qAero.n = obj.nsAp*2;
                    obj.qAero.group = ['AeroStates_' obj.partName];
                end


                if strcmp(obj.NBS_Master.aerodynamics,'WT')
                    
                    if tbf.FLAG_unsteady
                        obj.qAero.n = 0;
                        obj.qAero.group = ['AeroStates_' obj.partName];
                        
                        if strcmp(tbf.aeroForces, 'leishman')
                            obj.qAero.n = obj.nsAp*2;
                            obj.qAero.group = ['AeroStates_' obj.partName];
                        end
                        if strcmp(tbf.aeroForces, 'dynamic stall')
                            if strcmp(tbf.dsModel, 'oye')
                                obj.nQAero = 1;
                                obj.qAero.n = obj.nsAp*obj.nQAero;
                                obj.qAero.group = ['AeroStates_' obj.partName];
                            elseif strcmp(tbf.dsModel, 'larsen')
                                obj.nQAero = 4;
                                obj.qAero.n = obj.nsAp*obj.nQAero;
                                obj.qAero.group = ['AeroStates_' obj.partName];
                            end
                        end
                        if tbf.FLAG_dw
                            if strcmp(tbf.aeroForces, 'lookup 2D')
                                obj.nQAero = 0;
                                switch tbf.dwDetail
                                    case 'full'
                                        obj.qAero.n = obj.nsAp*4;
                                    case 'simple'
                                        obj.qAero.n = 1;
                                end
                            else                              
                                switch tbf.dwDetail
                                    case 'full'
                                        obj.qAero.n = obj.qAero.n + obj.nsAp*4;
                                    case 'simple'
                                        obj.qAero.n = obj.qAero.n + 1;
                                end
                            end
                        end
                    end
                end
             
            if ~isfield(obj.aeroData, 'aeroCoeff2D')
                obj.aeroData.aeroCoeff2D = [];
            end
            
            function mass_resampling(O)
                %function uses O.mps, O.msDiscrete, O.I_varTheta_ps_I and O.I_varTheta_discrete_I
                %to return the aggregated quantities O.ms, O.I_varTheta and O.massOffset_I

                if isempty(O.ms)
                    del_snode = cat(3,O.del_s(1),O.del_s(2:end)+O.del_s(1:end-1),O.del_s(end))/2;

                    O.msContinuous = O.mps.*del_snode;
                    O.ms = O.msContinuous + O.msDiscrete;
                end

                if isempty(O.massOffset_I) || sum(O.massOffset_I(:))==0
                    O.massOffset_I = zeros(3,1,O.ns);
                    O.massOffsetFlag = false;
                else
                    O.massOffsetFlag = true;
                end

                if isempty(O.I_varTheta)
                    O.I_varTheta_continuous_I = bsxfun(@times,O.I_varTheta_ps_I,del_snode);
                    O.I_varTheta = O.I_varTheta_continuous_I + O.I_varTheta_discrete_I;
                end
            end

        end

    end

%==========================================================================
%//////////////////////////////////////////////////////////////////////////
%==========================================================================

    methods

        function t = get.t(obj)
            t = obj.Parent.t;
        end

        function set.alpha_root(obj,val)
            %if alpha_root is a single value then
            %set the parameter and also update the root rotation matrix

            %%if alpha_root is a 1x1 cell containing a single value then
            %update the parameter alone
            if isa(val,'cell')
                obj.alpha_root = val{1};
            else
                obj.alpha_root = val;
                Rinfo = utility_functions.rTransform_projection_to_Rmat(obj.sweep_root,obj.dihedral_root,'alphaXZdihedral',obj.alpha_root);
                obj.R_A_W = {Rinfo.rotationMatrix};
            end
        end

        function set.sweep_root(obj,val)
            %see alpha_root set method
            if isa(val,'cell')
                obj.sweep_root = val{1};
            else
                obj.sweep_root = val;
                Rinfo = utility_functions.rTransform_projection_to_Rmat(obj.sweep_root,obj.dihedral_root,'alphaXZdihedral',obj.alpha_root);
                obj.R_A_W = {Rinfo.rotationMatrix};
            end
        end

        function set.dihedral_root(obj,val)
            %see alpha_root set method
            if isa(val,'cell')
                obj.dihedral_root = val{1};
            else
                obj.dihedral_root = val;
                Rinfo = utility_functions.rTransform_projection_to_Rmat(obj.sweep_root,obj.dihedral_root,'alphaXZdihedral',obj.alpha_root);
                obj.R_A_W = {Rinfo.rotationMatrix};
            end
        end

        function set.R_A_W(obj,Rmat)
            %set method for R_A_W

            %if R_A_W is a 3x3 matrix, then in addition the alpha_deg,
            %sweep_deg and dihedral_deg properties as well as the
            %R_A_W root rotation matrix will be updated

            %if R_A_W is a 3x3 matrix contained within a 1x1 cell
            %array then only R_A_W will be updated (special
            %functionality only used within this class)

            %blank entries in ASD_CELL will be ignored

            if isa(Rmat,'cell')
                assert(isequal(size(Rmat{1}),[3,3]) , 'property R_A_W must be a 3by3 rotation matrix');
                obj.R_A_W = Rmat{1};
            else
                assert(isequal(size(Rmat),[3,3]) , 'property R_A_W must be a 3by3 rotation matrix');
                obj.R_A_W = Rmat;
                Rinfo = utility_functions.rTransform_Rmat_to_projection(Rmat);
                obj.alpha_root = {Rinfo.alphaXZdihedral_deg};
                obj.sweep_root = {Rinfo.sweep_deg};
                obj.dihedral_root = {Rinfo.dihedral_deg};
%                 obj.ASD_CELL = {{alpha_deg,sweep_deg,dihed_deg}};
            end
        end

        function set.s(obj,val)
            assert(issorted(val),'Elements of s must consist of an ascending vector of spanwise locations')
            obj.s = reshape(val(:),1,1,[]);
            obj.ns = length(obj.s);
            obj.L = obj.s(end) - obj.s(1);
            obj.del_s = diff(obj.s);

            if numel(obj.distributed_force_global)==3, obj.distributed_force_global = zeros(3,1,obj.ns); end
            if numel(obj.distributed_force_local)==3, obj.distributed_force_local = zeros(3,1,obj.ns); end
            if numel(obj.distributed_moment_global)==3, obj.distributed_moment_global = zeros(3,1,obj.ns); end
            if numel(obj.distributed_moment_local)==3, obj.distributed_moment_local = zeros(3,1,obj.ns); end
        end

        function set.s_aero(obj,val)
            assert(issorted(val),'Elements of s_aero must consist of an ascending vector of spanwise locations')
            obj.s_aero = reshape(val(:),1,1,[]);
            obj.nsAp = length(obj.s_aero)-1;
            obj.s_aeroMid = (obj.s_aero(2:end)+obj.s_aero(1:end-1))/2;
            obj.del_s_aero = diff(obj.s_aero);
        end

        function set.StiffnessMatrix(obj,val)
            obj.StiffnessMatrix = val;
            szSM = size(val);
            if szSM(1)==3 && szSM(2)==3
                obj.StiffnessMatrix_kappaHalf = [obj.StiffnessMatrix obj.StiffnessMatrix*0];
                obj.StiffnessMatrix_shearHalf = obj.StiffnessMatrix_kappaHalf*0;
            elseif szSM(1)==6 && szSM(2)==6
                obj.StiffnessMatrix_kappaHalf = obj.StiffnessMatrix(1:3,:,:);
                obj.StiffnessMatrix_shearHalf = obj.StiffnessMatrix(4:end,:,:);
            else
                error('StiffnessMatrix property must have dimension 3x3 or 6x6');
            end
        end

        function set.DampingMatrix(obj,val)
            obj.DampingMatrix = val;
            szSM = size(val);
            if szSM(1)==3 && szSM(2)==3
                obj.DampingMatrix_kappaHalf = [obj.DampingMatrix obj.DampingMatrix*0];
                obj.DampingMatrix_shearHalf = obj.DampingMatrix_kappaHalf*0;
            elseif szSM(1)==6 && szSM(2)==6
                obj.DampingMatrix_kappaHalf = obj.DampingMatrix(1:3,:,:);
                obj.DampingMatrix_shearHalf = obj.DampingMatrix(4:end,:,:);
            else
                error('DampingMatrix property must have dimension 3x3 or 6x6');
            end
        end

        function set.I_varTheta(obj,val)
            %Avoid numerical issues if sectional inertia is torsional only
            if val(1,1,end)==0, val(1,1) = norm(val)*1e-6; end
            if val(3,3,end)==0, val(3,3) = norm(val)*1e-6; end

            obj.I_varTheta = val;
        end

        function set.tip_force_global(obj,val)
            obj.tip_force_global = reshape(val,3,1,1);

            obj.Pvec_appliedGlobal_G = obj.distributed_force_global;
            obj.Pvec_appliedGlobal_G(:,1,end) = obj.Pvec_appliedGlobal_G(:,1,end) + obj.tip_force_global;
        end

        function set.harmonic_load_global(obj,val)

                obj.harmonic_load_global = val;

         end

        function set.distributed_force_global(obj,val)
            obj.distributed_force_global = reshape(val,3,1,[]);

            obj.Pvec_appliedGlobal_G = obj.distributed_force_global;
            obj.Pvec_appliedGlobal_G(:,1,end) = obj.Pvec_appliedGlobal_G(:,1,end) + obj.tip_force_global;
        end

        function set.tip_force_local(obj,val)
            obj.tip_force_local = reshape(val,3,1,1);

            obj.Pvec_appliedLocal_I = obj.distributed_force_local;
            obj.Pvec_appliedLocal_I(:,1,end) = obj.Pvec_appliedLocal_I(:,1,end) + obj.tip_force_local;
        end

        function set.distributed_force_local(obj,val)
            obj.distributed_force_local = reshape(val,3,1,[]);

            obj.Pvec_appliedLocal_I = obj.distributed_force_local;
            obj.Pvec_appliedLocal_I(:,1,end) = obj.Pvec_appliedLocal_I(:,1,end) + obj.tip_force_local;
        end

        function set.tip_moment_global(obj,val)
            obj.tip_moment_global = reshape(val,3,1,1);

            obj.Mvec_appliedGlobal_G = obj.distributed_moment_global;
            obj.Mvec_appliedGlobal_G(:,1,end) = obj.Mvec_appliedGlobal_G(:,1,end) + obj.tip_moment_global;
        end

        function set.distributed_moment_global(obj,val)
            obj.distributed_moment_global = reshape(val,3,1,[]);

            obj.Mvec_appliedGlobal_G = obj.distributed_moment_global;
            obj.Mvec_appliedGlobal_G(:,1,end) = obj.Mvec_appliedGlobal_G(:,1,end) + obj.tip_moment_global;
        end

        function set.tip_moment_local(obj,val)
            obj.tip_moment_local = reshape(val,3,1,1);

            obj.Mvec_appliedLocal_I = obj.distributed_moment_local;
            obj.Mvec_appliedLocal_I(:,1,end) = obj.Mvec_appliedLocal_I(:,1,end) + obj.tip_moment_local;
        end

        function set.distributed_moment_local(obj,val)
            obj.distributed_moment_local = reshape(val,3,1,[]);

            obj.Mvec_appliedLocal_I = obj.distributed_moment_local;
            obj.Mvec_appliedLocal_I(:,1,end) = obj.Mvec_appliedLocal_I(:,1,end) + obj.tip_moment_local;
        end

        function set.Pvec_appliedGlobal_G(obj,val)
            if isempty(val), val = [0;0;0]; end
            obj.Pvec_appliedGlobal_G = val;
        end

        function set.Pvec_appliedLocal_I(obj,val)
            if isempty(val), val = [0;0;0]; end
            obj.Pvec_appliedLocal_I = val;
        end

        function set.Mvec_appliedGlobal_G(obj,val)
            if isempty(val), val = [0;0;0]; end
            obj.Mvec_appliedGlobal_G = val;
        end

        function set.Mvec_appliedLocal_I(obj,val)
            if isempty(val), val = [0;0;0]; end
            obj.Mvec_appliedLocal_I = val;
        end

    end

%==========================================================================
%//////////////////////////////////////////////////////////////////////////
%==========================================================================
    methods (Static) %Utility functions performing useful operations

    function R_A_W = return_R_A_W(varargin)

        reflect_centre = utility_functions.get_option(varargin,'reflect',false);
        SDA_CELL = utility_functions.get_option(varargin,'sweep_dihedral_alpha_deg',{0,0,0});
        [sweep_deg,dihedral_deg,alpha_deg] = deal(SDA_CELL{:});
        R_A_W = utility_functions.get_option(varargin,'R_A_W',[]);

        if isempty(R_A_W)
            Rinfo = rTransform_projection_to_Rmat(sweep_deg,dihedral_deg,'alphaXZdihedral',alpha_deg);
            R_A_W = Rinfo.rotationMatrix;
        end

        if reflect_centre
            R_A_W = diag([1,-1,1])*R_A_W*diag([-1,1,1]);                 %reflect R_A_W in XZ plane and reverse orientation of ex to maintain a right-handed system
        end

    end

    function Asample = sample(A,eval_idx,dim)
        %linearly sample matrix A at the indices eval_idx along the dimension 'dim'
        %example:
        %A = [1 2 3;2 4 6];
        %>> sample(A,[1.1 2.5 3],2)
        %    1.1    2.5    3.0
        %    2.2    5.0    6.0

        %---- extract the integer part of the evaluation indices
        eval_idx_floor  = floor(eval_idx);
        indices_fl  = {':',':',':'}; indices_fl{dim}  = eval_idx_floor;
        %---- modified integer part to handle special case eval_idx(i_) = size(A,dim)
        eval_idx_floor_ = min(eval_idx_floor,size(A,dim)-1);
        indices_fl_ = indices_fl;    indices_fl_{dim} = eval_idx_floor_;
        %---- extract the decimal part of the evaluation indices
        indices_dc = eval_idx - eval_idx_floor;
        dimensions_indices_dc = [1,1,1]; dimensions_indices_dc(dim) = length(eval_idx);

        indices_dc_rs = reshape(indices_dc,dimensions_indices_dc);

        Andim = length(size(A));
        Adiff = diff(A,1,dim);

        %linear sampling of A matrix
        Asample = A(indices_fl{1:Andim}) + Adiff(indices_fl_{1:Andim}).*indices_dc_rs;
    end

    function [MOMENT_xi , FORCE_xi] = material_law( KAPPA_I , ~ , KAPPA_0_I , TAU , ~ , TAU_0 , Linear_Stiffness_Matrix , ~)
        if isempty(TAU)
            xi = KAPPA_I - KAPPA_0_I;
            MOMENT_xi = -utility_functions.MultiProd_(Linear_Stiffness_Matrix{1}(1:3,1:3,:),xi);
%             cubicTerms = 5e-1.*utility_functions.MultiProd_(Linear_Stiffness_Matrix{1}(1:3,1:3,:),xi.^3);
%             quinticTerms = -7e-1*utility_functions.MultiProd_(Linear_Stiffness_Matrix{1}(1:3,1:3,:),xi.^5);
%             MOMENT_xi = MOMENT_xi + cubicTerms; %softening  
%             MOMENT_xi = MOMENT_xi - cubicTerms; %hardening
            FORCE_xi = zeros(3,0,0);
        else
            xi = [KAPPA_I - KAPPA_0_I ; TAU - TAU_0];
            MOMENT_xi = -utility_functions.MultiProd_(Linear_Stiffness_Matrix{1},xi);
            FORCE_xi  = -utility_functions.MultiProd_(Linear_Stiffness_Matrix{2},xi);
        end
    end

    function [MOMENT_dxidt , FORCE_dxidt] = damping_law( ~ , dKAPPA_dt_I , ~ , ~ , dTAU_dt , ~ , Linear_Damping_Matrix , ~)
        if isempty(dTAU_dt)
            dxidt = dKAPPA_dt_I;
            MOMENT_dxidt = -utility_functions.MultiProd_(Linear_Damping_Matrix{1}(1:3,1:3,:),dxidt);
            FORCE_dxidt  = zeros(3,0,0);
        else
            dxidt = [dKAPPA_dt_I ; dTAU_dt];
            MOMENT_dxidt = -utility_functions.MultiProd_(Linear_Damping_Matrix{1},dxidt);
            FORCE_dxidt  = -utility_functions.MultiProd_(Linear_Damping_Matrix{2},dxidt);
        end
    end

    end






    methods

        function [dW_dq_part,dM_dq_part,partInformationStruct] = f_flexPart_nonlinear(obj,tbf,Q,partInformationStruct,outputFormat,tidx,i_flex_part)
            %#ok<*PROPLC>          
            
            SimObject = obj.NBS_Master;
            flex_part_name = obj.partName;

            %--------------------------------------------------------------
            parentObj = obj.Parent;
            parentName = parentObj.partName;
            parentConnIdx = obj.connection_idx_ParentObj;
            Omega_G = partInformationStruct.(parentName).dvarTheta_dt_G;
%             beta_ = partInformationStruct.(parentName).azimuth;
            R_G_A = partInformationStruct.(parentName).E_G(:,:,parentConnIdx);
            dR_G_A_dt = partInformationStruct.(parentName).dE_dt_G(:,:,parentConnIdx);
            dR_G_A_dqe_Dim3x3x1xnqe = partInformationStruct.(parentName).dE_dq_G(:,:,parentConnIdx,:);
            d2R_G_A_dt2_star = partInformationStruct.(parentName).d2E_dt2_G_star(:,:,parentConnIdx);
%             dOmega_dt_G_star = partInformationStruct.(parentName).d2varTheta_dt2_G_star;
            wingRoot_offset_G = R_G_A*obj.wingRoot_offset_A;

            R_W_WE = obj.R_W_WE;
            R_W_WE_tr = obj.R_W_WE_tr;
            R_A_W = obj.R_A_W;
            TD = obj.TD;

            %%% specification of rotational transforms
            R_G_W = utility_functions.mult_Anm1_Bmpz(R_G_A,R_A_W);
            if isempty(dR_G_A_dqe_Dim3x3x1xnqe), dR_G_W_dqr_Dim3x3x1xnqr = zeros(3,3,1,0);
            else, dR_G_W_dqr_Dim3x3x1xnqr = utility_functions.MultiProd_(dR_G_A_dqe_Dim3x3x1xnqe,R_A_W); end
            R_A_G = R_G_A.';

            Rs_W_WE_tr = [1;1;1];
            for i_ = 1:3, [rR,~] = find(round(R_W_WE_tr)); Rv_W_WE_tr = rR; [~,cRn] = find(min(R_W_WE_tr,0)); Rs_W_WE_tr(cRn) = -1; end
            %--------------------------------------------------------------

            dqe_idx = partInformationStruct.(parentName).dq2nd_idx;

            th_idx = obj.th_idx.';
            si_idx = obj.si_idx.';
            ph_idx = obj.ph_idx.';
            Sx_idx = obj.Sx_idx.';
            Sy_idx = obj.Sy_idx.';
            Sz_idx = obj.Sz_idx.';

            dth_idx = obj.dth_idx.';
            dsi_idx = obj.dsi_idx.';
            dph_idx = obj.dph_idx.';
            dSx_idx = obj.dSx_idx.';
            dSy_idx = obj.dSy_idx.';
            dSz_idx = obj.dSz_idx.';
           
            nqs = obj.nqs;
            nqa = obj.nqa;
            nqr = size(dR_G_A_dqe_Dim3x3x1xnqe,4);

            dq2nd_idx = [...
                dth_idx(:);dsi_idx(:);dph_idx(:);dSx_idx(:);dSy_idx(:);dSz_idx(:);dqe_idx(:)];

            q2nd_idx = [...
                th_idx(:);si_idx(:);ph_idx(:);Sx_idx(:);Sy_idx(:);Sz_idx(:)];
%         
%             beta_ = beta_(1);
%             az_ib = 2*pi/SimObject.tbf.nB*(i_flex_part-1);
%             p1 = [1 cos(beta_ + az_ib) sin(beta_ + az_ib)];
%             p2 = [0 -sin(beta_ + az_ib) cos(beta_ + az_ib)];
% 
%                 q2nd_idx_b = q2nd_idx;
%                 dq2nd_idx_b = dq2nd_idx;
%                 p1_ = repmat(p1, [length(q2nd_idx_b) 1]);
%                 p2_ = repmat(p2, [length(dq2nd_idx_b) 1]);
%             
%             nQ = length(Q);
%             P = zeros(nQ(1),3,'like',Q);
%             P(q2nd_idx,:) = p1_;
%             P(dq2nd_idx,:) = Omega_G(1)*p2_;
%             
%             for i = 1:nQ
%             Q(i) = sum(P(i,:)\Q(i),1);
%             end
  
            [bools,dqg2nd_to_dq2nd_idx] = ismember(dq2nd_idx,SimObject.dqg2nd_idx); % dqg2nd_idx(dqg_to_dq_idx) = dq2nd_idx
            assert(all(bools),'Invalid Mapping');

            FLAG_shear = obj.FLAG_shear;

            s = obj.s;
            nszrs = s*0;
            ns = obj.ns;
            del_s = obj.del_s;
%             L = obj.L;
            ms = obj.ms;
            I_varTheta = obj.I_varTheta;
            KAPPA_0_I = obj.KAPPA_0_I;

            %manipulations of shape function vectors
            Ba_tr = obj.Ba_tr;
            B_th = obj.B_th; B_si = obj.B_si; B_ph = obj.B_ph;
            B_th_tr = obj.B_th_tr; B_si_tr = obj.B_si_tr; B_ph_tr = obj.B_ph_tr;
            dB_th = obj.dB_th; dB_si = obj.dB_si; dB_ph = obj.dB_ph;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Define Attitude Parmaeters + Spatial and Temporal derivatives
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %~~ Euler angles referencing intrinsic system
            th = utility_functions.mult_Anm1_Bmpz(Q(th_idx).',B_th)+obj.th0; %spanwise theta distribution      (1)x(1)x(ns)
            si = utility_functions.mult_Anm1_Bmpz(Q(si_idx).',B_si)+obj.si0; %spanwise psi distribution        (1)x(1)x(ns)
            ph = utility_functions.mult_Anm1_Bmpz(Q(ph_idx).',B_ph)+obj.ph0; %spanwise phi distribution        (1)x(1)x(ns)

            %~~ Angular temporal derivatives
            dth_dt = utility_functions.mult_Anm1_Bmpz(Q(dth_idx).',B_th); %                            (1)x(1)x(ns)
            dsi_dt = utility_functions.mult_Anm1_Bmpz(Q(dsi_idx).',B_si); %                            (1)x(1)x(ns)
            dph_dt = utility_functions.mult_Anm1_Bmpz(Q(dph_idx).',B_ph); %                            (1)x(1)x(ns)
            %~~ Angular spacial derivatives
            dth_ds = utility_functions.mult_Anm1_Bmpz(Q(th_idx).',dB_th)+obj.dth_ds0; %                  (1)x(1)x(ns)
            dsi_ds = utility_functions.mult_Anm1_Bmpz(Q(si_idx).',dB_si)+obj.dsi_ds0; %                  (1)x(1)x(ns)
            dph_ds = utility_functions.mult_Anm1_Bmpz(Q(ph_idx).',dB_ph)+obj.dph_ds0; %                  (1)x(1)x(ns)
            %~~ Angular dsdt derivatives
            d2th_dsdt = utility_functions.mult_Anm1_Bmpz(Q(dth_idx).',dB_th); %                        (1)x(1)x(ns)
            d2si_dsdt = utility_functions.mult_Anm1_Bmpz(Q(dsi_idx).',dB_si); %                        (1)x(1)x(ns)
            d2ph_dsdt = utility_functions.mult_Anm1_Bmpz(Q(dph_idx).',dB_ph); %                        (1)x(1)x(ns)
            %~~~~~~~~~~~~~~~~~~~~~
            if ~FLAG_shear
                [tau_x, tau_y, tau_z]= deal([]);
                [dtau_x_dt,dtau_y_dt,dtau_z_dt] = deal([]);
                [dtau_x_dqsx, dtau_y_dqsy, dtau_z_dqsz] = deal([]);
                [dtau_x_dq, dtau_y_dq, dtau_z_dq] = deal([]);
            else
                B_Sx = obj.B_Sx;   B_Sy = obj.B_Sy;   B_Sz = obj.B_Sz;
                B_Sx_tr = obj.B_Sx_tr;   B_Sy_tr = obj.B_Sy_tr;   B_Sz_tr = obj.B_Sz_tr;
                
                %~~ Shear terms
                tau_x = utility_functions.mult_Anm1_Bmpz(Q(Sx_idx).',B_Sx);
                tau_y = utility_functions.mult_Anm1_Bmpz(Q(Sy_idx).',B_Sy);
                tau_z = utility_functions.mult_Anm1_Bmpz(Q(Sz_idx).',B_Sz);
                %~~ Shear temporal derivatives
                dtau_x_dt = utility_functions.mult_Anm1_Bmpz(Q(dSx_idx).',B_Sx);
                dtau_y_dt = utility_functions.mult_Anm1_Bmpz(Q(dSy_idx).',B_Sy);
                dtau_z_dt = utility_functions.mult_Anm1_Bmpz(Q(dSz_idx).',B_Sz);
                %~~ Shear state derivatives
                dtau_x_dqsx = permute(B_Sx_tr,[1 4 3 2]);    % tau_x variation with respect to shear_x states
                dtau_y_dqsy = permute(B_Sy_tr,[1 4 3 2]);    % tau_y variation with respect to shear_y states
                dtau_z_dqsz = permute(B_Sz_tr,[1 4 3 2]);    % tau_z variation with respect to shear_z states
                nSx = numel(Sx_idx);
                nSy = numel(Sy_idx);
                nSz = numel(Sz_idx);
                dtau_x_dq = cat(4,zeros(1,1,ns,nqa),dtau_x_dqsx,zeros(1,1,ns,nSy+nSz),zeros(1,1,ns,nqr));    % tau_x variation with respect to all states
                dtau_y_dq = cat(4,zeros(1,1,ns,nqa),zeros(1,1,ns,nSx),dtau_y_dqsy,zeros(1,1,ns,nSz),zeros(1,1,ns,nqr));    % tau_y variation with respect to all states
                dtau_z_dq = cat(4,zeros(1,1,ns,nqa),zeros(1,1,ns,nSx+nSy),dtau_z_dqsz,zeros(1,1,ns,nqr));    % tau_z variation with respect to all states
            end
            
            %~~ Evaluate trigonometric terms
            st = sin(th); ss = sin(si); sp = sin(ph); %                                (1)x(1)x(ns)
            ct = cos(th); cs = cos(si); cp = cos(ph); %                                (1)x(1)x(ns)

            st_ss = st.*ss;
            st_cs = st.*cs;
            st_sp = st.*sp;
            st_cp = st.*cp;
            ct_ss = ct.*ss;
            ct_cs = ct.*cs;
            ct_sp = ct.*sp;
            ct_cp = ct.*cp;
            ss_sp = ss.*sp;
            ss_cp = ss.*cp;
            cs_sp = cs.*sp;
            cs_cp = cs.*cp;

            st_ss_sp = st.*ss.*sp;
            st_ss_cp = st.*ss.*cp;
            st_cs_sp = st.*cs.*sp;
            st_cs_cp = st.*cs.*cp;
            ct_ss_sp = ct.*ss.*sp;
            ct_ss_cp = ct.*ss.*cp;
            ct_cs_sp = ct.*cs.*sp;
            ct_cs_cp = ct.*cs.*cp;
            
            %list of all compound terms
            %st_ss,st_cs,st_sp,st_cp,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp
            %st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Specification of Intrinsic Orthonormal System
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dzeta_a_dt = [dth_dt;dsi_dt;dph_dt];
            dzeta_a_dt_tr = permute(dzeta_a_dt,[2 1 3]);%[dth_dt,dsi_dt,dph_dt];

            dR_G_W_dqe_Dim3x3x1xnqe = dR_G_W_dqr_Dim3x3x1xnqr;

            [E_WE,E_W,E_G,dE_dqa_G_Dim3x3xnsxnqa,dE_dqe_G_Dim3x3xnsxnqe,dE_dt_W] = static_method_groups.E_methods.get_E_group1(...
                ...
                st, ct,...
                st_ss,st_cs,st_sp,st_cp,...
                ct_ss,ct_cs,ct_sp,ct_cp,...
                ss_sp,ss_cp,cs_sp,cs_cp,...
                st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,...
                ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp,...
                ...
                Rs_W_WE_tr,Rv_W_WE_tr,TD,nszrs,ns,...
                R_A_W,R_G_W,...
                ...
                dzeta_a_dt_tr,Ba_tr,obj.temp_properties.yidx2,...
                ...
                dR_G_W_dqe_Dim3x3x1xnqe);

            ey_WE = E_WE(:,2,:);
            E_G_tr = permute(E_G,[2 1 3]);

            dE_dq_G_Dim3x3xnsxnq2nd = cat(4,dE_dqa_G_Dim3x3xnsxnqa,zeros(3,3,ns,nqs),dE_dqe_G_Dim3x3xnsxnqe);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Calculation of varTheta rotation vectors
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dvarTheta_dqe_G_Dim3x1x1xnqe = partInformationStruct.(parentName).dvarTheta_dq_G;

            dvarTheta_dt_root_G = partInformationStruct.(parentName).dvarTheta_dt_G;
            dvarTheta_dt_root_G_skew = partInformationStruct.(parentName).OmegaSkew_G;
            d2varTheta_dt2_root_G_star = partInformationStruct.(parentName).d2varTheta_dt2_G_star;

            [dvarTheta_dqa_G_Dim3x1xnsxnqa,dvarTheta_dt_G,d2varTheta_dt2_G_star] = static_method_groups.varTheta_methods.get_varTheta_group1(...
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
                d2varTheta_dt2_root_G_star);

            dvarTheta_dq_G_Dim3x1xnsxnq2nd = cat(4,dvarTheta_dqa_G_Dim3x1xnsxnqa , zeros(3,1,ns,nqs) , repmat(dvarTheta_dqe_G_Dim3x1x1xnqe,1,1,ns));
%             dvarTheta_dq_G_Dim1x3xnsxnq2nd = permute(dvarTheta_dq_G_Dim3x1xnsxnq2nd,[2 1 3 4]);
            dvarTheta_dq_G_Dim3xnq2ndxns = permute(dvarTheta_dq_G_Dim3x1xnsxnq2nd,[1 4 3 2]);
            dvarTheta_dq_G_Dimnq2ndx3xns = permute(dvarTheta_dq_G_Dim3x1xnsxnq2nd,[4 1 3 2]);

            dvarTheta_dt_G_skew = utility_functions.getSkewMat(dvarTheta_dt_G);
            d2varTheta_dt2_G_star_skew = utility_functions.getSkewMat(d2varTheta_dt2_G_star);

            dE_dt_G = utility_functions.MultiProd_(dvarTheta_dt_G_skew,E_G);
            d2E_dt2_G_star = utility_functions.MultiProd_(d2varTheta_dt2_G_star_skew,E_G) + utility_functions.MultiProd_(dvarTheta_dt_G_skew,dE_dt_G);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Calculation of Gamma position vectors
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rBarA_G = partInformationStruct.(parentName).Gamma_G;
            drBarA_dt_G = partInformationStruct.(parentName).dGamma_dt_G;
            d2rBarA_dt2_G_star = partInformationStruct.(parentName).d2Gamma_dt2_G_star;
            drBarA_G_dqr_G_Dim3x1x1xnqr = partInformationStruct.(parentName).dGamma_dq_G;

            Gamma_approx_lvl = obj.Gamma_approx_level;

            massOffset_I = obj.massOffset_I;
            FLAG_massOffset = obj.massOffsetFlag;

            Gamma_root_G = rBarA_G + wingRoot_offset_G;
            dGamma_dt_root_G = drBarA_dt_G + dR_G_A_dt*obj.wingRoot_offset_A;
            d2Gamma_dt2_root_G_star = d2rBarA_dt2_G_star + d2R_G_A_dt2_star*obj.wingRoot_offset_A;

            dGamma_dq_root_G_Dim3x1x1xnq2nd = cat(4,zeros(3,1,1,nqa+nqs) , drBarA_G_dqr_G_Dim3x1x1xnqr + utility_functions.MultiProd_(dR_G_A_dqe_Dim3x3x1xnqe,obj.wingRoot_offset_A));

            [Gamma_G,dGamma_dt_G,dGamma_dq_G_Dim3x1xnsxnq2nd,dGamma_m_dq_G_Dim3x1xnsxnq2nd,d2Gamma_dt2_G_star,d2Gamma_m_dt2_G_star] = static_method_groups.Gamma_methods.get_Gamma_group1(...
                SimObject.Gamma_int_fnc,del_s,ns,obj.root_idx,...
                E_W,dE_dt_G,d2E_dt2_G_star,E_G,dE_dq_G_Dim3x3xnsxnq2nd,...
                tau_x,tau_y,tau_z,...
                dtau_x_dt,dtau_y_dt,dtau_z_dt,...
                dtau_x_dq,dtau_y_dq,dtau_z_dq,...
                massOffset_I,FLAG_massOffset,...                               %dR_G_W_dqr456_931,FLAG_free_free,...
                FLAG_shear,Gamma_approx_lvl,...
                Gamma_root_G,dGamma_dt_root_G,dGamma_dq_root_G_Dim3x1x1xnq2nd,d2Gamma_dt2_root_G_star);

%             dGamma_m_dq_G_tr_Dim1x3xnsxnq2nd = permute(dGamma_m_dq_G_Dim3x1xnsxnq2nd,[2 1 3 4]);
            dGamma_m_dq_G_Dim3xnq2ndxns = permute(dGamma_m_dq_G_Dim3x1xnsxnq2nd,[1 4 3 2]);
            dGamma_m_dq_G_Dimnq2ndx3xns = permute(dGamma_m_dq_G_Dim3x1xnsxnq2nd,[4 1 3 2]);

            %%blah +- consider how d2Gamma_m_dt2_G_star and d2Gamma_dt2_G_star affect the folding wing tip definition

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Calculation of Kappa curvature vectors
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [KAPPA_I,dKAPPA_dt_I,dKAPPA_dqa_I_tr_Dim1x3xnsxnqa] = static_method_groups.kappa_methods.get_kappa_group1(...
                ...
                st,ct,sp,cp,dth_ds,dsi_ds,dph_ds,nszrs,R_W_WE,B_th,B_si,B_ph,dB_th,dB_si,dB_ph,st_sp,st_cp,ct_sp,ct_cp,obj.temp_properties.yidx2,TD,...
                dzeta_a_dt_tr,d2th_dsdt,d2si_dsdt,d2ph_dsdt);

            if ~FLAG_shear
                dTAU_dqs_tr = [];
                TAU = [];
                dTAU_dt = [];
            else
                TAU = [tau_x;tau_y;tau_z];
                dTAU_dt = [dtau_x_dt;dtau_y_dt;dtau_z_dt];
                dTAU_dqs = cat(4, dtau_x_dqsx.*[1;0;0] , dtau_y_dqsy.*[0;1;0] , dtau_z_dqsz.*[0;0;1]);
                dTAU_dqs_tr = permute(dTAU_dqs,[2 1 3 4]);       
            end
            
            TAU_0 = [0;0;0];
            Linear_Stiffness_Matrix = {obj.StiffnessMatrix_kappaHalf;obj.StiffnessMatrix_shearHalf};
            Linear_Damping_Matrix = {obj.DampingMatrix_kappaHalf;obj.DampingMatrix_shearHalf};

            Higher_Order_Arguments = [];

            [MOMENT_xi , FORCE_xi]       = obj.material_law( KAPPA_I , dKAPPA_dt_I , KAPPA_0_I , TAU , dTAU_dt , TAU_0 , Linear_Stiffness_Matrix , Higher_Order_Arguments);

            [MOMENT_dxidt , FORCE_dxidt] = obj.damping_law( KAPPA_I , dKAPPA_dt_I , KAPPA_0_I , TAU , dTAU_dt , TAU_0 , Linear_Damping_Matrix   , Higher_Order_Arguments);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% FlexPart_nonlinear Applied Loads from PvecApplied_G, MvecApplied_G, Gravitational Acceleration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isprop(tbf, 'operation')
                PvecWeight_G = SimObject.grav_acc.*SimObject.gravVec_G.*obj.ms;
            else
                switch tbf.operation
                    case 'fixed'
                        if tbf.FLAG_parked
                            PvecWeight_G = SimObject.grav_acc.*SimObject.gravVec_G.*obj.ms;
                        else
                            Omega_G = repmat([0; 0; 2*pi/SimObject.T], [1 1 obj.ns]);
                            PvecWeight_G = - utility_functions.crossn(2*obj.ms.*Omega_G, dGamma_dt_G, 1) - utility_functions.crossn(obj.ms.*Omega_G,utility_functions.crossn(Omega_G, Gamma_G, 1),1);
                        end
                    case {'forced', 'free'}
                        PvecWeight_G = SimObject.grav_acc.*SimObject.gravVec_G.*obj.ms;
                end
            end
            
            massOffset_G = utility_functions.MultiProd_(E_G,massOffset_I);
            
            PvecApplied_G = obj.Pvec_appliedGlobal_G + utility_functions.MultiProd_(E_G,obj.Pvec_appliedLocal_I) + PvecWeight_G;
            MvecApplied_G = obj.Mvec_appliedGlobal_G + utility_functions.MultiProd_(E_G,obj.Mvec_appliedLocal_I) + utility_functions.MultiProd_(utility_functions.getSkewMat(massOffset_G),PvecWeight_G);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% FlexPart_nonlinear Aerodynamic Loads PvecApplied_G, MvecApplied_G
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %PvecAero_pm_G: aerodynamic forces acting at quarter chord of each panel

            nsAp = obj.nsAp;

            if obj.isAero == true && isa(SimObject.aerodynamics,'char')

                Apm_idx = obj.Apm_idx;
                Apn_idx = obj.Apn_idx;

%                 if isempty(SimObject.aeroPartNames)
%                     aeroPartNames = [aeroPartNames {flex_part_name}];
%                 end

                switch SimObject.aerodynamics
                    case {'strip_steady','strip_unsteady','WT'}
                    E_G_pm = utility_functions.sample(E_G,Apm_idx,3);
                    dE_dt_G_pm = utility_functions.sample(dE_dt_G,Apm_idx,3);
                    d2E_dt2_G_star_pm = utility_functions.sample(d2E_dt2_G_star,Apm_idx,3);
                    % dGamma_dt_G_pm = sample(dGamma_dt_G,Apm_idx,3);
                    dvarTheta_dt_G_pm = utility_functions.sample(dvarTheta_dt_G,Apm_idx,3);
                    EAp_I_pm = obj.EAp_I_pm;
                    EAp_G_pm = utility_functions.MultiProd_(E_G_pm,EAp_I_pm);
                    dEAp_dt_G_pm = utility_functions.MultiProd_(dE_dt_G_pm,EAp_I_pm);
                    d2EAp_dt2_G_star_pm = utility_functions.MultiProd_(d2E_dt2_G_star_pm,EAp_I_pm);
                    %dexAp_dt_G_pm = dEAp_dt_G_pm(:,1,:);

                    %     dexAp_dt_G_pm = dEAp_dt_G_pm(:,1,:);
                    %     alphaCP = 3/4;
                    %     dGammaAlphaCP_dt_G_pm = bsxfun(@plus,drBarA_dt_G,...
                    %                            dGamma_dt_G_pm + bsxfun(@times,c_pm.*(alphaCP - beam_cntr_pm),dexAp_dt_G_pm));
                    partInformationStruct.(flex_part_name).aero_cntr_pm = obj.aero_cntr_pm;
                    partInformationStruct.(flex_part_name).nsAp = nsAp;
                    partInformationStruct.(flex_part_name).EAp_G_pm = EAp_G_pm;
                    partInformationStruct.(flex_part_name).dEAp_dt_G_pm = dEAp_dt_G_pm;
                    partInformationStruct.(flex_part_name).d2EAp_dt2_G_star_pm = d2EAp_dt2_G_star_pm;
                    partInformationStruct.(flex_part_name).Gamma_G_pm = utility_functions.sample(Gamma_G,Apm_idx,3);
                    partInformationStruct.(flex_part_name).d2Gamma_dt2_G_star_pm = utility_functions.sample(d2Gamma_dt2_G_star,Apm_idx,3);;
                    partInformationStruct.(flex_part_name).Gamma_G_pn = utility_functions.sample(Gamma_G,Apn_idx,3);
                    partInformationStruct.(flex_part_name).dGamma_dt_G_pm = utility_functions.sample(dGamma_dt_G,Apm_idx,3);
                    %partInformationStruct.(flex_part_name).dGammaAlphaCP_dt_G_pm = dGammaAlphaCP_dt_G_pm;
                    partInformationStruct.(flex_part_name).dvarTheta_dt_G_pm = dvarTheta_dt_G_pm;
                    partInformationStruct.(flex_part_name).chord = obj.c_pm;
                    partInformationStruct.(flex_part_name).ApWidth = obj.ApWidths_pm;
                    partInformationStruct.(flex_part_name).beam_cntr_pm = obj.beam_cntr_pm;
                    partInformationStruct.(flex_part_name).AIC = obj.AICs;
                    partInformationStruct.(flex_part_name).dGamma_dq_G_pm = utility_functions.sample(dGamma_dq_G_Dim3x1xnsxnq2nd,Apm_idx,3);
                    partInformationStruct.(flex_part_name).dvarTheta_dq_G_pm = utility_functions.sample(dvarTheta_dq_G_Dim3x1xnsxnq2nd,Apm_idx,3);
                    partInformationStruct.(flex_part_name).aeroData = obj.aeroData;
%                          F_SUM = utility_functions.sample(PvecApplied_G, Apm_idx,3);
%                          M_SUM = utility_functions.sample(MvecApplied_G, Apm_idx,3);
            
%                     partInformationStruct.(flex_part_name).F_SUM = F_SUM;
%                     partInformationStruct.(flex_part_name).M_SUM = M_SUM;

                    case {'VLM_steady'}

                    Gamma_G_pn = utility_functions.sample(Gamma_G,Apn_idx,3);
                    Gamma_G_pm = utility_functions.sample(Gamma_G,Apm_idx,3);
                    dGamma_dt_G_pm = utility_functions.sample(dGamma_dt_G,Apm_idx,3);
                    aero_cntr_pn = 0.25;
                    control_point_pm = 0.75;
                    beam_cntr_pm = obj.beam_cntr_pm;
                    beam_cntr_pn = obj.beam_cntr_pn;
                    chord_pm = obj.c_pm;
                    chord_pn = obj.c_pn;
                    E_G_pm = utility_functions.sample(E_G,Apm_idx,3);
                    E_G_pn = utility_functions.sample(E_G,Apn_idx,3);
                    EAp_I_pm = obj.EAp_I_pm;
                    EAp_I_pn = obj.EAp_I_pn;
                    EAp_G_pm = utility_functions.MultiProd_(E_G_pm,EAp_I_pm);
                    EAp_G_pn = utility_functions.MultiProd_(E_G_pn,EAp_I_pn);
                    dE_dt_G_pm = utility_functions.sample(dE_dt_G,Apm_idx,3);
                    dEAp_dt_G_pm = utility_functions.MultiProd_(dE_dt_G_pm,EAp_I_pm);
                    Pv = Gamma_G_pn + bsxfun(@times, (aero_cntr_pn - beam_cntr_pn).*chord_pn , EAp_G_pn(:,1,:) );
                    Pc = Gamma_G_pm + bsxfun(@times, (control_point_pm - beam_cntr_pm).*chord_pm , EAp_G_pm(:,1,:) );
                    PcNorm = EAp_G_pm(:,3,:);
                    dPc_dt = dGamma_dt_G_pm + bsxfun(@times, (control_point_pm - beam_cntr_pm).*chord_pm , dEAp_dt_G_pm(:,1,:) );

                    partInformationStruct.(flex_part_name).Pv = Pv;
                    partInformationStruct.(flex_part_name).Pc = Pc;
                    partInformationStruct.(flex_part_name).PcNorm = PcNorm;
                    partInformationStruct.(flex_part_name).dPc_dt = dPc_dt;
                    partInformationStruct.(flex_part_name).symmetric_plane_cell = [];
                end
            else
                    Apm_idx = obj.Apm_idx;
                    partInformationStruct.(flex_part_name).harmonic_load_global = obj.harmonic_load_global;
                    partInformationStruct.(flex_part_name).nsAp = nsAp;
                    partInformationStruct.(flex_part_name).dGamma_dq_G_pm = utility_functions.sample(dGamma_dq_G_Dim3x1xnsxnq2nd,Apm_idx,3);
                    partInformationStruct.(flex_part_name).dvarTheta_dq_G_pm = utility_functions.sample(dvarTheta_dq_G_Dim3x1xnsxnq2nd,Apm_idx,3);

            end

            %==========================================================================
            %Virtual work terms from rotational kinetic energy contributions
            [dW_dq_Kinetic_Rotation_ddqComponent,dW_dq_Kinetic_Rotation_remainder] = static_method_groups.dPi_dq_Kinetic_Rotation.get_dPi_dq_Kinetic_Rotation(...
                E_G, E_G_tr, I_varTheta,...
                dvarTheta_dq_G_Dim3xnq2ndxns, dvarTheta_dq_G_Dimnq2ndx3xns,...
                dvarTheta_dt_G,d2varTheta_dt2_G_star);
            %==========================================================================
            %Virtual work terms from translational kinetic energy contributions
            [dW_dq_Kinetic_Translation_ddqComponent,dW_dq_Kinetic_Translation_remainder] = static_method_groups.dPi_dq_Kinetic_Translation.get_dPi_dq_Kinetic_Translation(...
                ms, dGamma_m_dq_G_Dim3xnq2ndxns, dGamma_m_dq_G_Dimnq2ndx3xns,...
                d2Gamma_m_dt2_G_star);
            %==========================================================================
            %Virtual work terms from material deformation
            dPi_dq_KAPPA = utility_functions.MultiProd_(dKAPPA_dqa_I_tr_Dim1x3xnsxnqa , MOMENT_xi + MOMENT_dxidt);
            dPi_dq_SHEAR = utility_functions.MultiProd_(dTAU_dqs_tr     ,  FORCE_xi +  FORCE_dxidt);
            if isempty(dPi_dq_SHEAR)
                dW_dq_xi = cat(4, SimObject.intVal(s, dPi_dq_KAPPA , SimObject.int_fnc) , zeros(1,1,1,nqr));
            else
                dW_dq_xi = cat(4, SimObject.intVal(s,cat(4, dPi_dq_KAPPA , dPi_dq_SHEAR), SimObject.int_fnc) , zeros(1,1,1,nqr));
            end
            %==========================================================================
            %Virtual work terms from applied loads
            dPi_dq_Fapplied = utility_functions.dotn(PvecApplied_G,dGamma_dq_G_Dim3x1xnsxnq2nd,1);
            dPi_dq_Mapplied = utility_functions.dotn(MvecApplied_G,dvarTheta_dq_G_Dim3x1xnsxnq2nd,1);
            dW_dq_AppliedLoad = sum(dPi_dq_Fapplied + dPi_dq_Mapplied,3);
            %==========================================================================

            %add together the virtual work contributions specific to the current flexPart_nonlinear
            dW_dq_part = ...
                - dW_dq_Kinetic_Rotation_remainder(:)...
                - dW_dq_Kinetic_Translation_remainder(:)...
                + dW_dq_xi(:)...
                + dW_dq_AppliedLoad(:);
                
            dM_dq_part = ...
                + dW_dq_Kinetic_Rotation_ddqComponent...
                + dW_dq_Kinetic_Translation_ddqComponent;
            
            partInformationStruct.(flex_part_name).E_G = E_G;
            partInformationStruct.(flex_part_name).R_A_G = R_A_G;
            partInformationStruct.(flex_part_name).R_A_W = obj.R_A_W;
            partInformationStruct.(flex_part_name).dE_dt_G = dE_dt_G;
            partInformationStruct.(flex_part_name).d2E_dt2_G_star = d2E_dt2_G_star;
            partInformationStruct.(flex_part_name).dE_dq_G = dE_dq_G_Dim3x3xnsxnq2nd;
            partInformationStruct.(flex_part_name).dvarTheta_dq_G = dvarTheta_dq_G_Dim3x1xnsxnq2nd;
            partInformationStruct.(flex_part_name).dvarTheta_dt_G = dvarTheta_dt_G;
            partInformationStruct.(flex_part_name).Gamma_G = Gamma_G;
            partInformationStruct.(flex_part_name).dGamma_dt_G = dGamma_dt_G;
            partInformationStruct.(flex_part_name).d2Gamma_dt2_G_star = d2Gamma_dt2_G_star;
            partInformationStruct.(flex_part_name).dGamma_dq_G = dGamma_dq_G_Dim3x1xnsxnq2nd;
            partInformationStruct.(flex_part_name).q2nd_idx = q2nd_idx;
            partInformationStruct.(flex_part_name).dq2nd_idx = dq2nd_idx;
            partInformationStruct.(flex_part_name).dqg2nd_to_dq2nd_idx = dqg2nd_to_dq2nd_idx;
            partInformationStruct.(flex_part_name).ns = ns;
            partInformationStruct.(flex_part_name).d2varTheta_dt2_G_star = d2varTheta_dt2_G_star;
            if isempty(obj.qAero)
                partInformationStruct.(flex_part_name).qAero_idx = [];
            else
                partInformationStruct.(flex_part_name).qAero_idx = SimObject.StateInfo{'Index',obj.qAero.group}{:};
            end

            if isequal(outputFormat,'qoi')

                QOI_Container = utility_functions.get_field(SimObject,['QOI_Master.QOIcontainers_struct.flexParts_nonlinear.' flex_part_name]);

                %Gamma_G = bsxfun(@plus,rBarA_G + R_G_A*wingRoot_offset_A,MultiProd_(R_G_W,Gamma_W));
                QOI_Container.add_qoi('Gamma_G',tidx,Gamma_G,'1:ns','\Gamma#_{[G]}','m');
%                 QOI_Container.add_qoi('MOMENT_xi',tidx,MOMENT_xi,'1:ns','\MOMENT_xi','N');
                
                Gamma_A = rBarA_G + utility_functions.MultiProd_(R_A_G,Gamma_G);
                QOI_Container.add_qoi('Gamma_A',tidx,Gamma_A,'1:ns','\Gamma#_{[A]}','m');
                Gamma_G_pn = utility_functions.sample(Gamma_G, obj.Apn_idx, 3);
                QOI_Container.add_qoi('Gamma_G_pn',tidx,Gamma_G_pn, 1:nsAp+1,'\Gamma#_{[G]}','m');
                Gamma_m_G = Gamma_G + utility_functions.MultiProd_(E_G,massOffset_I);
                QOI_Container.add_qoi('Gamma_m_G',tidx,Gamma_m_G,'1:ns','\Gamma_m#_{[G]}','m');
                %dGamma_dt_G = squeeze(bsxfun(@plus,drBarA_dt_G,map_WtoG(dGamma_dt_W,R_G_W) + MultiProd_(dR_G_W_dt,Gamma_W)));
                QOI_Container.add_qoi('dGamma_dt_G',tidx,dGamma_dt_G,'1:ns','d\Gamma#/dt_{[G]}','m/s');
                dGamma_dt_A = utility_functions.MultiProd_(R_A_G,dGamma_dt_G);
                QOI_Container.add_qoi('dGamma_dt_A',tidx,dGamma_dt_A,'1:ns','d\Gamma#/dt_{[A]}','m/s');
                
                EAp_G = utility_functions.MultiProd_(E_G,obj.EAp_I);
                EAp_A = utility_functions.MultiProd_(R_A_G, EAp_G);
                QOI_Container.add_qoi('EAp_A',tidx,reshape(EAp_A,9,1,[]),'1:ns','EAp_{[A]}','');
                QOI_Container.add_qoi('EAp_I',tidx,reshape(obj.EAp_I,9,1,[]),'1:ns','EAp_{[A]}','');
                QOI_Container.add_qoi('ex_G',tidx,E_G(:,1,:),'1:ns','ex#_{[G]}','');
                QOI_Container.add_qoi('ey_G',tidx,E_G(:,2,:),'1:ns','ey#_{[G]}','');
                QOI_Container.add_qoi('ez_G',tidx,E_G(:,3,:),'1:ns','ez#_{[G]}','');
                QOI_Container.add_qoi('th',tidx,th*180/pi,'1:ns','\theta','deg');
                QOI_Container.add_qoi('si',tidx,si*180/pi,'1:ns','\psi','deg');
                QOI_Container.add_qoi('ph',tidx,ph*180/pi,'1:ns','\phi','deg');
                QOI_Container.add_qoi('dth_dt',tidx,dth_dt*180/pi,'1:ns','d\phi/dt','deg');
                QOI_Container.add_qoi('dsi_dt',tidx,dsi_dt*180/pi,'1:ns','d\phi/dt','deg');
                QOI_Container.add_qoi('dph_dt',tidx,dph_dt*180/pi,'1:ns','d\phi/dt','deg');
                QOI_Container.add_qoi('KAPPA_I',tidx,KAPPA_I,'1:ns','\Kappa#_{[I]}','rad/m');
                QOI_Container.add_qoi('PvecApplied_G',tidx,PvecApplied_G,'1:ns','P#_{Applied[G]}','N');
                QOI_Container.add_qoi('MvecApplied_G',tidx,MvecApplied_G,'1:ns','M#_{Applied[G]}','N');

                CoM_info_flexPart_nonlinear(1) = sum(ms);
                CoM_info_flexPart_nonlinear(2:4) = sum(ms.*Gamma_G,3)/sum(ms);
                CoM_G = CoM_info_flexPart_nonlinear(2:4);
 
               QOI_Container.add_qoi('CoM_G',tidx,CoM_G,'1:3','CoM#_{[G]}','m');

                QOI_Container.discretisationVariables.ns = obj.ns;
                QOI_Container.discretisationVariables.nt = SimObject.nt;
                QOI_Container.discretisationVariables.L = obj.L;
                QOI_Container.discretisationVariables.T = SimObject.t(end);

            end

        end

    end

end
