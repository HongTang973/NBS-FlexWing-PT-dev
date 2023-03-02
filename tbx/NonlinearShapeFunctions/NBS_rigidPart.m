classdef NBS_rigidPart < handle & matlab.mixin.Copyable

    properties
        NBS_Master
        Parent
        root_idx
        connection_idx_ParentObj
        address
        partName
        QOI_Container
    end
    
    properties (Dependent)
        t
    end
    
    properties
        Gamma_W
        connectionOffset_C = [0;0;0]
        R_C_W_0 %initial rotational offset from connection frame to rigidPart frame
    end
    
    
    properties
    
    s
        ns
        L
    s_aero
        s_aeroMid
        del_s_aero
        nAnodes
        Apn_idx
        Apm_idx
        ApWidths_pm
    sAp_idx_global
        
    EAp_I = eye(3) %TODO temporary assignment
        EAp_I_pm
    
    w
    h
    
    c
        c_pm
    aero_cntr = 0.5
    aero_cntr_pm 
    beam_cntr = 0.5
        beam_cntr_pm
    aeroData = struct();
    isAero
    AICs
    CrossSectionProfiles
    
    rMu_W %position of rigid part centre of mass, the origin of this vector lies at the root_idx
    Mu %mass of rigid part
    I_varTheta
    
    Gamma_approx_level = 0
    
    tip_force_global = [0;0;0]                                             %[N] globally applied tip force
    tip_force_local = [0;0;0]                                              %[N] locally applied tip force
    tip_moment_global = [0;0;0]                                            %[N] globally applied tip moment
    tip_moment_local = [0;0;0]                                             %[N] locally applied tip moment
                    tip_load                                               %[N] statement of all applied tip loads
                    Pvec_appliedGlobal_G = 0                               %[N] the set of all external global loads applied at the spanwise locations in 's'
                    Pvec_appliedLocal_I = 0                                %[N] the set of all external local loads applied at the spanwise locations in 's'
                    Mvec_appliedGlobal_G = 0                               %[N] the set of all external global moments applied at the spanwise locations in 's'
                    Mvec_appliedLocal_I = 0                                %[N] the set of all external local moments applied at the spanwise locations in 's'

    Rmat_Edraw
    
    end
    
    
properties
    zRot_Edraw_cell = {0}
    qAero
end

    properties (SetAccess = private)
                        nsAp                                               %[-] number of spanwise aero panels
    end

    
%==========================================================================
%//////////////////////////////////////////////////////////////////////////
%==========================================================================
    methods
    
    function obj = NBS_rigidPart(Master_Object,name,varargin)
        obj.Parent = utility_functions.get_option(varargin,'Parent',[]);
        obj.partName = name;
        obj.NBS_Master = Master_Object;
        obj.NBS_Master.rigidParts.(obj.partName) = obj;
        obj.NBS_Master.rigidParts_cell{end+1} = obj;
        obj.NBS_Master.allParts_cell{end+1} = obj;
        obj.NBS_Master.allParts_struct.(obj.partName) = obj;
    end
    
    function set_dependent_properties(obj)
        obj.ns = numel(obj.s);
        obj.L = obj.s(end) - obj.s(1);
        obj.nAnodes = numel(obj.s_aero);
        assert(isa(obj.isAero,'logical'),'property ''isAero'' must be set true or false');
        if obj.isAero, obj.NBS_Master.aeroPartNames{end+1} = obj.partName; end
        
        obj.s_aeroMid = (obj.s_aero(2:end)+obj.s_aero(1:end-1))/2;
        transfm_nodes2aero = utility_functions.sampleMat(obj.s,obj.s_aero);
        transfm_nodes2aeroMid = utility_functions.sampleMat(obj.s,obj.s_aeroMid);
        obj.Apn_idx = (1:obj.ns)*transfm_nodes2aero;
        obj.Apm_idx = (1:obj.ns)*transfm_nodes2aeroMid;
        
        obj.EAp_I_pm = repmat(obj.EAp_I,1,1,numel(obj.s_aeroMid));
        obj.c_pm = utility_functions.sample(obj.c,obj.Apm_idx,3);
        obj.beam_cntr = bsxfun(@plus, obj.beam_cntr , obj.s*0);
        obj.aero_cntr = bsxfun(@plus, obj.aero_cntr , obj.s*0);
        obj.beam_cntr_pm = utility_functions.sample(obj.beam_cntr,obj.Apm_idx,3);
        obj.aero_cntr_pm = utility_functions.sample(obj.aero_cntr,obj.Apm_idx,3);
        obj.del_s_aero = diff(obj.s_aero);
        obj.ApWidths_pm = abs(obj.del_s_aero.*obj.EAp_I_pm(2,2,:));
        
        %--------------------
        obj.tip_load = [obj.tip_force_global obj.tip_force_local obj.tip_moment_global obj.tip_moment_local];
        obj.Pvec_appliedGlobal_G = obj.Pvec_appliedGlobal_G + cat(3,zeros(3,1,obj.ns-1),obj.tip_force_global);
        obj.Pvec_appliedLocal_I = obj.Pvec_appliedLocal_I + cat(3,zeros(3,1,obj.ns-1),obj.tip_force_local);
        obj.Mvec_appliedGlobal_G = obj.Mvec_appliedGlobal_G + cat(3,zeros(3,1,obj.ns-1),obj.tip_moment_global);
        obj.Mvec_appliedLocal_I = obj.Mvec_appliedLocal_I + cat(3,zeros(3,1,obj.ns-1),obj.tip_moment_local);
        %--------------------
        if obj.isAero
        if strcmp(obj.NBS_Master.aerodynamics,'strip_unsteady')
            obj.qAero.n = (obj.nAnodes-1)*1;
            obj.qAero.group = ['AeroStates_' obj.partName];
        end
        end
    end
    
    function draw_part(obj,varargin)
        
        Tidx = utility_functions.get_option(varargin,'Tidx',obj.NBS_Master.nt);
        
        request_qoi_write = utility_functions.get_option(varargin,'qoiRequest',false);
        
        %optional intrinsic z rotation applied to E_draw triads
        %zRot_Edraw_cell = get_option(varargin,'zRot_Edraw',{zeros(1,1,obj.ns)});
        zRot_Edraw_cell_ = obj.zRot_Edraw_cell;
        assert(isa(zRot_Edraw_cell_,'cell'),'argument ''zRot_Edraw_cell'' must be a cell array');
        %if numel(zRot_Edraw_cell) == 1 then directly apply the rotations in zRot_Edraw_cell{1}
        %if numel(zRot_Edraw_cell) == 2 then directly apply zRot_Edraw_cell{1} rotations to the start of the part and zRot_Edraw_cell{2} rotations from the end
        if numel(zRot_Edraw_cell_) == 2
            zRot_Edraw = cat(3,reshape(zRot_Edraw_cell_{1},1,1,[]),...
                zeros(1,1,obj.ns-numel(zRot_Edraw_cell_{1})-numel(zRot_Edraw_cell_{2})),...
                reshape(zRot_Edraw_cell_{2},1,1,[]));
        elseif numel(zRot_Edraw_cell_) == 1
            zRot_Edraw = zRot_Edraw_cell_{1};
        end
        
        zRmat_Edraw = utility_functions.r_matrix([0;0;1],reshape(zRot_Edraw,1,1,[]));
        
        width = obj.w./cos(reshape(zRot_Edraw,1,1,[]));
        
        if request_qoi_write %<this is required for the Parent object, not NBS_rigidPart
            
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
        
        %Note on coordinate systems
        %[G] - the global coordinate system
        %[W] - the rigid part coordinate system
        %[Wdraw]
        %[C] - the connection coordinate system
        %[Cdraw]
        
        %get E triad data
        R_G_C = reshape([...
            obj.NBS_Master.get_qoiValue(obj.Parent,'ex_G','Tidx',Tidx,'Sidx',obj.connection_idx_ParentObj,'generate_QOIs',false);...
            obj.NBS_Master.get_qoiValue(obj.Parent,'ey_G','Tidx',Tidx,'Sidx',obj.connection_idx_ParentObj,'generate_QOIs',false);...
            obj.NBS_Master.get_qoiValue(obj.Parent,'ez_G','Tidx',Tidx,'Sidx',obj.connection_idx_ParentObj,'generate_QOIs',false)],3,3);
        
        R_C_W = obj.NBS_Master.get_qoiValue(obj,'R_C_W','Tidx',Tidx,'Sidx',1,'generate_QOIs',false);
        R_C_W = reshape(R_C_W,3,3);
        R_G_W = R_G_C*R_C_W;
        
        E_G = repmat(R_G_W , 1,1,obj.ns);
        
        E_Gdraw = utility_functions.MultiProd_(E_G,zRmat_Edraw);
        
%         E_draw_3D = repmat(R_G_C*obj.R_C_W , 1,1,obj.ns);
%         E_draw_3D_rotated = MultiProd_(E_draw_3D,zRmat_Edraw);
        
        E_Gdraw_flat = reshape(E_Gdraw,9,1,[]);
        
%         %get Gamma data
%         Gamma_reference = obj.NBS_Master.get_qoiValue(obj.Parent,'Gamma_G','Tidx',Tidx,'Sidx',obj.connection_idx_ParentObj);
%         Gamma_G = Gamma_reference + R_G_W*reshape(obj.Gamma_W ,3,[],1);
        Gamma_G = obj.NBS_Master.get_qoiValue(obj,'Gamma_G','Tidx',Tidx,'Sidx',':','generate_QOIs',false);
        
        %call the plotting function for this flexPart
        varargin = [varargin,...
            {'CrossSectionProfiles'},{obj.CrossSectionProfiles},...
            {'idx_ribs'},{[1 2]}];
        
        obj.NBS_Master.draw_genericPart(obj.s,E_Gdraw_flat,Gamma_G,width,obj.h,obj.beam_cntr,varargin{:});
    end
    
    end
    
    
    methods

        function set.s_aero(obj,val)
            assert(issorted(val),'Elements of s_aero must consist of an ascending vector of spanwise locations')
            obj.s_aero = reshape(val(:),1,1,[]);
            obj.nsAp = length(obj.s_aero)-1;
            obj.s_aeroMid = (obj.s_aero(2:end)+obj.s_aero(1:end-1))/2;
            obj.del_s_aero = diff(obj.s_aero);
        end
        
    end
    
%==========================================================================
%//////////////////////////////////////////////////////////////////////////
%==========================================================================
    methods (Static) %Utility functions performing useful operations
    
    
    end
    
    
    methods
        
        function t = get.t(obj)
            t = obj.NBS_Master.t;
        end
        
    end


    


    methods
        
        function [dW_dq_part,dM_dq_part,partInformationStruct] = f_rigidPart(obj,Q,partInformationStruct,outputFormat,tidx)
        %#ok<*PROPLC>
        
        SimObject = obj.NBS_Master;

        rigid_part_name = obj.partName;
        
        gravAccVec = SimObject.grav_acc*SimObject.gravVec_G;
        aerodynamics = SimObject.aerodynamics;
        dqg2nd_idx = SimObject.dqg2nd_idx;
        
        parentObj = obj.Parent;
        parentName = parentObj.partName;
        parentConnIdx = obj.connection_idx_ParentObj;
        
        
        dqe_idx = partInformationStruct.(parentName).dq2nd_idx;
        
            
            qCustom_idx  = [];
            dqCustom_idx = [];
            
            dq2nd_idx = dqe_idx;
            nq2nd = numel(dq2nd_idx);
            [bools,dqg2nd_to_dq2nd_idx] = ismember(dq2nd_idx,dqg2nd_idx);
            assert(all(bools),'Invalid Mapping');
            
            
            nqC = 0;
            R_C_W = obj.R_C_W_0;
            dR_C_W_dqC = zeros(3,3,1,0);
            depsilon_dqC_C = zeros(3,1,1,0);
            depsilon_dt_C = [0;0;0];
            dR_C_W_dt = 0;
            d2R_C_W_dt2_star = 0;
            dW_dq_epsilon = 0;
            
        
        ns = obj.ns;
        
        R_G_C = partInformationStruct.(parentName).E_G(:,:,parentConnIdx);
        dR_G_C_dt = partInformationStruct.(parentName).dE_dt_G(:,:,parentConnIdx);
        dR_G_C_dqe = partInformationStruct.(parentName).dE_dq_G(:,:,parentConnIdx,:);
        d2R_G_C_dt2_star = partInformationStruct.(parentName).d2E_dt2_G_star(:,:,parentConnIdx);
        
        R_G_W = R_G_C*R_C_W;
        dR_G_W_dt = dR_G_C_dt*R_C_W + R_G_C*dR_C_W_dt;
        dR_G_W_dqC = utility_functions.MultiProd_(R_G_C,dR_C_W_dqC);
        dR_G_W_dqe = utility_functions.MultiProd_(dR_G_C_dqe,R_C_W);
        dR_G_W_dq = cat(4,dR_G_W_dqC , dR_G_W_dqe);
        d2R_G_W_dt2_star = d2R_G_C_dt2_star*R_C_W + R_G_C*d2R_C_W_dt2_star;
        
        E_G = repmat(R_G_W,1,1,ns); E_G_tr = permute(E_G,[2 1 3]);
        dE_dt_G = repmat(dR_G_W_dt,1,1,ns);
        
        massOffset_W = obj.rMu_W;
        ms = obj.Mu;
        connectionOffset_C = obj.connectionOffset_C;
        I_varTheta = obj.I_varTheta;
        Gamma_W = obj.Gamma_W;
        Gamma_root_G = partInformationStruct.(parentName).Gamma_G(:,1,parentConnIdx) + R_G_C*connectionOffset_C;
        dGamma_dt_root_G = partInformationStruct.(parentName).dGamma_dt_G(:,1,parentConnIdx) + dR_G_C_dt*connectionOffset_C;
        d2Gamma_dt2_root_G_star = partInformationStruct.(parentName).d2Gamma_dt2_G_star(:,1,parentConnIdx) + d2R_G_C_dt2_star*connectionOffset_C;
        dGamma_dqe_root_G_Dim3x1x1xnqe = partInformationStruct.(parentName).dGamma_dq_G(:,:,parentConnIdx,:) + utility_functions.MultiProd_(dR_G_C_dqe,connectionOffset_C);
        dGamma_dqC_root_G_Dim3x1x1xnqC = zeros(3,1,1,nqC);
        dGamma_dq_root_G_Dim3x1x1xnq = cat(4, dGamma_dqC_root_G_Dim3x1x1xnqC , dGamma_dqe_root_G_Dim3x1x1xnqe);
        
        Gamma_G = Gamma_root_G + utility_functions.MultiProd_(R_G_W,Gamma_W);
        dGamma_dt_G = dGamma_dt_root_G + utility_functions.MultiProd_(dR_G_W_dt,Gamma_W);
        d2Gamma_dt2_G_star = d2Gamma_dt2_root_G_star + utility_functions.MultiProd_(d2R_G_W_dt2_star,Gamma_W);
        d2Gamma_m_dt2_G_star =  d2Gamma_dt2_root_G_star + utility_functions.MultiProd_(d2R_G_W_dt2_star,massOffset_W);
        dGamma_dq_G_Dim3x1xnsxnq2nd = dGamma_dq_root_G_Dim3x1x1xnq + utility_functions.MultiProd_(dR_G_W_dq,Gamma_W);
        dGamma_m_dq_G_Dim3x1x1xnq2nd = dGamma_dq_root_G_Dim3x1x1xnq + utility_functions.MultiProd_(dR_G_W_dq,massOffset_W);
        
        dvarTheta_dt_root_G = R_G_C*depsilon_dt_C + partInformationStruct.(parentName).dvarTheta_dt_G(:,1,parentConnIdx);
        d2varTheta_dt2_root_G_star = partInformationStruct.(parentName).d2varTheta_dt2_G_star(:,1,parentConnIdx);
        dvarTheta_dqe_root_G_Dim3x1x1xnqe = partInformationStruct.(parentName).dvarTheta_dq_G(:,:,parentConnIdx,:);
        dvarTheta_dqC_root_G_Dim3x1x1xnqC = utility_functions.MultiProd_(R_G_C,depsilon_dqC_C);
        dvarTheta_dq_root_G_Dim3x1x1xnq = cat(4 , dvarTheta_dqC_root_G_Dim3x1x1xnqC , dvarTheta_dqe_root_G_Dim3x1x1xnqe);
        
        dvarTheta_dt_G = repmat(dvarTheta_dt_root_G,1,1,ns);
        d2varTheta_dt2_G_star = repmat(d2varTheta_dt2_root_G_star,1,1,ns);
        dvarTheta_m_dq_G_Dim3x1x1xnq2nd = dvarTheta_dq_root_G_Dim3x1x1xnq;
        dvarTheta_dq_G_Dim3x1xnsxnq2nd = repmat(dvarTheta_dq_root_G_Dim3x1x1xnq,1,1,ns);
        
        dGamma_m_dq_G_Dim3xnq2ndx1 = permute(dGamma_m_dq_G_Dim3x1x1xnq2nd,[1 4 3 2]);
        dGamma_m_dq_G_Dimnq2ndx3x1 = permute(dGamma_m_dq_G_Dim3x1x1xnq2nd,[4 1 3 2]);
        [dW_dq_Kinetic_Translation_ddqComponent,dW_dq_Kinetic_Translation_remainder] = static_method_groups.dPi_dq_Kinetic_Translation.get_dPi_dq_Kinetic_Translation(...
            ms, dGamma_m_dq_G_Dim3xnq2ndx1, dGamma_m_dq_G_Dimnq2ndx3x1,...
            d2Gamma_m_dt2_G_star);
        
        dvarTheta_m_dq_G_Dim3xnq2ndx1 = permute(dvarTheta_m_dq_G_Dim3x1x1xnq2nd,[1 4 3 2]);
        dvarTheta_m_dq_G_Dimnq2ndx3x1 = permute(dvarTheta_m_dq_G_Dim3x1x1xnq2nd,[4 1 3 2]);
        [dW_dq_Kinetic_Rotation_ddqComponent,dW_dq_Kinetic_Rotation_remainder] = static_method_groups.dPi_dq_Kinetic_Rotation.get_dPi_dq_Kinetic_Rotation(...
            E_G, E_G_tr, I_varTheta,...
            dvarTheta_m_dq_G_Dim3xnq2ndx1, dvarTheta_m_dq_G_Dimnq2ndx3x1,...
            dvarTheta_dt_G,d2varTheta_dt2_G_star);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% RigidPart_nonlinear Applied Loads from PvecApplied_G, MvecApplied_G, Gravitational Acceleration
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        PvecWeight_G = gravAccVec.*obj.Mu;
        % massOffset_G = MultiProd_(E_G,massOffset_I);
        
        PvecApplied_G = obj.Pvec_appliedGlobal_G + utility_functions.MultiProd_(E_G,obj.Pvec_appliedLocal_I);
        MvecApplied_G = obj.Mvec_appliedGlobal_G + utility_functions.MultiProd_(E_G,obj.Mvec_appliedLocal_I);
        
        %==========================================================================
        %Virtual work terms from applied loads
        dPi_dq_Fapplied = utility_functions.dotn(PvecApplied_G,dGamma_dq_G_Dim3x1xnsxnq2nd,1);
        dPi_dq_Mapplied = utility_functions.dotn(MvecApplied_G,dvarTheta_dq_G_Dim3x1xnsxnq2nd,1);
        dPi_dq_Grav = utility_functions.dotn(PvecWeight_G,dGamma_m_dq_G_Dim3x1x1xnq2nd,1);
        dW_dq_AppliedLoad = sum(dPi_dq_Fapplied + dPi_dq_Mapplied,3) + dPi_dq_Grav;
        %==========================================================================
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        dW_dq_part = ...
            - dW_dq_Kinetic_Rotation_remainder(:)...
            - dW_dq_Kinetic_Translation_remainder(:)...
            + dW_dq_epsilon(:)...
            + dW_dq_AppliedLoad(:);
        
        dM_dq_part = ...
            + dW_dq_Kinetic_Rotation_ddqComponent...
            + dW_dq_Kinetic_Translation_ddqComponent;
        
        
        % partInformationStruct.(rigid_part_name).E_G = E_G;
        % partInformationStruct.(rigid_part_name).dE_dt_G = dE_dt_G;
        % partInformationStruct.(rigid_part_name).d2E_dt2_G_star = d2E_dt2_G_star;
        % partInformationStruct.(rigid_part_name).dE_dq_G = dE_dq_G_Dim3x3xnsxnq2nd;
        % partInformationStruct.(rigid_part_name).dvarTheta_dq_G = dvarTheta_dq_G_Dim3x1xnsxnq2nd;
        % partInformationStruct.(rigid_part_name).dvarTheta_dt_G = dvarTheta_dt_G;
        % partInformationStruct.(rigid_part_name).Gamma_G = Gamma_G;
        % partInformationStruct.(rigid_part_name).dGamma_dt_G = dGamma_dt_G;
        % partInformationStruct.(rigid_part_name).d2Gamma_dt2_G_star = d2Gamma_dt2_G_star;
        % partInformationStruct.(rigid_part_name).dGamma_dq_G = dGamma_dq_G_Dim3x1xnsxnq2nd;
        % partInformationStruct.(rigid_part_name).q2nd_idx = q2nd_idx;
        % partInformationStruct.(rigid_part_name).dq2nd_idx = dq2nd_idx;
        partInformationStruct.(rigid_part_name).dqg2nd_to_dq2nd_idx = dqg2nd_to_dq2nd_idx;
        % partInformationStruct.(rigid_part_name).ns = ns;
        % partInformationStruct.(rigid_part_name).d2varTheta_dt2_G_star = d2varTheta_dt2_G_star;
        
        
        
        
        
        if obj.isAero == true && ~isempty(aerodynamics)
            
            if isempty(SimObject.aeroPartNames)
                aeroPartNames = [aeroPartNames {rigid_part_name}];
            end
            
            if ismember(aerodynamics,{'strip_steady','strip_unsteady'})
                
                Apm_idx = obj.Apm_idx;
                nAmid = obj.nAnodes - 1;
                E_G_pm = utility_functions.sample(E_G,Apm_idx,3);
                dE_dt_G_pm = utility_functions.sample(dE_dt_G,Apm_idx,3);
                dvarTheta_dt_G_pm = utility_functions.sample(dvarTheta_dt_G,Apm_idx,3);
                EAp_I_pm = obj.EAp_I_pm;
                EAp_G_pm = utility_functions.MultiProd_(E_G_pm,EAp_I_pm);
                dEAp_dt_G_pm = utility_functions.MultiProd_(dE_dt_G_pm,EAp_I_pm);
                partInformationStruct.(rigid_part_name).aero_cntr_pm = obj.aero_cntr_pm;
                partInformationStruct.(rigid_part_name).nsAp = obj.nAnodes-1;
                partInformationStruct.(rigid_part_name).EAp_G_pm = EAp_G_pm;
                partInformationStruct.(rigid_part_name).dEAp_dt_G_pm = dEAp_dt_G_pm;
                partInformationStruct.(rigid_part_name).Gamma_G_pm = utility_functions.sample(Gamma_G,Apm_idx,3);
                partInformationStruct.(rigid_part_name).dGamma_dt_G_pm = utility_functions.sample(dGamma_dt_G,Apm_idx,3);
               %partInformationStruct.(rigid_part_name).dGammaAlphaCP_dt_G_pm = dGammaAlphaCP_dt_G_pm;
                partInformationStruct.(rigid_part_name).dvarTheta_dt_G_pm = dvarTheta_dt_G_pm;
                partInformationStruct.(rigid_part_name).chord = obj.c_pm;
                partInformationStruct.(rigid_part_name).ApWidth = obj.ApWidths_pm;
                partInformationStruct.(rigid_part_name).beam_cntr_pm = obj.beam_cntr_pm;
                partInformationStruct.(rigid_part_name).AIC = obj.AICs;
                partInformationStruct.(rigid_part_name).dGamma_dq_G_pm = utility_functions.sample(dGamma_dq_G_Dim3x1xnsxnq2nd,Apm_idx,3);
                partInformationStruct.(rigid_part_name).dvarTheta_dq_G_pm = utility_functions.sample(dvarTheta_dq_G_Dim3x1xnsxnq2nd,Apm_idx,3);
                partInformationStruct.(rigid_part_name).aeroData = obj.aeroData;
                if ~isempty(obj.qAero)
                    partInformationStruct.(rigid_part_name).qAero_idx = SimObject.StateInfo{'Index',obj.qAero.group}{:};
                else
                    partInformationStruct.(rigid_part_name).qAero_idx = [];
                end
                
            elseif ismember(aerodynamics,{'VLM_steady'})
                
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
                Pv = Gamma_G_pn + ((aero_cntr_pn - beam_cntr_pn).*chord_pn) .* EAp_G_pn(:,1,:) ;
                Pc = Gamma_G_pm + ((control_point_pm - beam_cntr_pm).*chord_pm) .* EAp_G_pm(:,1,:) ;
                PcNorm = EAp_G_pm(:,3,:);
                dPc_dt = dGamma_dt_G_pm + ((control_point_pm - beam_cntr_pm).*chord_pm) .* dEAp_dt_G_pm(:,1,:) ;
                
                partInformationStruct.(flex_part_name).Pv = Pv;
                partInformationStruct.(flex_part_name).Pc = Pc;
                partInformationStruct.(flex_part_name).PcNorm = PcNorm;
                partInformationStruct.(flex_part_name).dPc_dt = dPc_dt;
                partInformationStruct.(flex_part_name).symmetric_plane_cell = [];
                
            end
            
        end
        
        
        
        if isequal(outputFormat,'qoi')
            
            QOI_Container = utility_functions.get_field(SimObject,['QOI_Master.QOIcontainers_struct.rigidParts.' rigid_part_name]);
            
            %Gamma_G = bsxfun(@plus,rBarA_G + R_G_A*wingRoot_offset_A,MultiProd_(R_G_W,Gamma_W));
            QOI_Container.add_qoi('Gamma_G',tidx,Gamma_G,'1:ns','\Gamma#_{[G]}','m');
            QOI_Container.add_qoi('R_C_W',tidx,reshape(R_C_W,9,1,1),'1','R_C_W#','');
            
            CoM_G = reshape(Gamma_root_G + R_G_W*massOffset_W, 1,3,1);
            CoM_info_flexPart_nonlinear(1) = obj.Mu;
            CoM_info_flexPart_nonlinear(2:4) = CoM_G;
            QOI_Container.add_qoi('CoM_G',tidx,CoM_G,'1','CoM#_{[G]}','m');
            
            QOI_Container.add_qoi('ex_G',tidx,E_G(:,1,:),'1:ns','ex#_{[G]}','');
            QOI_Container.add_qoi('ey_G',tidx,E_G(:,2,:),'1:ns','ey#_{[G]}','');
            QOI_Container.add_qoi('ez_G',tidx,E_G(:,3,:),'1:ns','ez#_{[G]}','');
            QOI_Container.add_qoi('PvecApplied_I', tidx, obj.Pvec_appliedLocal_I,'1:ns','Pvec_I','');
            
            QOI_Container.discretisationVariables.ns = obj.ns;
            QOI_Container.discretisationVariables.nt = SimObject.nt;
            QOI_Container.discretisationVariables.T = SimObject.t(end);
            
        end
        
        end
    
    end





end