

function output = f(varargin)

if nargin == 0
    disp(' ')
    disp('f.m contains the system of equations governing the full aeroelastic system.')
    disp(' ')
    disp('Input Formats:')
    disp(' ')
    disp('<Q = StateVector, t = time>')
    disp(' ')
    disp('Static simulation             | residualVector = f( Q , MasterObject , ''static'' )')
    disp('                                t is assumed to be zero when evaluating expressions')
    disp(' ')
    disp('Dynamic simulation            | dQ/dt = f( t , Q , MasterObject , ''dQ'' )')
    disp(' ')
    disp('Dynamic residual              | dPI/dQ = f( t , Q , MasterObject , ''r'' )')
    disp('                                where dPI = total virtual work from all sources')
    disp(' ')
    disp('Dynamic residual derivative   | [ dPI/dQ , d2PI/dQ2 ] = f( t , Q , MasterObject , ''r_dr_dq'' )')
    disp('                                where d2PI/dQ2 is the nQ x nQ dynamic residual derivative matrix (used for example in the Newmark Beta Method)')
    disp(' ')
    disp('Jacobian                      | d[dQ/dt]/dQ = f( t , Q , MasterObject , ''jac'' )')
    disp('                                outputs the Jacobian for the first order system')
    disp(' ')
    disp('Output Quantities of Interest | no_assignment -> f( t , Q , MasterObject , ''qoi'' , tidx )')
    disp('                                update the output handle objects contained in MasterObject')
    disp('                                tidx is defined such that t = MasterObject.t(tidx)')
    return
end

%see comments above
tidx = [];
if strcmp(varargin{4},'static')
    t = 0;
    [qStatic,qStatic_idx,SimObject,outputFormat] = varargin{:};
    Q = SimObject.IC*0;
    Q(qStatic_idx) = qStatic;
else
    [t,Q,SimObject,outputFormat] = varargin{:};
    if isequal(outputFormat,'qoi'), tidx = varargin{5}; end
end
FLAG_static = SimObject.FLAG_static;
% Q =                                                                            %[temp,temp1,temp2,temp3,temp4,str,Struct,Cell,Table] = deal([]); %#ok<ASGLU> %workspace variables used only for debugging

global mult3d_mex
mult3d_mex = SimObject.mult3d_mex;
% int_fnc = SimObject.int_fnc;
aeroPartNames = SimObject.aeroPartNames;

% gravAccVec = SimObject.grav_acc*SimObject.gravVec_G;

%see class definition NBS_MasterObject.m for detailed parameter definitions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copy properties of SimObject structure to local variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aerodynamics = SimObject.aerodynamics;
if isa(SimObject.uVec_freeStream_G,'function_handle')
    uVec_freeStream_G = SimObject.uVec_freeStream_G(t);
else
    uVec_freeStream_G = SimObject.uVec_freeStream_G; %unit vector in free stream velocity direction
end

%%

V = SimObject.V; %m/s
rho = SimObject.rho; %kg/m^3
Vinf = V.*uVec_freeStream_G;
FLAG_free_free = SimObject.FLAG_free_free;
% Gamma_Integration_Function = SimObject.Gamma_int_fnc;
qg2nd_idx = SimObject.qg2nd_idx;
dqg2nd_idx = SimObject.dqg2nd_idx;
nqg2nd = SimObject.nqg2nd;
nQ = size(Q);

dQ_Aero = zeros(nQ(1), nQ(2), length(SimObject.aeroPartNames));
if isa(Q, 'AutoDiff')
    dQ_Aero = AutoDiff(dQ_Aero, Q);
end

% if FLAG_free_free, nqr = 6; else, nqr = 0; end

%=======remark on state vector notation========
% ---- 1st order system
% Q - state vector of 1st order system

% ---- 2nd order system
% q = [qa;qs;qr] - state vector of 2nd order system
% qa - attitude states (relating to th,si,ph)
% qs - shear states
% qr - rigid states
% qf = [qa;qs] - flexible states of 2nd order system
%==============================================

rT_idx  = SimObject.rT_idx;
rR_idx  = SimObject.rR_idx;
drT_idx = SimObject.drT_idx;
drR_idx = SimObject.drR_idx;

if FLAG_free_free
    
    rBarA_G = Q(rT_idx);
    Beta_G = Q(rR_idx);
    drBarA_dt_G = Q(drT_idx);
    Omega_G = Q(drR_idx);
    
    d2rBarA_dt2_G_star = [0;0;0];
    dOmega_dt_G_star = [0;0;0];
    BetaNorm_G = max(norm(Beta_G),1e-10);
    R_G_A = utility_functions.r_matrix(Beta_G./BetaNorm_G,BetaNorm_G);
    OmegaSkew_G = utility_functions.getSkewMat(Omega_G); dOmega_dtSkew_G_star = utility_functions.getSkewMat(dOmega_dt_G_star);
    dR_G_A_dt = OmegaSkew_G*R_G_A;
    d2R_G_A_dt2_star = dOmega_dtSkew_G_star*R_G_A + OmegaSkew_G*OmegaSkew_G*R_G_A;
    
    skewX = utility_functions.getSkewMat([1;0;0]);
    skewY = utility_functions.getSkewMat([0;1;0]);
    skewZ = utility_functions.getSkewMat([0;0;1]);
    
    dR_G_A_dqr_Dim3x3x1xnqr = utility_functions.MultiProd_(cat(4,zeros(3,3,1,3),skewX,skewY,skewZ) , R_G_A);
    dvarTheta_dqr_G_Dim3x1x1xnqr = reshape([zeros(3) eye(3)],3,1,1,6);
    drBarA_G_dqr_G_Dim3x1x1xnqr = reshape([eye(3) zeros(3)],3,1,1,6);
    
elseif ~isempty(SimObject.CUSTOM_free_states)
    
    custom_free_state_map = SimObject.CUSTOM_free_states(SimObject,Q,t);
    
    rBarA_G = custom_free_state_map.rBarA_G;
    drBarA_dt_G = custom_free_state_map.drBarA_dt_G;
    Beta_G = custom_free_state_map.Beta_G;
    R_G_A = custom_free_state_map.R_G_A;
    Omega_G = custom_free_state_map.Omega_G;
    dR_G_A_dqr_Dim3x3x1xnqr = custom_free_state_map.dR_G_A_dqr_Dim3x3x1xnqr;
    dvarTheta_dqr_G_Dim3x1x1xnqr = custom_free_state_map.dvarTheta_dqr_G_Dim3x1x1xnqr;
    drBarA_G_dqr_G_Dim3x1x1xnqr = custom_free_state_map.drBarA_G_dqr_G_Dim3x1x1xnqr;
    
    d2rBarA_dt2_G_star = [0;0;0];
    dOmega_dt_G_star = [0;0;0];
    BetaNorm_G = max(norm(Beta_G, 2), 1e-10);
    OmegaSkew_G = utility_functions.getSkewMat(Omega_G); dOmega_dtSkew_G_star = utility_functions.getSkewMat(dOmega_dt_G_star);
    dR_G_A_dt = OmegaSkew_G*R_G_A;
    d2R_G_A_dt2_star = dOmega_dtSkew_G_star*R_G_A + OmegaSkew_G*OmegaSkew_G*R_G_A;
    
else
    
    if ~isempty(SimObject.prescribedMotion_fnc)
        prescribed_motion = SimObject.prescribedMotion_fnc(SimObject,t);
        
        %         beta_ = prescribed_motion.beta_;
        Beta_G = [0;0;1];
        Omega_G = prescribed_motion.Omega_G;
        OmegaSkew_G = prescribed_motion.OmegaSkew_G;
        dOmega_dt_G_star = prescribed_motion.dOmega_dt_G;
        R_G_A = prescribed_motion.R_G_A;
        dR_G_A_dt = prescribed_motion.dR_G_A_dt;
        d2R_G_A_dt2_star = prescribed_motion.d2R_G_A_dt2;
        rBarA_G = prescribed_motion.rBarA_G;
        drBarA_dt_G = prescribed_motion.drBarA_dt_G;
        d2rBarA_dt2_G_star = prescribed_motion.d2rBarA_dt2_G;
    else
        Omega_G = zeros(3,1);
        if isempty(SimObject.R_G_A_0)
            R_G_A = eye(3);
        else
            R_G_A = SimObject.R_G_A_0;
        end
        dR_G_A_dt = zeros(3);
        d2R_G_A_dt2_star = zeros(3);
        OmegaSkew_G = zeros(3);
        d2rBarA_dt2_G_star = zeros(3,1);
        dOmega_dt_G_star = zeros(3,1);
        rBarA_G = zeros(3,1);
        drBarA_dt_G = zeros(3,1);
    end
    dR_G_A_dqr_Dim3x3x1xnqr = zeros(3,3,1,0);
    dvarTheta_dqr_G_Dim3x1x1xnqr = zeros(3,1,1,0);
    drBarA_G_dqr_G_Dim3x1x1xnqr = zeros(3,1,1,0);
end

output = [];

nflexParts_nonlinear = SimObject.nflexParts_nonlinear;
nrigidParts = SimObject.nrigidParts;

CoM_info_flexPart_nonlinear = zeros(4,nflexParts_nonlinear);
CoM_info_rigidPart = zeros(4,nrigidParts);

partInformationStruct = struct;

NBS_Master_partName = SimObject.partName;
partInformationStruct.(NBS_Master_partName).dq2nd_idx = [drT_idx(:);drR_idx(:)];
partInformationStruct.(NBS_Master_partName).connection_idx_ParentObj = 1;

partInformationStruct.(NBS_Master_partName).Gamma_G = rBarA_G;
partInformationStruct.(NBS_Master_partName).dGamma_dt_G = drBarA_dt_G;
partInformationStruct.(NBS_Master_partName).d2Gamma_dt2_G_star = d2rBarA_dt2_G_star;
partInformationStruct.(NBS_Master_partName).dGamma_dq_G = drBarA_G_dqr_G_Dim3x1x1xnqr;

% partInformationStruct.(NBS_Master_partName).azimuth = beta_;
partInformationStruct.(NBS_Master_partName).dvarTheta_dt_G = Omega_G;
partInformationStruct.(NBS_Master_partName).OmegaSkew_G = OmegaSkew_G;
partInformationStruct.(NBS_Master_partName).d2varTheta_dt2_G_star = dOmega_dt_G_star;
partInformationStruct.(NBS_Master_partName).dvarTheta_dq_G = dvarTheta_dqr_G_Dim3x1x1xnqr;

partInformationStruct.(NBS_Master_partName).E_G = R_G_A;
partInformationStruct.(NBS_Master_partName).dE_dt_G = dR_G_A_dt;
partInformationStruct.(NBS_Master_partName).dE_dq_G = dR_G_A_dqr_Dim3x3x1xnqr;
partInformationStruct.(NBS_Master_partName).d2E_dt2_G_star = d2R_G_A_dt2_star;

%TODO change calling order of all aircraft subparts to reflect object hierarchy

% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================

for i_flex_part = 1:nflexParts_nonlinear
    obj = SimObject.flexParts_nonlinear_cell{i_flex_part};
    [dW_dq_part,dM_dq_part,partInformationStruct] = f_flexPart_nonlinear(obj,Q,partInformationStruct,outputFormat,tidx); %TODO reduce to more general form like f_rigidPart()
    
    flex_part_name = obj.partName;
    %dq2nd_idx = partInformationStruct.(flex_part_name).dq2nd_idx;
    dqg2nd_to_dq2nd_idx = partInformationStruct.(flex_part_name).dqg2nd_to_dq2nd_idx;
    
    dW_dqg(dqg2nd_to_dq2nd_idx,1,i_flex_part) = dW_dq_part;
    dM_dqg(dqg2nd_to_dq2nd_idx,dqg2nd_to_dq2nd_idx,i_flex_part) = dM_dq_part;
end

% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================

for i_rigid_part = 1:nrigidParts
    obj = SimObject.rigidParts_cell{i_rigid_part};
    [dW_dq_part,dM_dq_part,partInformationStruct] = f_rigidPart(obj,Q,partInformationStruct,outputFormat,tidx);
    
    rigid_part_name = obj.partName;
    %dq2nd_idx = partInformationStruct.(rigid_part_name).dq2nd_idx;
    dqg2nd_to_dq2nd_idx = partInformationStruct.(rigid_part_name).dqg2nd_to_dq2nd_idx;
    
    dW_dqg(dqg2nd_to_dq2nd_idx,1,nflexParts_nonlinear + i_rigid_part) = dW_dq_part;
    dM_dqg(dqg2nd_to_dq2nd_idx,dqg2nd_to_dq2nd_idx,nflexParts_nonlinear + i_rigid_part) = dM_dq_part;
    
end

% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================


%Global Aerodynamic Models
%==========================================================================
if ~isempty(aerodynamics)
    
    switch aerodynamics
        case {'strip_steady','strip_unsteady'}
            %TODO pre allocate aero matrices after the first iteration
            
            [EAp_G_pm_global,...
                dEAp_dt_G_pm_global,...
                Gamma_G_pm_global,...
                dGamma_dt_G_pm_global,...
                dvarTheta_dt_G_pm_global,...
                chord_pm_global,ApWidth_pm_global,...
                beam_cntr_pm_global,...
                aero_cntr_pm_global,...
                AICs_global,...
                dGamma_dqg_G_pm_global,...
                dvarTheta_dqg_G_pm_global,...
                qAero_idx_global,...
                ] = deal([]);
            
            global_idx_counter = 0;
            
            for aeroPartName_cell = aeroPartNames
                aeroPartName = aeroPartName_cell{1};
                aero_cntr_pm_global = cat(3,aero_cntr_pm_global,partInformationStruct.(aeroPartName).aero_cntr_pm);
                aeroData_global = partInformationStruct.(aeroPartName).aeroData;
                nsAp = partInformationStruct.(aeroPartName).nsAp;
                EAp_G_pm_global = cat(3,EAp_G_pm_global,partInformationStruct.(aeroPartName).EAp_G_pm);
                dEAp_dt_G_pm_global = cat(3,dEAp_dt_G_pm_global,partInformationStruct.(aeroPartName).dEAp_dt_G_pm);
                Gamma_G_pm_global = cat(3,Gamma_G_pm_global,partInformationStruct.(aeroPartName).Gamma_G_pm);
                dGamma_dt_G_pm_global = cat(3,dGamma_dt_G_pm_global,partInformationStruct.(aeroPartName).dGamma_dt_G_pm);
                %dGammaAlphaCP_dt_G_pm_global = cat(3,dGammaAlphaCP_dt_G_pm_global,partInformationStruct.(aeroPartName).dGammaAlphaCP_dt_G_pm);
                dvarTheta_dt_G_pm_global = cat(3,dvarTheta_dt_G_pm_global,partInformationStruct.(aeroPartName).dvarTheta_dt_G_pm);
                chord_pm_global = cat(3,chord_pm_global,partInformationStruct.(aeroPartName).chord);
                ApWidth_pm_global = cat(3,ApWidth_pm_global,partInformationStruct.(aeroPartName).ApWidth);
                beam_cntr_pm_global = cat(3,beam_cntr_pm_global,partInformationStruct.(aeroPartName).beam_cntr_pm);
                AICs_global = cat(3,AICs_global,partInformationStruct.(aeroPartName).AIC);
                qAero_idx_global = cat(1,qAero_idx_global,partInformationStruct.(aeroPartName).qAero_idx(:));
                
                SimObject.allParts_struct.(aeroPartName).sAp_idx_global = 1:nsAp + global_idx_counter;
                global_idx_counter = global_idx_counter + nsAp;
                
                dqg2nd_to_dq2nd_idx = partInformationStruct.(aeroPartName).dqg2nd_to_dq2nd_idx;
                
                dGamma_dqg_G_pm_part = zeros(3,1,nsAp,nqg2nd);
                dGamma_dqg_G_pm_part(:,:,:,dqg2nd_to_dq2nd_idx) = partInformationStruct.(aeroPartName).dGamma_dq_G_pm;
                dGamma_dqg_G_pm_global = cat(3,dGamma_dqg_G_pm_global,dGamma_dqg_G_pm_part);
                dvarTheta_dqg_G_pm_part = zeros(3,1,nsAp,nqg2nd);
                dvarTheta_dqg_G_pm_part(:,:,:,dqg2nd_to_dq2nd_idx) = partInformationStruct.(aeroPartName).dvarTheta_dq_G_pm;
                dvarTheta_dqg_G_pm_global = cat(3,dvarTheta_dqg_G_pm_global,dvarTheta_dqg_G_pm_part);
                
            end
            
            aero_cntr = 1/4;
            alphaCP = 3/4;
            
            dexAp_dt_G_pm_global = dEAp_dt_G_pm_global(:,1,:);
            dGammaAlphaCP_dt_G_pm_global = ...%bsxfun(@plus,drBarA_dt_G,...
                (dGamma_dt_G_pm_global + chord_pm_global.*(alphaCP - beam_cntr_pm_global).*dexAp_dt_G_pm_global);
            
            aeroOffset_global = chord_pm_global.*(aero_cntr - beam_cntr_pm_global).*EAp_G_pm_global(:,1,:); %center of pressure offset from the beam line, +ve in ex direction
            aeroOffset_skew_global = utility_functions.getSkewMat(aeroOffset_global);
            
            switch aerodynamics
                
                %=================================quasi steady strip theory=============================V3qrt,xAp,yAp,zAp
                case 'strip_steady'
                    
                    Qaero = [];
                    V3qrt = dGammaAlphaCP_dt_G_pm_global;
                    xAp = EAp_G_pm_global(:,1,:);
                    yAp = EAp_G_pm_global(:,2,:);
                    zAp = EAp_G_pm_global(:,3,:);
                    dAoA = dvarTheta_dt_G_pm_global;
                    chord = chord_pm_global;
                    ApWidth = ApWidth_pm_global;
                    AIC = AICs_global;
                    C_D0 = 0;
                    qsteady = true;
                    aeroData = aeroData_global;
                    
                    [~,~,Fqc,Mqc,Drag,alpha_global] = aerodynamic_codes.aero_stripTheory_usteady_LeishmanIndicial(...
                        Qaero,rho,Vinf,V3qrt,xAp,yAp,zAp,dAoA,chord,ApWidth,AIC,C_D0,aeroData,qsteady);
                    
                case 'strip_unsteady'
                    
                    Qaero_idx = qAero_idx_global;
                    Qaero = reshape(Q(Qaero_idx),2,1,[]);
                    %             Vinf = V*uVec_freeStream_G;
                    V3qrt = dGammaAlphaCP_dt_G_pm_global;
                    xAp = EAp_G_pm_global(:,1,:);
                    yAp = EAp_G_pm_global(:,2,:);
                    zAp = EAp_G_pm_global(:,3,:);
                    dAoA = dvarTheta_dt_G_pm_global;
                    chord = chord_pm_global;
                    ApWidth = ApWidth_pm_global;
                    AIC = AICs_global;
                    C_D0 = 0;
                    qsteady = false;
                    aeroCoeff2D = [];
                    
                    [dQaero,~,Fqc,Mqc,Drag,alpha_global] = aerodynamic_codes.aero_stripTheory_usteady_LeishmanIndicial(...
                        Qaero,rho,Vinf,V3qrt,xAp,yAp,zAp,dAoA,chord,ApWidth,AIC,C_D0,aeroCoeff2D,qsteady);
                    
                    dQ_Aero(Qaero_idx,1) = dQaero(:);
                    
            end
            
        case 'WT'
            R_A_G = R_G_A.';
            %TODO pre allocate aero matrices after the first iteration
            [dQ_Aero, ...
                Fqc, Mqc, Drag,...
                Gamma_G_pm_global, Gamma_A_pm_global, ...
                dGamma_dqg_G_pm_global, ...
                dvarTheta_dqg_G_pm_global, EAp_G_pm_global,...
                aeroOffset_global, aeroOffset_skew_global, ApWidth_pm_global] ...
                = aero.wt_aero(SimObject,partInformationStruct,dQ_Aero,R_A_G,rBarA_G,Omega_G, t);
    end
    
    PvecAero_G_pm_global = Fqc + Drag;
    PvecAero_Gamma_G = Gamma_G_pm_global + EAp_G_pm_global(:,1,:).*aeroOffset_global;
    MvecAero_G_pm_global = Mqc + utility_functions.MultiProd_(aeroOffset_skew_global,(Fqc + Drag));
    
else
    
    [dGamma_dqg_G_pm_global,...
        dvarTheta_dqg_G_pm_global,...
        harmonic_load,...
        ] = deal([]);
    
    global_idx_counter = 0;
    
    for aeroPartName_cell = aeroPartNames
        
        aeroPartName = aeroPartName_cell{1};
        dqg2nd_to_dq2nd_idx = partInformationStruct.(aeroPartName).dqg2nd_to_dq2nd_idx;
        harmonic_load = partInformationStruct.(aeroPartName).harmonic_load_global;
        nsAp = partInformationStruct.(aeroPartName).nsAp;
        dGamma_dqg_G_pm_part = zeros(3,1,nsAp,nqg2nd);
        dGamma_dqg_G_pm_part(:,:,:,dqg2nd_to_dq2nd_idx) = partInformationStruct.(aeroPartName).dGamma_dq_G_pm;
        dGamma_dqg_G_pm_global = cat(3,dGamma_dqg_G_pm_global,dGamma_dqg_G_pm_part);
        dvarTheta_dqg_G_pm_part = zeros(3,1,nsAp,nqg2nd);
        dvarTheta_dqg_G_pm_part(:,:,:,dqg2nd_to_dq2nd_idx) = partInformationStruct.(aeroPartName).dvarTheta_dq_G_pm;
        dvarTheta_dqg_G_pm_global = cat(3,dvarTheta_dqg_G_pm_global,dvarTheta_dqg_G_pm_part);
        
        SimObject.allParts_struct.(aeroPartName).sAp_idx_global = 1:nsAp + global_idx_counter;
        global_idx_counter = global_idx_counter + nsAp;
    end
    
    if ~isempty(harmonic_load)
        PvecAero_G_pm_global = zeros(3,1,nsAp);
        harmonic_load_pm = sin(SimObject.ff_h*t).*harmonic_load;
        %     PvecAero_G_pm_global(:,1,end) = PvecAero_G_pm_global(:,1,end) + harmonic_load_pm;
        PvecAero_G_pm_global = PvecAero_G_pm_global +  harmonic_load_pm;
    else
        PvecAero_G_pm_global = [];
    end
    
    PvecAero_Gamma_G = [];
    MvecAero_G_pm_global = [];
    
end
%
% F_SUM = F_SUM + PvecAero_G_pm_global;
% M_SUM = M_SUM + MvecAero_G_pm_global;

%==========================================================================
%Virtual work terms from aerodynamic loads
if isempty(PvecAero_G_pm_global)
    dPi_dq_Faero_global = 0;
else
    dPi_dq_Faero_global = utility_functions.dotn(PvecAero_G_pm_global,dGamma_dqg_G_pm_global,1);
end
if isempty(MvecAero_G_pm_global)
    dPi_dq_Maero_global = 0;
else
    dPi_dq_Maero_global = utility_functions.dotn(MvecAero_G_pm_global,dvarTheta_dqg_G_pm_global,1);
end
dW_dq_AeroLoad_global = sum(dPi_dq_Faero_global + dPi_dq_Maero_global,3);
%==========================================================================


%==========================================================================
dW_dqg(:,:,end+1) = dW_dq_AeroLoad_global;
%==========================================================================

% Jacobian_flag = [];

dW_dqg_sum = sum(dW_dqg,3);

% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================

switch outputFormat
    case 'dQ'
        dQ = Q*0; %initialise 1st order state derivative
        
        %relate 1st and 2nd derivatives for second order variables
        if SimObject.FLAG_dynControl && SimObject.type == 3
            dQ(qg2nd_idx) = Q(dqg2nd_idx(1:end-1));
        else
            dQ(qg2nd_idx) = Q(dqg2nd_idx);
        end
        
        dM_dqg_sum = sum(dM_dqg,3);% + sum(dM_dqg_rigidPart,3);
        
        if ~FLAG_free_free
            dQ(dqg2nd_idx) = dM_dqg_sum\dW_dqg_sum;
        else
            %special 1st order derivative relation for rigid rotation
            Tmat_tr = utility_functions.T_matrix(Beta_G/BetaNorm_G,BetaNorm_G).'; %<- check transpose relation !!!!!!!!!!!!!!!!!!!!!!!!!
            qrR_idx = SimObject.StateInfo.qRigidR{4};
            dqrR_idx = SimObject.StateInfo.dqRigidR{4};
            dQ(qrR_idx) = Tmat_tr\Q(dqrR_idx);
            
            dQ(dqg2nd_idx) = dM_dqg_sum\dW_dqg_sum;
        end
        
        dQ = dQ+sum(dQ_Aero,3);
        output = dQ;
        %/////////////////////////
        %=========================
        
    case 'static'
        
        %         output = Q*0; %initialise 1st order static residual vector
        %         output(qg2nd_idx) = dW_dqg_sum;
        %         output = output + sum(Q_Aero,3);
        
        output = dW_dqg_sum;
        
        %    dW_dqg_static_sum = sum(dW_dqg_static,3);
        %    output = dW_dqg_static_sum;
        
    case 'r'
        
        d2q_dt2 = extraInput;
        residual = -M*d2q_dt2-w-K-C+dPiF_dq+dPiM_dq;
        output = residual;
        
    case 'r_dr_dq'
        %         d2q_dt2 = extraInput.d2q_dt2;
        %         beta2 = extraInput.beta2;
        %         gamma2 = extraInput.gamma2;
        %         delta_t = extraInput.delta_t;
        %         residual = -M*d2q_dt2-w-K-C+dPiF_dq+dPiM_dq;
        %         Jacobian_flag = 'dr_dq';
        
    case 'jac'
        
        %         Jacobian_flag = 'jac';
        
    case 'qoi'
        QOI_Container = utility_functions.get_field(SimObject,['QOI_Master.QOIcontainers_struct.',NBS_Master_partName]); %TODO read properties(QOI_Master) to get 'QOIcontainers_struct' string
        
        CoM_info = [CoM_info_flexPart_nonlinear , CoM_info_rigidPart];
        aircraftMass = sum(CoM_info(1,:));
        QOI_Container.add_qoi('aircraftMass',tidx,aircraftMass,'1','Aircraft Mass','kg');
        
        CoM_G = sum(bsxfun(@times,CoM_info(1,:),CoM_info(2:4,:)),2)/aircraftMass;
        QOI_Container.add_qoi('CoM_G',tidx,CoM_G,'1','CoM#_{[G]}','m');
        QOI_Container.add_qoi('R_G_A_flat',tidx,reshape(R_G_A,[9 1]),'1','R_{G,A}#','');
        %         QOI_Container.add_qoi('ForcesSum',tidx,F_SUM,'1:nsAp','ForcesSum','N','GlobalAeroQuantity',true);
        %         QOI_Container.add_qoi('MomentSum',tidx,M_SUM,'1:nsAp','MomentSum','Nm','GlobalAeroQuantity',true);
        QOI_Container.add_qoi('ex_G',tidx,R_G_A(:,1,:),'1','ex#_{[G]}','');
        QOI_Container.add_qoi('ey_G',tidx,R_G_A(:,2,:),'1','ey#_{[G]}','');
        QOI_Container.add_qoi('ez_G',tidx,R_G_A(:,3,:),'1','ez#_{[G]}','');
        
        nsAp = size(PvecAero_G_pm_global,3);
        
        if ~isempty(SimObject.aeroPartNames) && ~isempty(aerodynamics)
            %             QOI_Container.add_qoi('C_p', tidx,cp,'1','C_p','[]');
            %             QOI_Container.add_qoi('dCp_da', tidx,dcp_da,'1','dC_{p}/da','[]');
            %             QOI_Container.add_qoi('EAp_G_pm', tidx,reshape(EAp_G_pm_global,9,1,[]),'1:nsAp','EAp_{[G]}','[]','GlobalAeroQuantity',true);
            %             QOI_Container.add_qoi('EAp_A_pm', tidx,reshape(EAp_A_pm_global,9,1,[]),'1:nsAp','EAp_{[A]}','[]','GlobalAeroQuantity',true);
            QOI_Container.add_qoi('Gamma_G_pm', tidx,Gamma_G_pm_global,'1:nsAp','Gamma_{[G]}_pm','[]','GlobalAeroQuantity',true);
            QOI_Container.add_qoi('Gamma_A_pm', tidx,Gamma_A_pm_global,'1:nsAp','Gamma_{[A]}_pm','[]','GlobalAeroQuantity',true);
            QOI_Container.add_qoi('Aero_Forces_G' ,tidx,PvecAero_G_pm_global,'1:nsAp','AeroForce#_{[G]}','N','GlobalAeroQuantity',true);
            PvecAero_A_pm_global = pagemtimes(R_A_G,PvecAero_G_pm_global);
            QOI_Container.add_qoi('Aero_Forces_A' ,tidx,PvecAero_A_pm_global,'1:nsAp','AeroForce#_{[A]}','N','GlobalAeroQuantity',true);
            QOI_Container.add_qoi('Aero_ForcePerSpan_G', tidx,PvecAero_G_pm_global./ApWidth_pm_global,'1:nsAp','AeroForcePerSpan#_{[G]}','N/m','GlobalAeroQuantity',true);
            QOI_Container.add_qoi('Aero_Forces_Gamma_G', tidx, PvecAero_Gamma_G,'1:nsAp','AeroForce \Gamma#_{[G]}','m');
            QOI_Container.add_qoi('Aero_Moments_G',tidx, MvecAero_G_pm_global,'1:nsAp','AeroMoment#_{[G]}','Nm','GlobalAeroQuantity',true);
            MvecAero_A_pm_global = pagemtimes(R_A_G,MvecAero_G_pm_global);
            QOI_Container.add_qoi('Aero_Moments_A',tidx, MvecAero_A_pm_global,'1:nsAp','AeroMoment#_{[A]}','Nm','GlobalAeroQuantity',true);
            QOI_Container.add_qoi('Aero_MomentPerSpan_G', tidx, MvecAero_G_pm_global./ApWidth_pm_global,'1:nsAp','AeroMomentPerSpan#_{[G]}','N','GlobalAeroQuantity',true);
            QOI_Container.add_qoi('Net_Lift', tidx, sum(PvecAero_G_pm_global(3,1,:)),'1','NetLift','N');
            if exist('alpha_global','var')
                alpha_degrees = alpha_global*180/pi;
                QOI_Container.add_qoi('Angle_Of_Attack',tidx,reshape(alpha_degrees,1,1,[]),'1:nsAp','Angle Of Attack','deg');
            end
            if exist('CL','var')
                QOI_Container.add_generic_qoi('CL' ,tidx,reshape(CL,1,1,[]),'1:nsAp','CL','');
            end
        end
        
        QOI_Container.discretisationVariables.nsAp = nsAp;
        QOI_Container.discretisationVariables.nt = SimObject.nt;
        QOI_Container.discretisationVariables.ns = 1;
        
end




end

% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================

% TODO establish variable name convention:
%
% [variableName]_[sample points]_[coordinate system]_[extra suffix]
%
% sample points
% none: default, all ns points
% np: sampled at aero panel nodes
% mp: sampled at aero panel mid points
%
% coordinate system
% none: invariant quantity
% I: Intrinsic coordinate system
% WE: Euler angle reference system
% W: Wing root coordinate system
% A: Aircraft coordinate system
% G: Global coordinate system
%
% extra suffix
% none: -
% tr: transpose (dim1 <-> dim2) of indicated variable
% cell: cell array of indicated variable
%
% dimension
% [variableName, '_', dim1, dim2, dim3]: where dim1/2/3 is either a number or a two character string


%dq/dqr/dqr123/dqr456
