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

                                                                           %[temp,temp1,temp2,temp3,temp4,str,Struct,Cell,Table] = deal([]); %#ok<ASGLU> %workspace variables used only for debugging

global mult3d_mex
mult3d_mex = SimObject.mult3d_mex;
int_fnc = SimObject.int_fnc;
aeroPartNames = SimObject.aeroPartNames;
gravAccVec = SimObject.grav_acc*SimObject.gravVec_G;

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
V = SimObject.V; %m/s
rho = SimObject.rho; %kg/m^3
Vinf = V*uVec_freeStream_G;
FLAG_free_free = SimObject.FLAG_free_free;
Gamma_Integration_Function = SimObject.Gamma_int_fnc;
StateInfo = SimObject.StateInfo;
qg2nd_idx = SimObject.qg2nd_idx;
dqg2nd_idx = SimObject.dqg2nd_idx;
nqg2nd = SimObject.nqg2nd;
nQ = numel(Q);
dQ_Aero = zeros(nQ,1,SimObject.nParts);
if FLAG_free_free, nqr = 6; else, nqr = 0; end

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

rT_idx  = StateInfo{'Index','qRigidT'}{:};
rR_idx  = StateInfo{'Index','qRigidR'}{:};
drT_idx = StateInfo{'Index','dqRigidT'}{:};
drR_idx = StateInfo{'Index','dqRigidR'}{:};

if FLAG_free_free
    
    rBarA_G = Q(rT_idx);
    Beta_G = Q(rR_idx);
    drBarA_dt_G = Q(drT_idx);
    Omega_G = Q(drR_idx);
    
    d2rBarA_dt2_G_star = [0;0;0];
    dOmega_dt_G_star = [0;0;0];
    BetaNorm_G = max(norm(Beta_G),1e-10);
    R_G_A = r_matrix(Beta_G./BetaNorm_G,BetaNorm_G);
    OmegaSkew_G = getSkewMat(Omega_G); dOmega_dtSkew_G_star = getSkewMat(dOmega_dt_G_star);
    dR_G_A_dt = OmegaSkew_G*R_G_A;
    d2R_G_A_dt2_star = dOmega_dtSkew_G_star*R_G_A + OmegaSkew_G*OmegaSkew_G*R_G_A;
    
    skewX = getSkewMat([1;0;0]);
    skewY = getSkewMat([0;1;0]);
    skewZ = getSkewMat([0;0;1]);
    
    dR_G_A_dqr_Dim3x3x1xnqr = MultiProd_(cat(4,zeros(3,3,1,3),skewX,skewY,skewZ) , R_G_A);
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
    BetaNorm_G = max(norm(Beta_G),1e-10);
    OmegaSkew_G = getSkewMat(Omega_G); dOmega_dtSkew_G_star = getSkewMat(dOmega_dt_G_star);
    dR_G_A_dt = OmegaSkew_G*R_G_A;
    d2R_G_A_dt2_star = dOmega_dtSkew_G_star*R_G_A + OmegaSkew_G*OmegaSkew_G*R_G_A;
    
else
    
    if ~isempty(SimObject.prescribedMotion_fnc)
        prescribed_motion = SimObject.prescribedMotion(SimObject,t);
        
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
        R_G_A = eye(3);
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

partInformationStruct.(NBS_Master_partName).dvarTheta_dt_G = Omega_G;
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
    [dW_dq_part,dM_dq_part,partInformationStruct] = f_flexPart_nonlinear(obj,Q,partInformationStruct,outputFormat,tidx,nqr,dqg2nd_idx,R_G_A,dR_G_A_dqr_Dim3x3x1xnqr,dR_G_A_dt,d2R_G_A_dt2_star,dvarTheta_dqr_G_Dim3x1x1xnqr,Omega_G,OmegaSkew_G,dOmega_dt_G_star,rBarA_G,drBarA_dt_G,d2rBarA_dt2_G_star,drBarA_G_dqr_G_Dim3x1x1xnqr); %TODO reduce to more general form like f_rigidPart()

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
    
aero_cntr = 1/4;
alphaCP = 3/4;

if ismember(aerodynamics,{'strip_steady','strip_unsteady'})
    %TODO pre allocate aero matrices after the first iteration
    
    [EAp_G_pm_global,...
        dEAp_dt_G_pm_global,...
        Gamma_G_pm_global,...
        dGamma_dt_G_pm_global,...
        dvarTheta_dt_G_pm_global,...
        chord_pm_global,ApWidth_pm_global,...
        beam_cntr_pm_global,...
        AICs_global,...
        dGamma_dqg_G_pm_global,...
        dvarTheta_dqg_G_pm_global,...
        qAero_idx_global,...
        aeroCoeff2D_global] = deal([]);
    
    global_idx_counter = 0;
    
    for aeroPartName_cell = aeroPartNames
        
        aeroPartName = aeroPartName_cell{1};
       
        aeroCoeff2D_global = cat(3,aeroCoeff2D_global,partInformationStruct.(aeroPartName).aeroCoeff);
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
    
    dexAp_dt_G_pm_global = dEAp_dt_G_pm_global(:,1,:);
    dGammaAlphaCP_dt_G_pm_global = ...%bsxfun(@plus,drBarA_dt_G,...
                           (dGamma_dt_G_pm_global + bsxfun(@times,chord_pm_global.*(alphaCP - beam_cntr_pm_global),dexAp_dt_G_pm_global));
    
    aeroOffset_global = bsxfun(@times, chord_pm_global.*bsxfun(@plus,aero_cntr, -beam_cntr_pm_global) , EAp_G_pm_global(:,1,:)); %center of pressure offset from the beam line, +ve in ex direction
    aeroOffset_skew_global = getSkewMat(aeroOffset_global);
    
    switch aerodynamics
        %=================================quasi steady strip theory=============================V3qrt,xAp,yAp,zAp
        case 'strip_steady'
            
            Qaero = [];
            V3qrt = dGammaAlphaCP_dt_G_pm_global;
            xAp = EAp_G_pm_global(:,1,:);
            yAp = EAp_G_pm_global(:,2,:);
            zAp = EAp_G_pm_global(:,3,:);
            Omega = dvarTheta_dt_G_pm_global;
            chord = chord_pm_global;
            ApWidth = ApWidth_pm_global;
            AIC = AICs_global;
            C_D0 = 0;
            qsteady = true;
            aeroCoeff2D = aeroCoeff2D_global;
            
            [dQaero,Qaero,Fqc,Mqc,Drag,alpha_global] = aero_stripTheory_usteady_LeishmanIndicial(Qaero,rho,Vinf,V3qrt,xAp,yAp,zAp,Omega,chord,ApWidth,AIC,C_D0,qsteady,aeroCoeff2D);
            
        case 'strip_unsteady'
            
            Qaero_idx = qAero_idx_global;
            Qaero = reshape(Q(Qaero_idx),2,1,[]);
            Vinf = V*uVec_freeStream_G;
            V3qrt = dGammaAlphaCP_dt_G_pm_global;
            xAp = EAp_G_pm_global(:,1,:);
            yAp = EAp_G_pm_global(:,2,:);
            zAp = EAp_G_pm_global(:,3,:);
            Omega = dvarTheta_dt_G_pm_global;
            chord = chord_pm_global;
            ApWidth = ApWidth_pm_global;
            AIC = AICs_global;
            C_D0 = 0;
            qsteady = false;
            
            [dQaero,Qaero,Fqc,Mqc,Drag,alpha_global] = aero_stripTheory_usteady_LeishmanIndicial(Qaero,rho,Vinf,V3qrt,xAp,yAp,zAp,Omega,chord,ApWidth,AIC,C_D0,qsteady);
            
             dQ_Aero(Qaero_idx,1) = dQaero(:);
            
    end
    
end

if ismember(aerodynamics,{'VLM_steady'})
    
    nAeroParts = numel(aeroPartNames);
    
    %TODO pre allocate the cell arrays populated from partInformationStruct
    %[] = deal(cell(1,1,nAeroParts));
    
    for i_ = 1:nAeroParts
        aeroPartName = aeroPartNames{i_};
        
%         nsAp = partInformationStruct.(aeroPartName).nsAp;
%         EAp_G_pm_global = cat(3,EAp_G_pm_global,partInformationStruct.(aeroPartName).EAp_G_pm);
%         dEAp_dt_G_pm_global = cat(3,dEAp_dt_G_pm_global,partInformationStruct.(aeroPartName).dEAp_dt_G_pm);
%         dGamma_dt_G_pm_global = cat(3,dGamma_dt_G_pm_global,partInformationStruct.(aeroPartName).dGamma_dt_G_pm);
%         %dGammaAlphaCP_dt_G_pm_global = cat(3,dGammaAlphaCP_dt_G_pm_global,partInformationStruct.(aeroPartName).dGammaAlphaCP_dt_G_pm);
%         dvarTheta_dt_G_pm_global = cat(3,dvarTheta_dt_G_pm_global,partInformationStruct.(aeroPartName).dvarTheta_dt_G_pm);
%         chord_pm_global = cat(3,chord_pm_global,partInformationStruct.(aeroPartName).chord);
%         ApWidth_pm_global = cat(3,ApWidth_pm_global,partInformationStruct.(aeroPartName).ApWidth);
%         beam_cntr_pm_global = cat(3,beam_cntr_pm_global,partInformationStruct.(aeroPartName).beam_cntr_pm);
%         AICs_global = cat(3,AICs_global,partInformationStruct.(aeroPartName).AIC);
%         
%         dqg2nd_to_dq2nd_idx = partInformationStruct.(aeroPartName).dqg2nd_to_dq2nd_idx;
%         
%         dGamma_dqg_G_pm_part = zeros(3,1,nsAp,nqg2nd);
%         dGamma_dqg_G_pm_part(:,:,:,dqg2nd_to_dq2nd_idx) = partInformationStruct.(aeroPartName).dGamma_dq_G_pm;
%         dGamma_dqg_G_pm_global = cat(3,dGamma_dqg_G_pm_global,dGamma_dqg_G_pm_part);
%         dvarTheta_dqg_G_pm_part = zeros(3,1,nsAp,nqg2nd);
%         dvarTheta_dqg_G_pm_part(:,:,:,dqg2nd_to_dq2nd_idx) = partInformationStruct.(aeroPartName).dvarTheta_dq_G_pm;
%         dvarTheta_dqg_G_pm_global = cat(3,dvarTheta_dqg_G_pm_global,dvarTheta_dqg_G_pm_part);


        Pv_cell{i_} = partInformationStruct.(aeroPartName).Pv;
        Pc_cell{i_} = partInformationStruct.(aeroPartName).Pc;
        PcNorm_cell{i_} = partInformationStruct.(aeroPartName).PcNorm;
        dPc_dt_cell{i_} = partInformationStruct.(aeroPartName).dPc_dt;
        symmetric_plane_cell{i_} = partInformationStruct.(aeroPartName).symmetric_plane_cell;
        
    end
    
%     dexAp_dt_G_pm_global = dEAp_dt_G_pm_global(:,1,:);
%     dGammaAlphaCP_dt_G_pm_global = bsxfun(@plus,drBarA_dt_G,...
%                            dGamma_dt_G_pm_global + bsxfun(@times,chord_pm_global.*(alphaCP - beam_cntr_pm_global),dexAp_dt_G_pm_global));
%     
%     aeroOffset_global = bsxfun(@times, chord_pm_global.*bsxfun(@plus,aero_cntr, -beam_cntr_pm_global) , EAp_G_pm_global(:,1,:)); %center of pressure offset from the beam line, +ve in ex direction
%     aeroOffset_skew_global = getSkewMat(aeroOffset_global);
    
    
    
    
    
    
    switch aerodynamics
            
        case 'VLM_steady'
            
            %variables that already exist
            %rho
            %uVec_freeStream
            %Vinf
            %dPc_dt_cell
            %Pv_cell
            %Pc_cell
            %PcNorm_cell
            %symmetric_plane_cell
            
            [Fqc,Mqc] = aero_VLM_qsteady(rho,uVec_freeStream_G,V,dPc_dt_cell,Pv_cell,Pc_cell,PcNorm_cell,symmetric_plane_cell,L);
    end
    
end

PvecAero_G_pm_global = Fqc;
PvecAero_Gamma_G = Gamma_G_pm_global + EAp_G_pm_global(:,1,:).*aeroOffset_global;
MvecAero_G_pm_global = Mqc + MultiProd_(aeroOffset_skew_global,Fqc);

else
    
    PvecAero_G_pm_global = [];
    PvecAero_Gamma_G = [];
    MvecAero_G_pm_global = [];
    
end

%==========================================================================
%Virtual work terms from aerodynamic loads
if isempty(PvecAero_G_pm_global)
    dPi_dq_Faero_global = 0;
else
    dPi_dq_Faero_global = dotn(PvecAero_G_pm_global,dGamma_dqg_G_pm_global,1);
end
if isempty(MvecAero_G_pm_global)
    dPi_dq_Maero_global = 0;
else
    dPi_dq_Maero_global = dotn(MvecAero_G_pm_global,dvarTheta_dqg_G_pm_global,1);
end
dW_dq_AeroLoad_global = sum(dPi_dq_Faero_global + dPi_dq_Maero_global,3);
%==========================================================================


%==========================================================================
dW_dqg(:,:,end+1) = dW_dq_AeroLoad_global;
%==========================================================================

Jacobian_flag = [];

dW_dqg_sum = sum(dW_dqg,3);

% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================





switch outputFormat
    
    case 'dQ'

    dQ = Q*0; %initialise 1st order state derivative
    
    %relate 1st and 2nd derivatives for second order variables
    dQ(qg2nd_idx) = Q(dqg2nd_idx);
    
    dM_dqg_sum = sum(dM_dqg,3);% + sum(dM_dqg_rigidPart,3);
    
    
    if ~FLAG_free_free
%         dQ(1:nq) = Q(nq+1:2*nq);
%         dQ(nq+1:2*nq) = dM_dqg\dW_dqg;
        dQ(dqg2nd_idx) = dM_dqg_sum\dW_dqg_sum;
    else
        %special 1st order derivative relation for rigid rotation
        Tmat_tr = T_matrix(Beta_G/BetaNorm_G,BetaNorm_G).'; %<- check transpose relation !!!!!!!!!!!!!!!!!!!!!!!!!
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
    
    % output = Q*0; %initialise 1st order static residual vector
    
    output = dW_dqg_sum;
    
%    dW_dqg_static_sum = sum(dW_dqg_static,3);
%    output = dW_dqg_static_sum;
    
    case 'r'
    
    d2q_dt2 = extraInput;
    residual = -M*d2q_dt2-w-K-C+dPiF_dq+dPiM_dq;
    output = residual;

    case 'r_dr_dq'
    d2q_dt2 = extraInput.d2q_dt2;
    beta2 = extraInput.beta2;
    gamma2 = extraInput.gamma2;
    delta_t = extraInput.delta_t;
    residual = -M*d2q_dt2-w-K-C+dPiF_dq+dPiM_dq;
    Jacobian_flag = 'dr_dq';
    
    case 'jac'
    
    Jacobian_flag = 'jac';
    
    case 'qoi'
    
    QOI_Container = get_field(SimObject,['QOI_Master.QOIcontainers_struct.',NBS_Master_partName]); %TODO read properties(QOI_Master) to get 'QOIcontainers_struct' string
    CoM_info = [CoM_info_flexPart_nonlinear , CoM_info_rigidPart];
    aircraftMass = sum(CoM_info(1,:));
    QOI_Container.add_qoi('aircraftMass',tidx,aircraftMass,'1','Aircraft Mass','kg');
    CoM_G = sum(bsxfun(@times,CoM_info(1,:),CoM_info(2:4,:)),2)/aircraftMass;
    QOI_Container.add_qoi('CoM_G',tidx,CoM_G,'1','CoM#_{[G]}','m');
    QOI_Container.add_qoi('R_G_A_flat',tidx,reshape(R_G_A,[9 1]),'1','R_{G,A}#','');
    
    QOI_Container.add_qoi('ex_G',tidx,R_G_A(:,1,:),'1','ex#_{[G]}','');
    QOI_Container.add_qoi('ey_G',tidx,R_G_A(:,2,:),'1','ey#_{[G]}','');
    QOI_Container.add_qoi('ez_G',tidx,R_G_A(:,3,:),'1','ez#_{[G]}','');
    
    nsAp = size(PvecAero_G_pm_global,3);
    
    if obj.isAero == true && ~isempty(aerodynamics)
        QOI_Container.add_qoi('Aero_Forces_G' ,tidx,PvecAero_G_pm_global,'1:nsAp','AeroForce#_{[G]}','N','GlobalAeroQuantity',true);
        QOI_Container.add_qoi('Aero_ForcePerSpan_G' ,tidx,PvecAero_G_pm_global.*ApWidth_pm_global,'1:nsAp','AeroForcePerSpan#_{[G]}','N/m','GlobalAeroQuantity',true);
        QOI_Container.add_qoi('Aero_Forces_Gamma_G' ,tidx,PvecAero_Gamma_G,'1:nsAp','AeroForce \Gamma#_{[G]}','m');
        QOI_Container.add_qoi('Aero_Moments_G' ,tidx,MvecAero_G_pm_global,'1:nsAp','AeroMoment#_{[G]}','Nm');
        QOI_Container.add_qoi('Aero_MomentPerSpan_G' ,tidx,MvecAero_G_pm_global.*ApWidth_pm_global,'1:nsAp','AeroMomentPerSpan#_{[G]}','N');
        QOI_Container.add_qoi('Net_Lift' ,tidx,sum(PvecAero_G_pm_global(3,1,:)),'1','NetLift','N');
        if exist('alpha','var')
            alpha_degrees = alpha_global*180/pi;
            QOI_Container.add_qoi('Angle_Of_Attack' ,tidx,reshape(alpha_degrees,1,1,[]),'1:nAp','Angle Of Attack','deg');
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





function [MOMENT_xi , FORCE_xi] = material_law( KAPPA_I , ~ , KAPPA_0_I , TAU , ~ , TAU_0 , Linear_Stiffness_Matrix , ~)
if isempty(TAU)
    xi = [KAPPA_I - KAPPA_0_I];
    MOMENT_xi = -MultiProd_(Linear_Stiffness_Matrix{1}(1:3,1:3,:),xi);
    FORCE_xi = zeros(3,0,0);
else
    xi = [KAPPA_I - KAPPA_0_I ; TAU - TAU_0];
    MOMENT_xi = -MultiProd_(Linear_Stiffness_Matrix{1},xi);
    FORCE_xi  = -MultiProd_(Linear_Stiffness_Matrix{2},xi);
end
end

function [MOMENT_dxidt , FORCE_dxidt] = damping_law( ~ , dKAPPA_dt_I , ~ , ~ , dTAU_dt , ~ , Linear_Damping_Matrix , ~)
if isempty(dTAU_dt)
    dxidt = dKAPPA_dt_I;
    MOMENT_dxidt = -MultiProd_(Linear_Damping_Matrix{1}(1:3,1:3,:),dxidt);
    FORCE_dxidt  = zeros(3,0,0);
else
    dxidt = [dKAPPA_dt_I ; dTAU_dt];
    MOMENT_dxidt = -MultiProd_(Linear_Damping_Matrix{1},dxidt);
    FORCE_dxidt  = -MultiProd_(Linear_Damping_Matrix{2},dxidt);
end
end





% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================
% ======================================================================================================================================================================================================================================





%%

 %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%
  %========================= nested functions ============================%
   %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%

    function M_W = map_WEtoW(M_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD)
        %M_W = MultiProd_(R_W_WE,MultiProd_(M_WE,R_W_WE_tr,[1 2]),[1 2]);%(3)x(3)x(ns)
        
        %if the rotation R_W_WE is modulo 90 degrees use following code
        %much quicker than MultiProd_ line
        smw = Rs_W_WE_tr; pmw = Rv_W_WE_tr;
        M_W = M_WE;
        for ii_ = 1:3
            for jj_ = 1:3
                M_W(ii_,jj_,:) = smw(ii_)*smw(jj_)*M_WE(pmw(ii_),pmw(jj_),:);
            end
        end
        if ~isempty(TD)
            M_W = MultiProd_(M_W,TD);
        end
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function M_W_flat = map_WEtoW_stackDim2(M_WE_flat,Rs_W_WE_tr,Rv_W_WE_tr,TD)
        %mapping performing the transform M_W = R_W_WE*M_WE*R_W_WE_tr
        %where M_WE consists of a group of 3x3xns matrices concatinated along dim 2
        %rotation R_W_WE must consist of combined 90 degree rotations
        %    i.e. all entries of R_W_WE are 1 or -1
        
        smw = Rs_W_WE_tr; pmw = Rv_W_WE_tr; ns = size(M_WE_flat,3);
        M_W_flat = M_WE_flat;
        for ii_ = 1:3
            for jj_ = 1:3
                M_W_flat(ii_+3*(jj_-1),:,:) = smw(ii_)*smw(jj_)*M_WE_flat(pmw(ii_)+3*(pmw(jj_)-1),:,:);
            end
        end

        if ~isempty(TD)
            M_W = reshape(M_W_flat,3,[],ns);
            for I = 1:size(M_W_flat,2)
                M_W(:,(1:3)+3*(I-1),:) = MultiProd_(M_W(:,(1:3)+3*(I-1),:),TD);
            end
            M_W_flat = reshape(M_W,size(M_W_flat));
        end
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function M_G = map_WtoG(M_W,R_G_W)
        M_G = MultiProd_(R_G_W,M_W);
        %M_G = mult_Anm1_Bmpz(R_G_W,M_W);
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function M_A = map_WtoA(M_W,R_A_W)
        M_A = mult_Anm1_Bmpz(R_A_W,M_W);
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function M_G = map_AtoG(M_A,R_G_A)
        M_G = mult_Anm1_Bmpz(R_G_A,M_A);
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function M_G = map_WE_G(M_WE,R_G_W,Rs_W_WE_tr,Rv_W_WE_tr,TD)
        M_W = map_WEtoW(M_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD);
        M_G = map_WtoG(M_W,R_G_W);
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function M_G_tr = map_GtoW(M_W_tr,R_G_W_tr)
        %M_A_tr = MultiProd_(M_W_tr,R_A_W_tr,[1 2]);
        %M_G_tr = MultiProd_(M_W_tr,R_G_W_tr);
        M_G_tr = mult_Anmz_Bmp1(M_W_tr,R_G_W_tr);
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    

%%


    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function int_val = intVal(x,y,int_fnc,dim)
        if nargin < 4
            dim = 3;
        end
        fullValOnly = true;
        int_val = int_fnc(x,y,fullValOnly,dim);
    end
    
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%
    function M_W2_cell = map_WEtoW_cell(M_WE_cell,Rs_W_WE_tr,Rv_W_WE_tr,TD)
        %M_W = MultiProd_(R_W_WE,MultiProd_(M_WE,R_W_WE_tr,[1 2]),[1 2]);%(3)x(3)x(ns)
        
        %if the rotation R_W_WE is modulo 90 degrees use following code
        %much quicker than MultiProd_ line
        smw = Rs_W_WE_tr; pmw = Rv_W_WE_tr;
        M_W2_cell = M_WE_cell;
        for ii_ = 1:3
            for jj_ = 1:3
                %M_W2_cell{ii_,jj_,:} = smw(ii_)*smw(jj_)*M_WE_cell{pmw(ii_),pmw(jj_),:};
                M_W2_cell(ii_,jj_) = {smw(ii_)*smw(jj_)*M_WE_cell{pmw(ii_),pmw(jj_)}};
                if ~isempty(TD)
                    M_WE_cell{ii_,jj_} = MultiProd_(M_W2_cell(ii_,jj_),TD);
                end
            end
        end
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function d3E_dzeta3_W_cell = get_d3E_dzeta3_W(ct,nszrs,ns,Rs_W_WE_tr,Rv_W_WE_tr,st_ss,st_cs,st_sp,st_cp,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp,TD)
    
        DDDe_WE = cell(3);

        DDDe_WE{1,1} = [%DDDex1_WE = [
            -ct_ss_sp,       -st_cs_sp,       -st_ss_cp;...
            -st_cs_sp,       -ct_ss_sp,        ct_cs_cp;...
            -st_ss_cp,        ct_cs_cp,       -ct_ss_sp;...
            -st_cs_sp,       -ct_ss_sp,        ct_cs_cp;...
            -ct_ss_sp, ss_cp-st_cs_sp, cs_sp-st_ss_cp;...
             ct_cs_cp, cs_sp-st_ss_cp, ss_cp-st_cs_sp;...
            -st_ss_cp,        ct_cs_cp,       -ct_ss_sp;...
             ct_cs_cp, cs_sp-st_ss_cp, ss_cp-st_cs_sp;...
            -ct_ss_sp, ss_cp-st_cs_sp, cs_sp-st_ss_cp];
        
        DDDe_WE{2,1} = [%DDDex2_WE = [
            -ct_cs_sp,        st_ss_sp,       -st_cs_cp;...
             st_ss_sp,       -ct_cs_sp,       -ct_ss_cp;...
            -st_cs_cp,       -ct_ss_cp,       -ct_cs_sp;...
             st_ss_sp,       -ct_cs_sp,       -ct_ss_cp;...
            -ct_cs_sp, cs_cp+st_ss_sp,-ss_sp-st_cs_cp;...
            -ct_ss_cp,-ss_sp-st_cs_cp, cs_cp+st_ss_sp;...
            -st_cs_cp,       -ct_ss_cp,       -ct_cs_sp;...
            -ct_ss_cp,-ss_sp-st_cs_cp, cs_cp+st_ss_sp;...
            -ct_cs_sp, cs_cp+st_ss_sp,-ss_sp-st_cs_cp];
        
        DDDe_WE{3,1} = [%DDDex3_WE = [%<~>faster to define blocks of zeros?
            -st_sp, nszrs, ct_cp;...
             zeros(1,3,ns)        ;...
             ct_cp, nszrs,-st_sp;...
             zeros(3,3,ns)        ;...
             ...
             ...
             ct_cp, nszrs,-st_sp;...
             zeros(1,3,ns)        ;...
            -st_sp, nszrs, ct_cp];
        
        DDDe_WE{1,2} = [%DDDey1_WE = [
             st_ss,-ct_cs, nszrs;...
            -ct_cs, st_ss, nszrs;...
             zeros(1,3,ns);
            -ct_cs, st_ss, nszrs;...
             st_ss,-ct_cs, nszrs;...
             zeros(4,3,ns)...
             ...
             ...
             ];
    
        DDDe_WE{2,2} = [%DDDey2_WE = [
             st_cs, ct_ss, nszrs;...
             ct_ss, st_cs, nszrs;...
             zeros(1,3,ns);...
             ct_ss, st_cs, nszrs;...
             st_cs, ct_ss, nszrs;...
             zeros(4,3,ns);...
             ...
             ...
             ];
        
        DDDe_WE{3,2} = [%DDDey3_WE = [
            -ct, zeros(1,2,ns);zeros(8,3,ns)];
        
        %DDDe_WE{1,3} = [];
        %DDDe_WE{2,3} = [];
        %DDDe_WE{3,3} = [];
        DDDe_WE{1,3} = [
             ct_ss_cp,        st_cs_cp,       -st_ss_sp;...
             st_cs_cp,        ct_ss_cp,        ct_cs_sp;...
            -st_ss_sp,        ct_cs_sp,        ct_ss_cp;...
             st_cs_cp,        ct_ss_cp,        ct_cs_sp;...
             ct_ss_cp, ss_sp+st_cs_cp,-cs_cp-st_ss_sp;...
             ct_cs_sp,-cs_cp-st_ss_sp, ss_sp+st_cs_cp;...
            -st_ss_sp,        ct_cs_sp,        ct_ss_cp;...
             ct_cs_sp,-cs_cp-st_ss_sp, ss_sp+st_cs_cp;...
             ct_ss_cp, ss_sp+st_cs_cp,-cs_cp-st_ss_sp];
        
        DDDe_WE{2,3} = [
             ct_cs_cp,       -st_ss_cp,       -st_cs_sp;...
            -st_ss_cp,        ct_cs_cp,       -ct_ss_sp;...
            -st_cs_sp,       -ct_ss_sp,        ct_cs_cp;...
            -st_ss_cp,        ct_cs_cp,       -ct_ss_sp;...
             ct_cs_cp, cs_sp-st_ss_cp, ss_cp-st_cs_sp;...
            -ct_ss_sp, ss_cp-st_cs_sp, cs_sp-st_ss_cp;...
            -st_cs_sp,       -ct_ss_sp,        ct_cs_cp;...
            -ct_ss_sp, ss_cp-st_cs_sp, cs_sp-st_ss_cp;...
             ct_cs_cp, cs_sp-st_ss_cp, ss_cp-st_cs_sp];
    
        DDDe_WE{3,3} = [
             st_cp, nszrs, ct_sp;...
             zeros(1,3,ns);...
             ct_sp, nszrs, st_cp;...
             zeros(3,3,ns);...
             ...
             ...
             ct_sp, nszrs, st_cp;...
             zeros(1,3,ns);...
             st_cp, nszrs, ct_sp];
        d3E_dzeta3_W_cell = map_WEtoW_cell(DDDe_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD);

    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function d2xi = get_d2xi(obj,st,ct,sp,cp,dth_ds,dsi_ds,dth_dt,dsi_dt,dph_dt,d2th_dsdt,d2si_dsdt,d2ph_dsdt,nszrs,yidx3,B,dB,BB,BdB,dBB,ns,R_W_WE_tr,st_sp,st_cp,ct_sp,ct_cp,yidx2)
        d2xi_dzeta2_IE = [%                                                  d2xi_dzeta2(:,1) = d2kappax_[dth[dth;dsi;dph] ; dsi[dth;dsi;dph] ; dph[dth;dsi;dph]]
            -dsi_ds.*ct_sp              , dsi_ds.*st , dsi_ds.*ct_cp              ;
             nszrs                       , nszrs      , nszrs                       ;
            -dsi_ds.*st_cp              , nszrs      ,-dsi_ds.*st_sp              ;
            %
             nszrs                       , nszrs      , nszrs                       ;
             nszrs                       , nszrs      , nszrs                       ;
             nszrs                       , nszrs      , nszrs                       ;
            %
            -dsi_ds.*st_cp              , nszrs      ,-dsi_ds.*st_sp              ;
             nszrs                       , nszrs      , nszrs                       ;
            -dsi_ds.*ct_sp - dth_ds.*cp , nszrs      , dsi_ds.*ct_cp - dth_ds.*sp];
        
        d2xi_dzetadzetaPrime_IE = [%-sp = d2[xi_1]/d[ph]d[thPrime]
             nszrs                       , nszrs      , nszrs                       ;
             nszrs                       , nszrs      , nszrs                       ;
            -sp                          , nszrs      , cp                          ;
            %
            -st_sp                      ,-ct         , st_cp                      ;
             nszrs                       , nszrs      , nszrs                       ;
             ct_cp                      , nszrs      , ct_sp                      ;
            %
             nszrs                       , nszrs      , nszrs                       ;
             nszrs                       , nszrs      , nszrs                       ;
             nszrs                       , nszrs      , nszrs];
         
        d2xi_dzetaPrimedzeta_IE = [
             nszrs                       , nszrs      , nszrs                       ;
            -st_sp                      ,-ct         , st_cp                      ;
             nszrs                       , nszrs      , nszrs                       ;
            %
             nszrs                       , nszrs      , nszrs                       ;
             nszrs                       , nszrs      , nszrs                       ;
             nszrs                       , nszrs      , nszrs                       ;
            %
            -sp                          , nszrs      , cp                          ;
             ct_cp                      , nszrs      , ct_sp                      ;
             nszrs                       , nszrs      , nszrs];
         
        d2xi_dzeta2 = mult_Anmz_Bmp1(d2xi_dzeta2_IE,R_W_WE_tr);
        d2xi_dzetadzetaPrime = mult_Anmz_Bmp1(d2xi_dzetadzetaPrime_IE,R_W_WE_tr);
        d2xi_dzetaPrimedzeta = mult_Anmz_Bmp1(d2xi_dzetaPrimedzeta_IE,R_W_WE_tr);
        %d2xi_dzetaPrimedzeta = multitransp(d2xi_dzetadzetaPrime);
        %------------------------------------------------------------------
        d2xi_dtdzeta_temp = bsxfun(@times,[dth_dt;dsi_dt;dph_dt;d2th_dsdt;d2si_dsdt;d2ph_dsdt],[reshape(d2xi_dzeta2,3,9,ns);reshape(d2xi_dzetaPrimedzeta,3,9,ns)]);
        d2xi_dtdzeta = reshape(sum(d2xi_dtdzeta_temp,1),3,3,ns);
        d2xi_dtdzetaPrime_temp = bsxfun(@times,[dth_dt;dsi_dt;dph_dt],reshape(d2xi_dzetadzetaPrime,3,9,ns));
        d2xi_dtdzetaPrime = reshape(sum(d2xi_dtdzetaPrime_temp,1),3,3,ns);
        
        temp10 = d2xi_dtdzeta(yidx2,:,:);
        temp20 = d2xi_dtdzetaPrime(yidx2,:,:);
        
        d2xi_dtdq = bsxfun(@times,B,temp10) + bsxfun(@times,dB,temp20);
        %------------------------------------------------------------------
        
        temp1 = d2xi_dzeta2(yidx3,:,:);
        temp2 = d2xi_dzetadzetaPrime(yidx3,:,:);
        temp3 = d2xi_dzetaPrimedzeta(yidx3,:,:);
        temp4 = bsxfun(@times,BdB,temp2);
        temp5 = bsxfun(@times,dBB,temp3);
       
        d2xi_dq2 = bsxfun(@times,BB,temp1) + temp4 + temp5;
        %------------------------------------------------------------------
        
        
        %dzeta_dt_ = [repmat(dth_ds,n*nq,1,1);repmat(dsi_ds,m*nq,1,1);repmat(dph_ds,o*nq,1,1)];
        %d2zeta_dsdt_ = [repmat(d2th_dsdt,n*nq,1,1);repmat(d2si_dsdt,m*nq,1,1);repmat(d2ph_dsdt,o*nq,1,1)];
        %dzetadtB = dzeta_dt_.*onsB;
        %d2zetadsdtB = d2zeta_dsdt_.*onsB;
        %dzetadtdB = dzeta_dt_.*onsdB;
       
        %d2xi_dtdq = bsxfun(@times,dzetadtB,temp1) + bsxfun(@times,d2zetadsdtB+dzetadtdB,temp2);
        %------------------------------------------------------------------
        
        %d2xi = {d2xi_dzeta2,d2xi_dzetadzetaPrime,d2xi_dq2,d2xi_dtdq};
        d2xi.d2xi_dzeta2 = d2xi_dzeta2;
        d2xi.d2xi_dzetadzetaPrime = d2xi_dzetadzetaPrime;
        d2xi.d2xi_dq2 = d2xi_dq2;
        d2xi.d2xi_dtdq = d2xi_dtdq;
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    function option_value = get_option(VARARGIN,option_flag,default_value)
        %checks for optional_flag amoung the input arguments in VARARGIN
        %if flag is present, return the following argument as its value
        %if flag is not present, return the default value
        [member,index] = find(strcmp(option_flag,VARARGIN));
        if member
            option_value = VARARGIN{index+1};
        else
            option_value = default_value;
        end        
    end
    
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

        
