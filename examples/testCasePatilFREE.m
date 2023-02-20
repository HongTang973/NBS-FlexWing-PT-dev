function O = testCasePatilFREE(tipFORCE, qRigid)
%==========================================================================
%>>Master Level Object<<%
O = NBS_Master;
%--------------------------------------------------------------------------
%//////Flight Condition
O.V = 30; %m/s                                                             free stream airspeed
O.rho = 0.0881; %kg/m^3                                                    air density
O.aerodynamics = []; %                                         aerodynamic model
O.uVec_freeStream_G = [1;0;0]; %                                           free stream unit velocity vector
%--------------------------------------------------------------------------
O.grav_acc = 9.807*0;
O.gravVec_G = [0;0;-1];
%--------------------------------------------------------------------------
O.qRigidT = qRigid;
O.qRigidR = qRigid;
%--------------------------------------------------------------------------
O.plotBounds = [[-1 1];[-1 1];[-1 1]]*40;
%==========================================================================
%>>Add Flexible Patil/Hodges Starboard Wing

O_SBWing = parameter_sets.FlexPart_Library.flexPart_Patil_Hodges_16m_halfWing(O,O,'SBWing');
O_SBWing.alpha_root = 0; % the root angle of attack of the wing in degrees
O_SBWing.sweep_root = 0; % the root sweep angle of the wing in degrees
O_SBWing.wingRoot_offset_A = [0;0;0];

%//////Shape Functions
O_SBWing.shape_class_bend = 'chebyshev_1st';
O_SBWing.shape_class_twist = 'chebyshev_1st';
O_SBWing.shape_BCs_bend  = [0 1;1 1;1 1];
O_SBWing.shape_BCs_twist = [0 1;1 1;1 1];

O_SBWing.qth.n = 8;
O_SBWing.qsi.n = 6;
O_SBWing.qph.n = 6;

O_SBWing.tip_force_local = [0;0;tipFORCE];
O_SBWing.populate_shape_set('PLOT',false,'setName','Starboard Wing');

%==========================================================================
%>>Add Flexible Patil/Hodges Port Wing

O_PTWing = parameter_sets.FlexPart_Library.flexPart_Patil_Hodges_16m_halfWing(O,O,'PTWing');
O_PTWing.alpha_root = 0; O_PTWing.sweep_root = 0;
O_PTWing.wingRoot_offset_A = [0;0;0];

%//////Shape Functions
O_PTWing.shape_class_bend = 'chebyshev_1st';
O_PTWing.shape_class_twist = 'chebyshev_1st';
O_PTWing.shape_BCs_bend  = [0 1;1 1;1 1];
O_PTWing.shape_BCs_twist = [0 1;1 1;1 1];

O_PTWing.qth.n = 8;
O_PTWing.qsi.n = 6;
O_PTWing.qph.n = 6;

O_PTWing.tip_force_local = [0;0;tipFORCE];
O_PTWing.update_R_A_W('reflect',true);
O_PTWing.populate_shape_set('PLOT',false,'setName','Port Wing');

%==========================================================================
%>>Add Fuselage

% O_Fuselage = parameter_sets.RigidPart_Library.rigidPart_ExampleFuselage(O,O,'Fuselage');
% O_Fuselage.connection_idx_ParentObj = 1;
% O_Fuselage.connectionOffset_C = [0;0;0];
% O_Fuselage.R_C_W_0 = utility_functions.r_matrix([0;0;-1],pi/2);
% O_Fuselage.set_dependent_properties();

%==========================================================================
%>>Add HTP

% O_HTP = parameter_sets.RigidPart_Library.rigidPart_ExampleHTP(O,O,'HTP');
% O_HTP.connection_idx_ParentObj = 1;
% O_HTP.connectionOffset_C = [1;0;0]*O_Fuselage.s(end);
% O_HTP.R_C_W_0 = O_SBWing.R_A_W.'*utility_functions.r_matrix([0;1;0],0);
% O_HTP.R_C_W_0 = utility_functions.r_matrix([0;1;0],0);
% O_HTP.set_dependent_properties();

%==========================================================================
set_dependent_properties(O);
%============================ ==============================================
end