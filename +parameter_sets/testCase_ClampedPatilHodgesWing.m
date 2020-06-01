function O = testCase_ClampedPatilHodgesWing()

%==========================================================================
%>>Master Level Object<<%
O = NBS_Master;
%--------------------------------------------------------------------------
%//////Flight Condition
O.V = 30; %m/s                                                             airspeed
O.rho = 0.0881; %kg/m^3                                                    air density
O.aerodynamics = 'strip_steady';
O.uVec_freeStream_G = [1;0;0];
%--------------------------------------------------------------------------
O.grav_acc = 9.807*1;
O.gravVec_G = [0;0;-1];

%--------------------------------------------------------------------------
O.qRigidT = [];
O.qRigidR = [];

%//////Initial Condition
% O.ICStruct.qth.val = zeros(1,8);        O.ICStruct.dqth.val = [];
% O.ICStruct.qsi.val = zeros(1,8);        O.ICStruct.dqsi.val = [];
% O.ICStruct.qph.val = zeros(1,8);        O.ICStruct.dqph.val = [];
% O.ICStruct.qSx.val = zeros(1,0);        O.ICStruct.dqSx.val = [];
% O.ICStruct.qSy.val = zeros(1,0);        O.ICStruct.dqSy.val = [];
% O.ICStruct.qSz.val = zeros(1,0);        O.ICStruct.dqSz.val = [];
% O.ICStruct.qxT.val = zeros(1,0);        O.ICStruct.dqxT.val = [];
% O.ICStruct.qyT.val = zeros(1,0);        O.ICStruct.dqyT.val = [];
% O.ICStruct.qzT.val = zeros(1,0);        O.ICStruct.dqzT.val = [];
% O.ICStruct.qxR.val = zeros(1,0);        O.ICStruct.dqxR.val = [];
% O.ICStruct.qyR.val = zeros(1,0);        O.ICStruct.dqyR.val = [];
% O.ICStruct.qzR.val = zeros(1,0);        O.ICStruct.dqzR.val = [];
% 
% O.ICStruct.qAeroStates.val = [];



% qTheta          = zeros(1,8);        dqTheta          = 0*qTheta;
% qPsi            = zeros(1,8);        dqPsi            = 0*qPsi;
% qPhi            = zeros(1,8);        dqPhi            = 0*qPhi;
% qTau_x          = zeros(1,0);        dqTau_x          = 0*qTau_x;
% qTau_y          = zeros(1,0);        dqTau_y          = 0*qTau_y;
% qTau_z          = zeros(1,0);        dqTau_z          = 0*qTau_z;
% qRig_x          = [];                dqRig_x          = 0*qRig_x;
% qRig_y          = [];                dqRig_y          = 0*qRig_y;
% qRig_z          = [];                dqRig_z          = 0*qRig_z;
% qRig_varTheta_x = [];                dqRig_varTheta_x = 0*qRig_varTheta_x;
% qRig_varTheta_y = [];                dqRig_varTheta_y = 0*qRig_varTheta_y;
% qRig_varTheta_z = [];                dqRig_varTheta_z = 0*qRig_varTheta_z;
% 
% %--------------------------------
% O.IC = [...
%     qTheta  , qPsi    , qPhi   ,...
%     qTau_x  , qTau_y  , qTau_z ,...
%     qRig_x  , qRig_y  , qRig_z ,...
%     qRig_varTheta_x   , qRig_varTheta_y , qRig_varTheta_z ,...
%     ...
%     dqTheta , dqPsi   , dqPhi   ,...
%     dqTau_x , dqTau_y , dqTau_z ,...
%     dqRig_x , dqRig_y , dqRig_z ,...
%     dqRig_varTheta_x  , dqRig_varTheta_y , dqRig_varTheta_z...
%     ];
%--------------------------------

% nqth  = length(qTheta); nqsi  = length(qPsi);   nqph  = length(qPhi);
% nqTx = length(qTau_x); nqTy = length(qTau_y); nqTz = length(qTau_z);
% O.nmo = [nqth nqsi nqph nqTx nqTy nqTz];
%--------------------------------------------------------------------------

%==========================================================================

%==========================================================================
%>>Add Flexible Patil/Hodges Half-Wing

O_halfWing = parameter_sets.FlexPart_Library.flexPart_Patil_Hodges_16m_halfWing(O,O,'halfWing');
O_halfWing.wingRoot_offset_A = [0;0;0];

%//////Applied Loads
O_halfWing.tip_force_global  = [0;0;0];
O_halfWing.tip_force_local   = [0;0;0];
O_halfWing.tip_moment_global = [0;0;0];
O_halfWing.tip_moment_local  = [0;0;0];
%//////Shape Functions
O_halfWing.shape_class_bend = 'chebyshev_1st';
O_halfWing.shape_class_twist = 'chebyshev_1st';
O_halfWing.shape_class_shear = 'chebyshev_1st';
O_halfWing.shape_class_exten = 'chebyshev_1st';
O_halfWing.shape_BCs_bend  = [0 1;1 1;1 1];
O_halfWing.shape_BCs_twist = [0 1;1 1;1 1];
O_halfWing.shape_BCs_shear = [1 1;1 1;1 1];
O_halfWing.shape_BCs_exten = [1 1;1 1;1 1];

O_halfWing.qth.n = 4; O_halfWing.qth.group = 'qth1';
O_halfWing.qsi.n = 4; O_halfWing.qsi.group = 'qsi1';
O_halfWing.qph.n = 4; O_halfWing.qph.group = 'qph1';
O_halfWing.qSx.n = 0; O_halfWing.qSx.group = 'qSx1';
O_halfWing.qSy.n = 0; O_halfWing.qSy.group = 'qSy1';
O_halfWing.qSz.n = 0; O_halfWing.qSz.group = 'qSz1';

%set_dependent_properties(O_halfWing);
O_halfWing.populate_shape_set('PLOT',false);

% %------modify shape sets to allow for kink discontinuity-------------------
% sob = O_halfWing.shapeObject_bend;
% sot = O_halfWing.shapeObject_twist;
% 
% sob.addCustomFunctions(obj,s_custom,y_custom,dy_ds_custom,varargin)

%--------------------------------------------------------------------------

% O.flexParts_nonlinear.halfWing = O_halfWing;
% O.flexParts_nonlinear_cell = [O.flexParts_nonlinear_cell ; {O_halfWing}];
%==========================================================================

O.plotBounds = [[-8 8];[0 16];[-8 8]]*1;
set_dependent_properties(O);

end