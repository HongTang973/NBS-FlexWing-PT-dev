function O = testCase_JA_preBendExample()

%==========================================================================
%>>Master Level Object<<%
O = NBS_Master;
%--------------------------------------------------------------------------
%//////Flight Condition
O.V = 30; %m/s                                                             airspeed
O.rho = 0.0881; %kg/m^3                                                    air density
% O.aerodynamics = 'strip_steady';
O.uVec_freeStream_G = [1;0;0];
%--------------------------------------------------------------------------
O.grav_acc = 9.807;
O.gravVec_G = [0;0;-1];
%--------------------------------------------------------------------------
O.qRigidT = [];
O.qRigidR = [0]; %#ok<NBRAK>
O.CUSTOM_free_states = @user_functions.custom_free_states;
%==========================================================================

%==========================================================================
%>>Add Blade JA_preBendExample

O_Blade = parameter_sets.FlexPart_Library.flexPart_JA_preBendExample(O,O,'JA_preBendExample');
O_Blade.wingRoot_offset_A = [0;1;0];

%//////Applied Loads
O_Blade.tip_force_global  = [0;0;0];
O_Blade.tip_force_local   = [0;0;0];
O_Blade.tip_moment_global = [0;0;0];
O_Blade.tip_moment_local  = [0;0;0];
%//////Shape Functions
O_Blade.shape_class_bend = 'chebyshev_1st';
O_Blade.shape_class_twist = 'chebyshev_1st';
O_Blade.shape_class_shear = 'chebyshev_1st';
O_Blade.shape_class_exten = 'chebyshev_1st';
O_Blade.shape_BCs_bend  = [0 1;1 1;1 1];
O_Blade.shape_BCs_twist = [0 1;1 1;1 1];
O_Blade.shape_BCs_shear = [1 1;1 1;1 1];
O_Blade.shape_BCs_exten = [1 1;1 1;1 1];

O_Blade.qth.n = 4; O_Blade.qth.group = 'qth1';
O_Blade.qsi.n = 4; O_Blade.qsi.group = 'qsi1';
O_Blade.qph.n = 4; O_Blade.qph.group = 'qph1';
O_Blade.qSx.n = 0; O_Blade.qSx.group = 'qSx1';  %TODO - Unit test the case where we have no shear/extension states
O_Blade.qSy.n = 0; O_Blade.qSy.group = 'qSy1';
O_Blade.qSz.n = 0; O_Blade.qSz.group = 'qSz1';

%set_dependent_properties(O_halfWing);
O_Blade.populate_shape_set('PLOT',false);

%==========================================================================

O.plotBounds = [[-8 8];[-16 16];[-16 8]]*1;
set_dependent_properties(O);

end