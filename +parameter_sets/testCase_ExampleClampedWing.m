function O = testCase_ExampleClampedWing()

%==========================================================================
%>>Master Level Object<<%
O = NBS_Master;
%--------------------------------------------------------------------------
%//////Flight Condition
O.V = 30; %m/s                                                             free stream airspeed
O.rho = 0.0881; %kg/m^3                                                    air density
O.aerodynamics = 'strip_unsteady'; %                                         aerodynamic model
O.uVec_freeStream_G = [1;0;0]; %                                           free stream unit velocity vector
%--------------------------------------------------------------------------
O.grav_acc = 9.807;
O.gravVec_G = [0;0;-1];

O.plotBounds = [-0.5 0.5;0 1;-0.5 0.5]*16;

%==========================================================================
%>>Add Flexible Patil/Hodges Half-Wing

O_halfWing = parameter_sets.FlexPart_Library.flexPart_ExampleWing(O,O,'halfWing');
O_halfWing.alpha_root = 4; % the root angle of attack of the wing in degrees
O_halfWing.wingRoot_offset_A = [0;0;0];

%//////Applied Loads
O_halfWing.tip_force_global  = [0;0;0];
O_halfWing.tip_force_local   = [0;0;0];
O_halfWing.tip_moment_global = [0;0;0];
O_halfWing.tip_moment_local  = [0;0;0];
%//////Shape Functions
O_halfWing.shape_class_bend = 'chebyshev_1st';
O_halfWing.shape_class_twist = 'chebyshev_1st';
O_halfWing.shape_BCs_bend  = [0 1;1 1;1 1];
O_halfWing.shape_BCs_twist = [0 1;1 1;1 1];

O_halfWing.qth.n = 8;
O_halfWing.qsi.n = 6;
O_halfWing.qph.n = 6;

O_halfWing.populate_shape_set('PLOT',false);

%==========================================================================
set_dependent_properties(O);
%==========================================================================

end