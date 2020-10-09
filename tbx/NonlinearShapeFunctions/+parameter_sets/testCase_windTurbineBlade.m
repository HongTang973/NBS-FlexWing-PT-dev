function O = testCase_windTurbineBlade()

%==========================================================================
%>>Master Level Object<<%
O = NBS_Master;
%--------------------------------------------------------------------------
%//////Flight Condition
O.V = 30; %m/s                                                             airspeed
O.rho = 0.0881; %kg/m^3                                                    air density
%O.aerodynamics = 'strip_steady';
O.uVec_freeStream_G = [1;0;0];
%--------------------------------------------------------------------------
O.grav_acc = 9.807*0;
O.gravVec_G = [0;0;-1];

%--------------------------------------------------------------------------
O.qRigidT = [];
O.qRigidR = [];
%==========================================================================

%==========================================================================
%>>Add Wind Turbine Blade

O_turbineBlade = parameter_sets.FlexPart_Library.flexPart_windTurbineBlade(O,O,'turbineBlade');
O_turbineBlade.wingRoot_offset_A = [0;0;0];

%//////Applied Loads
O_turbineBlade.tip_force_global  = [0;0;0];
O_turbineBlade.tip_force_local   = [0;0;0];
O_turbineBlade.tip_moment_global = [0;0;0];
O_turbineBlade.tip_moment_local  = [0;0;0];
%//////Shape Functions
O_turbineBlade.shape_class_bend = 'chebyshev_1st';
O_turbineBlade.shape_class_twist = 'chebyshev_1st';
O_turbineBlade.shape_class_shear = 'chebyshev_1st';
O_turbineBlade.shape_class_exten = 'chebyshev_1st';
O_turbineBlade.shape_BCs_bend  = [0 1;1 1;1 1];
O_turbineBlade.shape_BCs_twist = [0 1;1 1;1 1];
O_turbineBlade.shape_BCs_shear = [1 1;1 1;1 1];
O_turbineBlade.shape_BCs_exten = [1 1;1 1;1 1];

O_turbineBlade.qth.n = 4; O_turbineBlade.qth.group = 'qth1';
O_turbineBlade.qsi.n = 4; O_turbineBlade.qsi.group = 'qsi1';
O_turbineBlade.qph.n = 4; O_turbineBlade.qph.group = 'qph1';
O_turbineBlade.qSx.n = 0; O_turbineBlade.qSx.group = 'qSx1';
O_turbineBlade.qSy.n = 0; O_turbineBlade.qSy.group = 'qSy1';
O_turbineBlade.qSz.n = 0; O_turbineBlade.qSz.group = 'qSz1';

%set_dependent_properties(O_turbineBlade);
O_turbineBlade.populate_shape_set('PLOT',false);
%==========================================================================

O.plotBounds = [[-0.5 0.5];[0 1];[-0.5 0.5]]*1.2*O_turbineBlade.L;
set_dependent_properties(O);

end