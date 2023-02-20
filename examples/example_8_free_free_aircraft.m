%%
qRigid = [];
initialDeflection = 1;

[x, fval, ~, ~] = fzero(@(x)tipDEF(x,qRigid,initialDeflection),0);
O_canti = testCasePatilFREE(x,qRigid);
O_static = runSim(0,1,'analysisType','static','fromObject', O_canti, 'suppressIter', true);
Q_predef = O_static.Q(:,end);

%%
qRigid = [0 0 0];
O_free = testCasePatilFREE(0,qRigid);
O_free.IC = [Q_predef; zeros(12,1)];
O_dynamic = runSim(0,30,'analysisType','dynamic','fromObject', O_free, 'delta_t', 0.01);
% O_dynamic.generate_video();

flex_part_name = O_dynamic.flexParts_nonlinear_cell{1}.partName;
gamma_dim3xnsxnt_tip = O_dynamic.get_qoiValue(flex_part_name,...
    'Gamma_G','Tidx','1:nt','Sidx','ns',...
    'generate_QOIs',true,...
    'display',false);
t = O_dynamic.t;

figure()
plot(t, squeeze(gamma_dim3xnsxnt_tip(3,:,:)), 'DisplayName', 'Tip')
hold on

gamma_dim3xnsxnt_root = O_dynamic.get_qoiValue(flex_part_name,...
    'Gamma_G','Tidx','1:nt','Sidx','1',...
    'generate_QOIs',true,...
    'display',false);
plot(t, squeeze(gamma_dim3xnsxnt_root(3,:,:)), 'DisplayName', 'Root')

flex_part_name = O_dynamic.flexParts_nonlinear_cell{1}.partName;
COM_G_SB = O_dynamic.get_qoiValue(flex_part_name,...
    'CoM_G','Tidx','1:nt','Sidx','1:3',...
    'generate_QOIs',true,...
    'display',false);
SB_mass = sum(O_dynamic.flexParts_nonlinear.(flex_part_name).ms);
flex_part_name = O_dynamic.flexParts_nonlinear_cell{2}.partName;
COM_G_PT = O_dynamic.get_qoiValue(flex_part_name,...
    'CoM_G','Tidx','1:nt','Sidx','1:3',...
    'generate_QOIs',true,...
    'display',false);

% rigidQ_T = O_dynamic.Q(O_dynamic.rT_idx,:);
% plot(t, rigidQ_T(3,:))
CoM_G = sum(bsxfun(@times,[repmat(SB_mass, [1 1 size(t,2)]) ; repmat(PT_mass, [1 1 size(t,2)])],[COM_G_SB; COM_G_PT]),1)/(PT_mass + SB_mass);
figure()
plot(t, squeeze(CoM_G(1,3,:)))

PT_mass = sum(O_dynamic.flexParts_nonlinear.(flex_part_name).ms);
CoM_info = [repmat(SB_mass, [1 size(t,2)]).', squeeze(COM_G_SB).'; repmat(PT_mass, [1 size(t,2)]).', squeeze(COM_G_PT).'];
aircraftMass = SB_mass + PT_mass;
CoM_G = sum(bsxfun(@times,CoM_info(:,1),CoM_info(:,2:4)),2)/aircraftMass;

plot(t, CoM_G, 'DisplayName', 'CoM_G')
legend('-DynamicLegend')

function residual = tipDEF(tipFORCE,qRigid,initialDeflection)
O = testCasePatilFREE(tipFORCE,qRigid);
O_static = runSim(0,1,'analysisType','static','fromObject', O, 'suppressIter', true);
flex_part_name = O_static.flexParts_nonlinear_cell{1}.partName;
gamma_dim3xnsxnt = O_static.get_qoiValue(flex_part_name,...
    'Gamma_G','Tidx','1:nt','Sidx','ns',...
    'generate_QOIs',true,...
    'display',false);

residual = gamma_dim3xnsxnt(3,1,2) - initialDeflection;
end
