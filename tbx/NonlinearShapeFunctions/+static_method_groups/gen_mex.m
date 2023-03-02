clear functions
folder = fileparts(which('static_method_groups.gen_mex'));

%% ct methods
%map_WEtoW
M_WE = coder.typeof(0,[3,3,inf],[1 1 1]);  %3x3xns
Rs_W_WE_tr = coder.typeof(0,[3,1],[0 0]);  %3x1
Rv_W_WE_tr = coder.typeof(0,[3,1],[0 0]);  %3x1
TD = coder.typeof(0,[3,3],[0 0]);  %3x3

%map_WEtoW_stackDim2
M_WE_flat = coder.typeof(0,[3,3,inf],[1 1 1]);  %3x3xns
% Rs_W_WE_tr = coder.typeof(0,[3,1],[0 0]);  %3x1
% Rv_W_WE_tr = coder.typeof(0,[3,1],[0 0]);  %3x1
% TD = coder.typeof(0,[3,3],[0 0]);  %3x3


%% E methods
% get_E_group1(...
st = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ct = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
st_ss = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
st_cs = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
st_sp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
st_cp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ct_ss = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ct_cs = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ct_sp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ct_cp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ss_sp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ss_cp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
cs_sp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
cs_cp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
st_ss_sp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
st_ss_cp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
st_cs_sp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
st_cs_cp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ct_ss_sp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ct_ss_cp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ct_cs_sp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ct_cs_cp = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
% Rs_W_WE_tr = coder.typeof(0,[3,1],[0 0]);  %3x1
% Rv_W_WE_tr = coder.typeof(0,[3,1],[0 0]);  %3x1
% TD = coder.typeof(0,[3,3],[0 0]);  %3x3
nszrs = coder.typeof(0,[1,1,inf],[0 0 1]);  %1x1xns
ns = coder.typeof(0,1,0);  %ns
R_A_W = coder.typeof(0,[3,3], [0 0]);  %3x3
R_G_W = coder.typeof(0,[3,3], [0 0]);  %3x3
dzeta_a_dt_tr = coder.typeof(0,[1,inf,inf],[0 1 1]);  %1x3xns
Ba_tr = coder.typeof(0,[1, inf, inf],[0 1 1]);  %1xnaxns
yidx2 = coder.typeof(0,[1, inf],[0 1]);  %1x1xns
dR_G_W_dqe_3x3x1xnqe = coder.typeof(0,[3, 3, 1, inf],[0 0 0 1]);  %3x3x1xnqe
FLAG_free_free = coder.typeof(true,1,0);  %1

%% Gamma methods

%%
coder.extrinsic('mtimesx');
cfg = coder.config('mex');
cfg.TargetLang = 'C++';
cfg.InlineBetweenUserFunctions = 'Readability';
cfg.EnableAutoParallelization = 1;
cfg.IntegrityChecks = false;
cfg.ResponsivenessChecks = false;

%  cfg.DynamicMemoryAllocation = 'off';
% 
vars = {'-o' fullfile(folder,'mex') '-d' fullfile(folder,'codegen','mex','mex')...
    '-config' 'cfg' '-report'};

%% ct methods 
% vars = [vars(:)' {which('static_method_groups.map_WEtoW')} {'-args'} ...
%     {{M_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD}}];
% vars = [vars(:)' {which('static_method_groups.map_WEtoW_stackDim2')} {'-args'} ...
%     {{M_WE_flat,Rs_W_WE_tr,Rv_W_WE_tr,TD}}];

%% E methods
vars = [vars(:)' {which('static_method_groups.get_E_group1')} {'-args'} ...
    {{...
    ...
    st, ct,...
    st_ss,st_cs,st_sp,st_cp,...
    ct_ss,ct_cs,ct_sp,ct_cp,...
    ss_sp,ss_cp,cs_sp,cs_cp,...
    st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,...
    ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp,...
    ...
    Rs_W_WE_tr,Rv_W_WE_tr,TD,nszrs,ns,...
    R_A_W, R_G_W,...
    ...
    dzeta_a_dt_tr,Ba_tr,yidx2,...
    ...
    dR_G_W_dqe_3x3x1xnqe, FLAG_free_free}}];

%% Gamma methods

tic;
codegen(vars{:})
toc;