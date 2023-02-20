function SimObject = runSim(t1,t2,varargin)

if nargin == 0
    disp(' ')
    disp('function call of the form ''SimObject = runSim(t1,t2,options)''')
    disp('where t1 and t2 are the start and end times (ignored for static sim)')
    disp(' ')
    disp('additional keyword/value pairs specified in ''options'' argument')
    disp(' ')
    disp('optional keywords (default values marked with ''*''):')
    disp('        analysisType: [''dynamic''*/''static''] - specify dynamic or static analysis to be run')
    disp('        solver:       [''ode15s''*/''NewmarkBeta''] - choose between ode15s and Newmark Beta solvers (ignored for static sim)')
    disp('        intFnc:       [@functionName] - specify integration function to be used; Expected form, intVal = functionName(x,y)')
    disp('                                        where: size(x) = [1 nx], size(y) = [nfunc nx], size(intVal) = [nfunc nx]')
    disp('        profileCode:  [false*/true] - run code profiler if true')
    disp('        fromObject:   [object of type NBS_MasterObject] - inherit model and parameters from existing object (e.g. useful for re-starts)')
    return
end

%mtimesx('LOOPS');
% addpath('./utility_functions');
% addpath('./static_method_groups');
% addpath('./aerodynamic_codes');
%addpath(genpath('.')) <- add all subfolders to search path

%TODO - add generic inputParser;
clear f.m
analysisType = utility_functions.get_option(varargin,'analysisType','dynamic');
solver = utility_functions.get_option(varargin,'solver','ode15s');
intFnc = utility_functions.get_option(varargin,'intFnc',[]);
profileCode = utility_functions.get_option(varargin,'profileCode',false);
SimObj = utility_functions.get_option(varargin,'fromObject',[]);
if strcmpi(analysisType,'dynamic')
    displayWaitBar = utility_functions.get_option(varargin,'waitBar',true);
else
    displayWaitBar = false;
end
fHandle = utility_functions.get_option(varargin,'fHandle',@f);
testCase = utility_functions.get_option(varargin,'testCase',[]);

if strcmp(analysisType,'dynamic')
delta_t = utility_functions.get_option(varargin,'delta_t',min((t2-t1)/400,1/20));
    tsteps = t1:delta_t:t2; %time steps at which output is requested
    t_end = tsteps(end); t_temp = tsteps(1);
else
    displayWaitBar = false;
end

if profileCode
    profile on;
end

if displayWaitBar, waitBar = waitbar(0,'Initialising','CreateCancelBtn',@cancelFnc); end

if isempty(SimObj)
    if ~isempty(testCase)

        SimObject = parameter_sets.(['testCase_' testCase]);

    else

    SimObject = parameter_sets.testCase_ClampedPatilHodgesWing();

    end

else
    SimObject = SimObj;
end

if ~isempty(intFnc), SimObject.int_fnc = intFnc; end

if isequal(analysisType,'none')
    return;
end

if isequal(analysisType,'dynamic') %//dynamic analysis/////////////////////
SimObject.FLAG_static = false;
switch solver
    case 'NewmarkBeta'

tic,
[t,u] = Newmark_Beta_Solver_2ndOrder(tsteps,SimObject.IC.',SimObject);
runTime = toc; disp(['runTime: ' num2str(runTime)]);
SimObject.runTime = runTime;
SimObject.t = t.';
SimObject.Q = u;

    case 'ode15s'
                                                                           % %uncomment this code to return an estimate of the system jacobian
                                                                           % [J_est,err] = jacobianest(fHandle,SimObject.IC(:),{t1,SimObject.IC(:),SimObject,'dQ'});disp('done');J_est(abs(J_est)<1e-9)=0;
                                                                           % sort(imag(eig(J_est))/(2*pi)),
tic

options = odeset('OutputFcn',@outputFunction,'BDF','off','relTol',1e-5);
[t,u] = ode15s(fHandle, tsteps, SimObject.IC, options, SimObject,'dQ');

% options = odeset('OutputFcn',@outputFunction,'relTol',1e-5);
% [t,u] = ode45(fHandle, tsteps, SimObject.IC, options, SimObject,'dQ');

runTime = toc; disp(['   runTime: ' num2str(runTime)]);
SimObject.runTime = runTime;
SimObject.t = t.';
SimObject.Q = u.';

end

elseif isequal(analysisType,'static') %/////////////////////////////////static analysis/////////////////////
    SimObject.FLAG_static = true;
    options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter','MaxIter',1e3,'MaxFunctionEvaluations',1000,'OutputFcn',@getQ_iter);

    suppressIter = utility_functions.get_option(varargin,'suppressIter',false);
    if suppressIter, options.Display = 'none'; end

    Q_iter = [];

tic,

State_idx_static = [SimObject.qg2nd_idx(:) ; SimObject.qg1st_idx(:)];
x0 = SimObject.IC(State_idx_static);

[x,fval,exitFLAG,output] = fsolve(fHandle, x0, options,State_idx_static,SimObject,analysisType);

SimObject.Q = SimObject.IC.*[0,0];
SimObject.Q(State_idx_static,2) = x;

runTime = toc; %disp(['   runTime: ' num2str(runTime)]);
SimObject.runTime = runTime;
SimObject.t = [0 1];
SimObject.temp_properties.Q_iter = Q_iter;
SimObject.temp_properties.exitFLAG = exitFLAG;
SimObject.temp_properties.fval = fval;
SimObject.temp_properties.output = output.message;
end

if profileCode
profile off
SimObject.profileStructure = profile('info');
profile viewer
end




%% //////////////////////////////////////////////////////// Post Processing
if displayWaitBar, waitbar(1,waitBar,'Post Processing'); end

% profile on;
SimObject.initialise_QOI_Master();
SimObject.initialise_qois();
SimObject.QOI_Master.write_QOI_values('default','display',false);
profile off;
p = profile('info');

rootAddress = pwd;
for i_ = 1:numel(p.FunctionTable)
    fileAddress = p.FunctionTable(i_).FileName;
    if contains(fileAddress,rootAddress)
        fileAddress = strrep(fileAddress,pwd,'.');
        SimObject.archive(fileAddress);
    end
end

if displayWaitBar, close(waitBar); end
%//////////////////////////////////////////////////////////////////////////




%>>>>>>>>>>>>>>>>>>>>>>nested functions
function out = getQ_iter(Q,optimValues,state,state_idx,SimObject,outputFormat)

    if isequal(state,'iter')
        Q_iter = [Q_iter Q];
    end
    out = 0;

end

function out = outputFunction(t,Q,flag,~,~)
    waitBarSteps = 0.2; t_temp = -waitBarSteps;
    if displayWaitBar && ~isempty(t) && (t(1)-t_temp) >= waitBarSteps
        %disp(['t: ' num2str(t(1))]);
        waitbar(t(1)/t_end,waitBar,['Simulating System.  t: ' num2str(t(1))]);
        t_temp = t_temp + waitBarSteps;
    end

    out = 0;
end

end


function cancelFnc(~,~)
obj = gcbo;
if isprop(obj,'Style') %true if callback initiated from Cancel button
    delete(obj.Parent);
else                   %else callback initiated from X button
    delete(obj);
end
end
