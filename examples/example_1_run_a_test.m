%% example_1_run_a_test
%
%
%

%Make a model
% ModelDefinition = parameter_sets.testCase_ClampedPatilHodgesWing;
ModelDefinition = parameter_sets.testCase_JA_preBendExample;

%Run an analysis
AnalysisOutput = runSim(0, 1, 'analysisType', 'static', 'fromObject', ModelDefinition);

%Post process the results
AnalysisOutput.generate_2dplot('JA_preBendExample', {'t',1,'1','1:nt'}, {'Gamma_G', 3, 'ns', '1:nt'});
