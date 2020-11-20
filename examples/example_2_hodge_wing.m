%% example_2_hodges_wing
%
%
%

%Make a model
ModelDefinition = parameter_sets.testCase_ClampedPatilHodgesWing;

%Run an analysis
AnalysisOutput = runSim(0, 1, 'analysisType', 'static', 'fromObject', ModelDefinition);

%Post process the results
part_name = AnalysisOutput.flexParts_nonlinear_cell{1}.partName;
AnalysisOutput.generate_2dplot(part_name, {'t',1,'1','1:nt'}, {'Gamma_G', 3, 'ns', '1:nt'});
