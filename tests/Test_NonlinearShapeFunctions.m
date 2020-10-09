classdef Test_NonlinearShapeFunctions < matlab.unittest.TestCase
    %Test_NonlinearShapeFunctions Runs the test cases for the 
    
    properties (TestParameter)
        %List of files in ../examples/ that start with "example"
        ExampleBasedTest = getExampleScripts;
    end
    
    %   - examples
    methods (Test)
        function example_run_example_scripts(obj, ExampleBasedTest)
            %runExamples Attempts to run each of the example files in the
            %'../examples/' directory.
                        
            run(ExampleBasedTest);
            
        end
    end
end

%Generating parametric test inputs
function example_files = getExampleScripts
%getExampleScripts Returns a cell-array containing the full file path to
%the example scripts in the '\tbx\AwiWttToolbox\workflow' folder.

tbx_path    = fileparts(fileparts(mfilename('fullpath')));
example_loc = fullfile(tbx_path, '\examples');

%What is in the folder?
contents = dir(example_loc);
ext = cell(1, numel(contents));
for ii = 1 : numel(contents)
    [~, ~, ext{ii}] = fileparts(contents(ii).name);
end

%Only retain '.m' files
contents = contents(ismember(ext, '.m'));

%Only retain files that begin with 'example'
contents      = contents(startsWith({contents.name}, 'example'));
example_files = fullfile(example_loc, {contents.name});

end

