function modify_sandbox_path(folders, option)
%modify_sandbox_path Add/remove sandbox folders from the MATLAB path.
%
% Original code taken from the Mathworks Toolbox Tools File Exchange code:
% https://www.mathworks.com/matlabcentral/fileexchange/60070-toolbox-tools
%
% Copyright 2016 The Mathworks, Inc.
%
% Edits by Christopher Szczyglowski, University of Bristol, 2020

assert(iscellstr(folders), ['Expected the folders to be provided ', ...
    'as a cell-array of strings.']); %#ok<ISCLSTR>
option = validatestring(option, {'add', 'remove'});

%Construct full file paths
package_directory       = fileparts(mfilename('fullpath'));
package_sub_directories = fullfile(package_directory, folders);

%Include all subdirectories in the 'tbx' folder.
idx_tbx = contains(package_sub_directories, 'tbx');
if nnz(idx_tbx) == 1
    tbx_folder = package_sub_directories{idx_tbx};
    contents = dir(tbx_folder);
    idx_directory = ...
        [contents.isdir] & ...
        ~strcmp({contents.name}, '.' ) & ...
        ~strcmp({contents.name}, '..' );
    contents = contents(idx_directory);
    tbx_sub_directories = strcat({contents.folder}, filesep, {contents.name})';
elseif nnz(idx_tbx) == 0
    tbx_sub_directories = [];
else
    error(['Ambigious match for the `tbx` directory. Expected there ', ...
        'to be only one folder in the package directory containing ', ...
        'the string `tbx`.']);
end
folder_path = [package_sub_directories ; tbx_sub_directories];

% Capture path
oldPathList = path();

% Add toolbox directory to saved path
userPathList = userpath();
if isempty( userPathList )
    userPathCell = cell( [0 1] );
else
    userPathCell = textscan( userPathList, '%s', 'Delimiter', ';' );
    userPathCell = userPathCell{:};
end
savedPathList = pathdef();
savedPathCell = textscan( savedPathList, '%s', 'Delimiter', ';' );
savedPathCell = savedPathCell{:};
savedPathCell = setdiff( savedPathCell, userPathCell, 'stable' );

switch option
    case 'add'
        
        savedPathCell = [folder_path; savedPathCell];
        path_fcn = @addpath;
        
    case 'remove'
        
        savedPathCell = setdiff( savedPathCell, folder_path, 'stable' );
        path_fcn = @rmpath;
        
end

path( sprintf( '%s;', userPathCell{:}, savedPathCell{:} ) )
savepath()

% Restore path plus toolbox directory
path( oldPathList )
path_fcn( sprintf( '%s;', folder_path{:} ) )

end

