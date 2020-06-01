function [address,fileName] = splitFilePath(filePath)
    
    folder_sym_idx = strfind(filePath,filesep);
    
    if ~isempty(folder_sym_idx)
        folder_sym_last_idx = folder_sym_idx(end);
        address = filePath(1:folder_sym_last_idx-1);
        fileName = filePath(folder_sym_last_idx+1:end);
        
    else
        
        address = ['.' filesep];
        fileName = filePath;
        
    end

end