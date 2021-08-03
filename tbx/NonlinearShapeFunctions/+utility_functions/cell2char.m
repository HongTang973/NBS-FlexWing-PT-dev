function charAry = cell2char(cell,varargin)

%concatenate a cell array of mixed character, string, and double entries into a single character array
%use 'parseChars' keyword to specify separating characters between concatenated entries
%
%e.g. cell2char({'t =',3,"seconds"},'parseChars',' ') -> 't = 3 seconds'

parseChars = utility_functions.get_option(varargin,'parseChars','');

if ~strcmp(class(cell),'cell'), cell = {cell}; end

charAry = '';
for i_ = 1:numel(cell)
    cell_entry = cell{i_};
    switch class(cell_entry)
        case 'char'
            charAry = [charAry parseChars cell_entry];
        case 'string'
            charAry = [charAry parseChars cell_entry{:}];
        case 'double'
            charAry = [charAry parseChars num2str(cell_entry)];
    end
end

if ~isempty(cell), charAry = charAry(numel(parseChars)+1 : end); end

end

