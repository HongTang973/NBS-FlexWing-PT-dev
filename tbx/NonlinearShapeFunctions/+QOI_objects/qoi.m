classdef qoi < handle
    
    properties
        parentObject
        qoiName
        value
        address
        dim1
        plotName_str
        units_str = ''
    end
    properties (Dependent)
        is_empty
        is_default
        dimension
        dimension_str
        %sidx
        %tidx
    end
    properties
        Tidx
        Sidx
        Tidx_request, Tidx_request_str = 'default'                         % The initial temporal request string and corresponding indices
        Sidx_request, Sidx_request_str = 'default'                         % The initial spatial request string and corresponding indices
        Tidx_append                                                        % Temporal indices to append to the current index set in order to fulfil user request
        Sidx_append                                                        % Spatial indices to append to the current index set in order to fulfil user request
        
        Tidx_str
        Sidx_str
        Tidx_default_str = '[1 nt]'
        Sidx_default_str
        
        counter
    end
    
    methods %constructor method
        function obj = qoi(parentObj,qoiName,sidx_default,dim1,plotName_str,units_str)
            obj.parentObject = parentObj;
            obj.qoiName = qoiName;
            obj.Sidx_default_str = sidx_default;
            obj.dim1 = dim1;
            obj.plotName_str = plotName_str;
            if ~isempty(units_str)
                obj.units_str = [' (' units_str ')'];
            end
        end
    end
    
%     methods
%         function obj = appendValue(obj,value)
%             obj.value = cat(3,obj.value,value);
%         end
%     end
    
    methods% get/set methods
        function val = get.is_empty(obj)
            val = isempty(obj.value);
        end
        function val = get.is_default(obj)
            if isequal(obj.Tidx , obj.parentObject.StringToMat(obj.Tidx_default_str)) && isequal(obj.Sidx , obj.parentObject.StringToMat(obj.Sidx_default_str))
                val = true;
                obj.Tidx_str = obj.Tidx_default_str;
                obj.Sidx_str = obj.Sidx_default_str;
            else
                val = false;
            end
        end
        function val = get.dimension(obj)
            val = [...
                size(obj.value,1),...
                size(obj.value,2),...
                size(obj.value,3)];
        end
        function val = get.dimension_str(obj)
            dim = obj.dimension;
            val = [num2str(dim(1)) 'x' num2str(dim(2)) 'x' num2str(dim(3))];
        end
    end
    
    methods %methods for adding data to the qoi object
        function preProcessRequest(obj,mirror_append_flag)
            %set the various Tidx and Sidx quantities in response to the request strings Sidx_request_str & Tidx_request_str
            
            %if the requested temporal or spatial index string = 'default' then set it to its default string.
            if strcmp(obj.Sidx_request_str,'default'), obj.Sidx_request_str = obj.Sidx_default_str; obj.Sidx_str = obj.Sidx_default_str; end
            if strcmp(obj.Tidx_request_str,'default'), obj.Tidx_request_str = obj.Tidx_default_str; obj.Tidx_str = obj.Tidx_default_str; end
            
            obj.Tidx_request = obj.parentObject.StringToMat(obj.Tidx_request_str);
            obj.Sidx_request = obj.parentObject.StringToMat(obj.Sidx_request_str);
            
            if strcmp(mirror_append_flag,'mirror')
                [~, value_Tidx_delete] = setdiff(obj.Tidx,obj.Tidx_request); %find all of the existing Tidx that are not part of the Tidx_request
                [~, value_Sidx_delete] = setdiff(obj.Sidx,obj.Sidx_request); %find all of the existing Sidx that are not part of the Sidx_request
            elseif strcmp(mirror_append_flag,'append')
                %-----------------------------
                Sidx_append_request = obj.Sidx_request;
                if all(ismember(Sidx_append_request,obj.Sidx))
                    obj.Sidx_request_str = obj.Sidx_str;
                else
                    obj.Sidx_request_str = [ '[ ' obj.Sidx_str ' , ' obj.Sidx_request_str ' ]'];
                end
                obj.Sidx_request = obj.parentObject.StringToMat(obj.Sidx_request_str);
                %-----------------------------
                Tidx_append_request = obj.Tidx_request;
                if all(ismember(Tidx_append_request,obj.Tidx))
                    obj.Tidx_request_str = obj.Tidx_str;
                else
                    obj.Tidx_request_str = [ '[ ' obj.Tidx_str ' , ' obj.Tidx_request_str ' ]'];
                end
                obj.Tidx_request = obj.parentObject.StringToMat(obj.Tidx_request_str);
                %-----------------------------
                value_Tidx_delete = [];
                value_Sidx_delete = [];
            else
                error('mirror_append_flag argument to qoi.preProcessRequest method must be assigned a value in the set {''mirror''/''append''}');
            end
            obj.Tidx_append = setdiff(obj.Tidx_request,obj.Tidx); %find all of the Tidx_request that are not part of the existing Tidx
            obj.Sidx_append = setdiff(obj.Sidx_request,obj.Sidx); %find all of the Sidx_request that are not part of the existing Sidx
            
            if ~isempty(obj.Sidx_append) %if more s indices are requested than currently present then obj.value must be re-written in its entirety
                value_Tidx_delete = ':'; value_Sidx_delete = ':'; %i.e. delete all S and T entries in obj.value
                obj.Tidx_append = obj.Tidx_request; obj.Sidx_append = obj.Sidx_request; %Tidx_append and Sidx_append then become the entire requested set
            end
            
            obj.value(:,value_Sidx_delete,:) = [];
            obj.value(:,:,value_Tidx_delete) = [];
            obj.Tidx(value_Tidx_delete) = [];
            obj.Sidx = obj.Sidx_request;
            
            obj.counter = numel(obj.Tidx) + 1;
            
            append_value = zeros(...
                obj.dim1,...
                numel(obj.Sidx_request),...
                numel(obj.Tidx_append));
            
            if isempty(obj.value), obj.value = []; end
            if isempty(append_value), append_value = []; end
%             if isempty(obj.Tidx), obj.value = []; end
%             if isempty(obj.), obj.value = []; end
            
            obj.value = cat(3,obj.value,append_value);
            if ~isempty(obj.Tidx_append)
                obj.Tidx(numel(obj.Tidx)+numel(obj.Tidx_append)) = 0;
            end
            
        end
        
        function addValue(obj,tidx,value)
            if ismember(tidx,obj.Tidx_append)
                obj.value(:,:,obj.counter) = value(:,obj.Sidx_request);
                obj.Tidx(obj.counter) = tidx;
                obj.counter = obj.counter + 1;
            end
        end
        
        function postProcess_qoi_request(obj)
            [Tidx_srt,Tsrt_idx] = sort(obj.Tidx);
            [Sidx_srt,Ssrt_idx] = sort(obj.Sidx);
            obj.value = obj.value(:,Ssrt_idx,Tsrt_idx);
            obj.Tidx = Tidx_srt;
            obj.Sidx = Sidx_srt;
            obj.Tidx_str = obj.Tidx_request_str;
            obj.Sidx_str = obj.Sidx_request_str;
        end
    end
    
    methods
        function O = set.Tidx_request_str(O,str)
            if strcmp(str,':')
                str = O.Tidx_default_str;
            end
            O.Tidx_request_str = str;
        end
        function O = set.Sidx_request_str(O,str)
            if strcmp(str,':')
                str = O.Sidx_default_str;
            end
            O.Sidx_request_str = str;
        end
    end
    
end

