classdef QOI_Container < handle
    %a mater object to contain individual qoi objects from a particular
    %aircraft part
    
    properties
        QOI_Master
        partObject
        qoiStruct
        qoi_cell = {}
        address

        Tidx_append_all
        
        %qoiRequest_unfiltered %requested qoi information regardless of current qoi data
        TableView
        plotRequestView
        discretisationVariables %<- variables that the user may utilise in the GUI, e.g. ns, nt, nsAp etc...
        
        qoi_request_partLevel %<- generic request at the aircraft part level
    end
    properties (Dependent)
        qoi_request_GUItable %<- property directly written by TableView object
        qoi_request_GUIplot %<- property directly written by plotRequestView object
        %GUI_table_data
        containerName
        qoiNames
        nqoi
        qoiFlat
        %qoiRequest_filtered %requested qoi information different to that already contained in qoi's
        %qoiSidx_ref
        %qoiDimension
        %qoiPopulated
    end
    
    properties%these will be hidden
        qoi_Generic %qois for s, t, s_aero etc.
        qoi_Generic_cell = {}
        qoi_All_cell
        qoi_All_array
    end
    properties (Dependent)
        qoiNames_Generic
        nqoi_Generic
        qoiFlat_Generic
        qoiNames_all
        qoi_dim1_all
    end
    
    events (NotifyAccess = protected)
        EVENT_updateTable
    end
    
    methods
        function obj = QOI_Container(parentObj,address,partObj)
            %TODO checks to implement: qoiAddress must be string
            %qoiAddress must point to a valid location
            obj.QOI_Master = parentObj;
            obj.address = address;
            obj.partObject = partObj;
            partObj.QOI_Container = obj;
            
            obj.add_generic_qoi('t',reshape(obj.partObject.t,1,1,[]),'-','1:nt','t','s');
            if ismember('s',properties(obj.partObject))
                obj.add_generic_qoi('s',reshape(obj.partObject.s,1,[],1),'1:ns','-','s','m');
            end
            if ismember('s_aero',properties(obj.partObject)) && ~isempty(obj.partObject.nsAp)
                obj.add_generic_qoi('sAp',reshape(obj.partObject.s_aeroMid,1,[],1),'1:nsAp','-','sAp','m');
                obj.discretisationVariables.nsAp = obj.partObject.nsAp;
            end
        end
        
        function add_generic_qoi(obj,qoiName,value,sidx_default,tidx_default,plotName_str,units_str)
            qoi_ = QOI_objects.qoi(obj,qoiName,sidx_default,size(value,1),plotName_str,units_str);
            qoi_.value = value;
            qoi_.Tidx_default_str = tidx_default;
            qoi_.Sidx = 1:size(value,2);
            qoi_.Tidx = 1:size(value,3);
            obj.qoi_Generic.(qoiName) = qoi_;
            obj.qoi_Generic_cell = [obj.qoi_Generic_cell;qoi_];
        end
        
        function add_qoi(obj,qoiName,tidx,value,dim2_str,plotName_str,units_str,varargin)
            if size(value,2) == 1, value = permute(value,[1 3 2]); end
            %add a qoi value to the current QOI_Container
            if isfield(obj.qoiStruct,qoiName) %if the qoi already exists then
                qoi_ = obj.qoiStruct.(qoiName);
                qoi_.addValue(tidx,value); %pass 'value' to the qoi (using qoi method 'addValue')
                
            else %if the qoi does not exist in the QOI_Container
                qoi_ = QOI_objects.qoi(obj,qoiName,dim2_str,size(value,1),plotName_str,units_str); %create the qoi object
                qoi_.address = [obj.address '.qoiStruct.' qoiName]; %ensure qoi.address is constent with its membership of QOI_Container
                obj.qoiStruct.(qoiName) = qoi_;
                obj.qoi_cell = [obj.qoi_cell;qoi_];
               %qoi_ = qoi_.appendValue(value); %append 'value' to the qoi (using qoi method 'appendValue')
            end
            
            isGlobalAeroQuantity = get_option(varargin,'GlobalAeroQuantity',false);
            
            if isGlobalAeroQuantity && isa(obj.partObject,'NBS_Master') %if the qoi is a global aerodynamic quantity and obj is the Master object
                NBS_Master_obj = obj.partObject;
                
                for pt_obj_cell = NBS_Master_obj.allParts_cell; %for each child object in the model
                    pt_obj = pt_obj_cell{1};
                    if ~isa(pt_obj,'NBS_Master') && pt_obj.isAero
                        sAp_idx_global = pt_obj.sAp_idx_global; %retrieve the global indices that reference the child object
                        pt_obj.QOI_Container.add_qoi(qoiName,tidx,value(:,sAp_idx_global),dim2_str,plotName_str,units_str); %extract the section of the global aero quantity corresponding to the object and add the section to its qoi list
                    end
                end
            end
            
        end
        
        function pass_request_to_qois(obj,qoiRequestSource,custom_qoiRequest,custom_mirror_append_flag)
            
            if nargin > 2
                assert(strcmp(qoiRequestSource,'custom'),'Warning: custom arguments are ignored if ''custom'' is not specified as qoiRequestSource');
            end
            
            %qoiRequest of the form {[true/false],'qoiName','sidx','tidx'}
            
            %retrieve qoiRequest based on qoiRequestSource argument
            %also specify whether to overwrite the current qoi data or append the request to the current qoi data
            switch qoiRequestSource
                case 'GUI'
                    qoiRequest = obj.qoi_request_GUItable; %get qoiRequest from GUI table
                    mirror_append_flag = 'mirror';
                case 'GUIplot'
                    qoiRequest = obj.qoi_request_GUIplot; %get qoiRequest from GUI 2d plot pane
                    mirror_append_flag = 'append';
                case 'default'
                    qoiRequest = [cell(numel(obj.qoiNames),1) obj.qoiNames]; %construct minimal qoiRequest input
                    mirror_append_flag = 'mirror';
                case 'partLevel'
                    qoiRequest = obj.qoi_request_partLevel; %fetch the qoiRequest that is already written to this QOI_Container
                    mirror_append_flag = 'append';
                case 'systemLevel'
                    qoiRequest = obj.QOI_Master.qoi_request_systemLevel; %fetch the qoiRequest that is already written to this QOI_Container
                    mirror_append_flag = 'append';
                case 'custom'
                    qoiRequest = custom_qoiRequest;
                    mirror_append_flag = custom_mirror_append_flag;
                %case 'systemLevel'
                %    qoiRequest = obj.QOI_Master.qoi_request_systemLevel;
            end
            
            if ~isempty(qoiRequest)
                %remove any qoiRequests for qois that are part of the generic set
                qoiRequestNames = qoiRequest(:,2);
                for qoiName_Generic_ = {obj.qoiNames_Generic}
                    idx = strcmp(qoiRequestNames,qoiName_Generic_);
                    qoiRequest(idx,:) = [];
                end
            end
            
            Tidx_eval = [];
            
            for i_ = 1:size(qoiRequest,1)
                request_ = qoiRequest(i_,:);
                if ~ismember(request_{2},obj.qoiNames)
                    continue;
                end
                qoi_ = obj.qoiStruct.(request_{2});
                if request_{1}
                qoi_.Sidx_request_str = request_{3};
                qoi_.Tidx_request_str = request_{4};
                else
                qoi_.Sidx_request_str = 'default';
                qoi_.Tidx_request_str = 'default';
                end
                    
                qoi_.preProcessRequest(mirror_append_flag);
                
                Tidx_eval = [Tidx_eval,reshape(qoi_.Tidx_append,1,[],1)];
            end
            
            obj.Tidx_append_all = unique(Tidx_eval);
        end
        
        function postProcess_qoi_request(obj)
            for i_ = 1:numel(obj.qoiNames)
                qoiName_ = obj.qoiNames{i_};
                qoi_ = obj.qoiStruct.(qoiName_);
                qoi_.postProcess_qoi_request();
            end
            obj.notify('EVENT_updateTable'); %once all qois in this QOI_Container have been updated, push event to update GUI table
        end
        
        function mat = StringToMat(obj,string)
            vars = obj.discretisationVariables;
            varNames = fields(vars);
            for j_ = 1:numel(varNames)
                varName_ = eval(['varNames{' num2str(j_) '}']);
                eval([varName_ '= vars.(varName_);']);
            end
            
            if isa(string,'double')
                mat = unique(string); %pass through argument if it is already numeric
            elseif strcmp(string,':')
                mat = ':';
            elseif isempty(string)
                mat = [];
            else
                mat = unique(eval(string));
            end
        end
        
        
        
        
        
        function [] = generate_2dplot(obj,axisHandle,qoiNames,components,sidx_str,tidx_str,varargin)
            
            %note: if axisHandle is blank, data will be drawn to a new figure
            
            while true
                
                idx_x = strcmp(obj.qoiNames_all,qoiNames(1));
                idx_y = strcmp(obj.qoiNames_all,qoiNames(2));
                
                qoiObject_x = obj.qoi_All_array(idx_x);
                qoiObject_y = obj.qoi_All_array(idx_y);
                
                [~,sidx_x] = ismember(obj.StringToMat(sidx_str{1}),qoiObject_x.Sidx); if strcmp(qoiObject_x.Sidx_default_str,'-'), sidx_x = 1; end
                [~,sidx_y] = ismember(obj.StringToMat(sidx_str{2}),qoiObject_y.Sidx); if strcmp(qoiObject_y.Sidx_default_str,'-'), sidx_y = 1; end
                [~,tidx_x] = ismember(obj.StringToMat(tidx_str{1}),qoiObject_x.Tidx); if strcmp(qoiObject_x.Tidx_default_str,'-'), tidx_x = 1; end
                [~,tidx_y] = ismember(obj.StringToMat(tidx_str{2}),qoiObject_y.Tidx); if strcmp(qoiObject_y.Tidx_default_str,'-'), tidx_y = 1; end
                
                if ismember(0,[sidx_x,tidx_x,sidx_y,tidx_y])
                    %this section should only trigger when the user is requesting a plot via the GUI
                    appendChoice = questdlg('Append plot data to current qoi set?', ...
                        'Append Data','Yes','No','No');
                    switch appendChoice
                        case 'Yes'
                            obj.plotRequestView.process_QOI_append_request();
                        case 'No'
                            return;
                    end
                    
                else
                    
                    break;
                    
                end
                
            end
            
            
            assert( numel(sidx_x)==1 || numel(tidx_x)==1 || numel(sidx_y)==1 || numel(tidx_y)==1 , 'Either x or y plot data must be a vector')
            assert( numel(sidx_x)==numel(sidx_y) || numel(tidx_x)==numel(tidx_y) , 'xdata and ydata must share an equal length dimension')
            
            if numel(sidx_x)==numel(sidx_y)
                xdata = permute(qoiObject_x.value(components{1},sidx_x,tidx_x),[2 3 1]);
                ydata = permute(qoiObject_y.value(components{2},sidx_y,tidx_y),[2 3 1]);
            end
            if numel(tidx_x)==numel(tidx_y)
                xdata = permute(qoiObject_x.value(components{1},sidx_x,tidx_x),[3 2 1]);
                ydata = permute(qoiObject_y.value(components{2},sidx_y,tidx_y),[3 2 1]);
            end
            
            if ~isempty(axisHandle)
                hsp = axisHandle;
            else
                figure(); hsp = subplot(1,1,1);
            end
            plot(hsp,xdata,ydata,varargin{:});
            xlabel_str = [qoiObject_x.plotName_str qoiObject_x.units_str];
            ylabel_str = [qoiObject_y.plotName_str qoiObject_y.units_str];
            if qoiObject_x.dim1>1
                xlabel_str = strrep(xlabel_str,'#',['(' num2str(components{1}) ')']);
            else
                xlabel_str = strrep(xlabel_str,'#','');
            end
            if qoiObject_y.dim1>1
                ylabel_str = strrep(ylabel_str,'#',['(' num2str(components{2}) ')']);
            else
                ylabel_str = strrep(ylabel_str,'#','');
            end
            xlabel(hsp,xlabel_str); ylabel(hsp,ylabel_str);
            title_str = {...
                ['xdata: ' qoiObject_x.qoiName ', Component: ' num2str(components{1}) ', Sidx: ' sidx_str{1} ', Tidx: ' tidx_str{1}];
                ['ydata: ' qoiObject_y.qoiName ', Component: ' num2str(components{2}) ', Sidx: ' sidx_str{2} ', Tidx: ' tidx_str{2}]};
            title(hsp,title_str,'Interpreter','none');

        end
        
        
        
        
        
        
        
        
    end
    
    methods
        function val = get.qoiNames(obj)
            val = fields(obj.qoiStruct);
        end
        function val = get.qoiNames_Generic(obj)
            val = fields(obj.qoi_Generic);
        end
        function val = get.nqoi(obj)
            val = numel(obj.qoiNames);
        end
        function val = get.nqoi_Generic(obj)
            val = numel(obj.qoiNames_Generic);
        end
        function qoi_All_cell = get.qoi_All_cell(obj)
            qoi_All_cell = {obj.qoi_All_array};
        end
        function qoi_All_array = get.qoi_All_array(obj)
            qoi_All_array = [obj.qoi_Generic_cell;obj.qoi_cell];
        end
        function qoiNames_all = get.qoiNames_all(obj)
            qoiNames_all = {obj.qoi_All_array.qoiName}.';
        end
        function qoi_dim1_all = get.qoi_dim1_all(obj)
            qoi_dim1_all = {obj.qoi_All_array.dim1}.';
        end
%         function val = get.qoiFlat(obj)
%             qoiAttributes = properties(obj.qoiStruct.(obj.qoiNames{1}));
%             CELL = cell(obj.nqoi,numel(qoiAttributes));
%             %CELL(1,:) = qoiAttributes.';
%             for i_ = 1:obj.nqoi
%                 qoi = obj.qoiStruct.(obj.qoiNames{i_});
%                 for j_ = 1:numel(qoiAttributes)
%                     CELL{i_,j_} = qoi.(qoiAttributes{j_});
%                 end
%             end
%             
%             TABLE = cell2table(CELL);
%             TABLE.Properties.VariableNames = qoiAttributes;
%             
%             val = TABLE;
%         end
        
        function val = get.qoiFlat(obj)
            val = obj.format_qoiFlat(obj.qoiStruct);
        end
        function val = get.qoiFlat_Generic(obj)
            val = obj.format_qoiFlat(obj.qoi_Generic);
        end
        function qoiFlat = format_qoiFlat(~,qoiStructure)
            qoiNames = fields(qoiStructure);
            nqoi = numel(qoiNames);
            qoiAttributes = properties(qoiStructure.(qoiNames{1}));
            CELL = cell(nqoi,numel(qoiAttributes));
            %CELL(1,:) = qoiAttributes.';
            for i_ = 1:nqoi
                qoi = qoiStructure.(qoiNames{i_});
                for j_ = 1:numel(qoiAttributes)
                    CELL{i_,j_} = qoi.(qoiAttributes{j_});
                end
            end
            
            TABLE = cell2table(CELL);
            TABLE.Properties.VariableNames = qoiAttributes;
            
            qoiFlat = TABLE;
        end
        function val = get.qoi_request_GUItable(obj)
            val = obj.TableView.QOI_request;
        end
        function val = get.qoi_request_GUIplot(obj)
            val = obj.plotRequestView.QOI_append_request;
        end
        function val = get.containerName(obj)
            idx = strfind(obj.address,'.');
            val = obj.address(idx(end)+1 : end);
        end
        
    end
    
end

