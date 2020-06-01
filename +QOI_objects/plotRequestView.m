classdef plotRequestView < handle % UI_Container
    
    properties %(Access = private)
        QOI_Container
%         xdataListener
%         ydataListener
        xdata_selection_ui
        ydata_selection_ui
        xdata_dim1_ui
        ydata_dim1_ui
        xdata_sidx_ui
        ydata_sidx_ui
        xdata_tidx_ui
        ydata_tidx_ui
        %Listeners
        Command_cell
        generatingCode_ui
        NBS_Master_varName
    end
    
    properties
        qoiName_x
        qoiName_y
        qoiObject_x
        qoiObject_y
    end
    
    properties (Dependent)
        QOI_append_request
    end
    
%     events
%         EVENT_PlotRequestChange
%     end
    
    methods
        function obj = plotRequestView(UIparent,QOI_Container,NBS_Master_varName)
            
            plotRequestPanel = uipanel(UIparent,'Title','2d Plot');%,'FontSize',12,'BackgroundColor','white','Position',[20 20 440 321]);
            %uicontrol('Parent',uipanel,'Style','text','String','Select desired quantities of interest to write to object');
            
            
            obj.QOI_Container = QOI_Container; %<- QOI_Container
            
            obj.NBS_Master_varName = NBS_Master_varName;
            
            %tHBox = uix.HBoxFlex('Parent',plotRequestPanel);
            
            tVBox = uix.VBoxFlex('Parent',plotRequestPanel);
            tGrid = uix.Grid('Parent',tVBox);
            
            %tGrid = uix.Grid('Parent',plotRequestPanel);%,'Widths',[-1;-1],'Heights',[-1;-1]);
            uicontrol('Parent',tGrid,'Style','text');
            uicontrol('Parent',tGrid,'Style','text','String','xdata');
            uicontrol('Parent',tGrid,'Style','text','String','ydata');
            
            uicontrol('Parent',tGrid,'Style','text','String','QOI Name');
            obj.xdata_selection_ui = uicontrol('Parent',tGrid,'Style','popupmenu','String',obj.QOI_Container.qoiNames_all,'Callback',@obj.xSelectionChange);
            obj.ydata_selection_ui = uicontrol('Parent',tGrid,'Style','popupmenu','String',obj.QOI_Container.qoiNames_all,'Callback',@obj.ySelectionChange);
            
            uicontrol('Parent',tGrid,'Style','text','String','Component');
            obj.xdata_dim1_ui = uicontrol('Parent',tGrid,'Style','popupmenu','String',1,'Callback',@obj.update_generatingCode_ui);
            obj.ydata_dim1_ui = uicontrol('Parent',tGrid,'Style','popupmenu','String',1,'Callback',@obj.update_generatingCode_ui);
            
            uicontrol('Parent',tGrid,'Style','text','String','Spatial Indices');
            obj.xdata_sidx_ui = uicontrol('Parent',tGrid,'Style','edit','Callback',@obj.update_generatingCode_ui);
            obj.ydata_sidx_ui = uicontrol('Parent',tGrid,'Style','edit','Callback',@obj.update_generatingCode_ui);
            
            uicontrol('Parent',tGrid,'Style','text','String','Temporal Indices');
            obj.xdata_tidx_ui = uicontrol('Parent',tGrid,'Style','edit','Callback',@obj.update_generatingCode_ui);
            obj.ydata_tidx_ui = uicontrol('Parent',tGrid,'Style','edit','Callback',@obj.update_generatingCode_ui);
            
            uicontrol('Parent',tGrid,'Style','text');
            uicontrol('Parent',tGrid,'Style','pushbutton','String','Append Data','Callback',@obj.process_QOI_append_request);
            uicontrol('Parent',tGrid,'Style','pushbutton','String','Plot','Callback',@obj.gather_data_and_plot);
           %uicontrol('Parent',tGrid,'Style','pushbutton','String','Generate Code','Callback',@obj.update_generatingCode_ui);
            
            tGrid.Widths = [-0.5;-1;-1;-1;-1];
            tGrid.Heights = [22;22;22];
            tGrid.Padding = 10;
            tGrid.Spacing = 10;
            
            obj.generatingCode_ui = uicontrol('Parent',tVBox,'Style','edit','String','','Callback',@obj.update_generatingCode_ui);
            obj.generatingCode_ui.FontSize = 9;
            obj.generatingCode_ui.FontName = 'Monospaced';
            obj.generatingCode_ui.Visible = 'off';
            
            
            tVBox.Heights = [110;22];
            
            obj.xSelectionChange();
            obj.ySelectionChange();
            QOI_Container.plotRequestView = obj;
            
%             %create a listener
%             plotRequest_listener = event.listener(obj,'EVENT_PlotRequestChange',@obj.update_generatingCode_ui);
%             obj.Listeners = plotRequest_listener;
            
        end
    end
    
    methods
        
        function xSelectionChange(obj,~,~)
            idx = obj.xdata_selection_ui.Value;
            obj.qoiName_x = obj.QOI_Container.qoiNames_all{idx};
            obj.qoiObject_x = obj.QOI_Container.qoi_All_array(idx);
            obj.xdata_dim1_ui.String = 1:obj.QOI_Container.qoi_dim1_all{idx};
            obj.xdata_dim1_ui.Value = 1;
            if isequal(obj.qoiObject_x.Sidx_default_str,'-')
                obj.xdata_sidx_ui.String = '1';
                obj.xdata_sidx_ui.Enable = 'off';
            else
                obj.xdata_sidx_ui.String = '';
                obj.xdata_sidx_ui.Enable = 'on';
            end
            if isequal(obj.qoiObject_x.Tidx_default_str,'-')
                obj.xdata_tidx_ui.String = '1';
                obj.xdata_tidx_ui.Enable = 'off';
            else
                obj.xdata_tidx_ui.String = '';
                obj.xdata_tidx_ui.Enable = 'on';
            end
            obj.update_generatingCode_ui();
        end
        
        function ySelectionChange(obj,~,~)
            idx = obj.ydata_selection_ui.Value;
            obj.qoiName_y = obj.QOI_Container.qoiNames_all{idx};
            obj.qoiObject_y = obj.QOI_Container.qoi_All_array(idx);
            obj.ydata_dim1_ui.String = 1:obj.QOI_Container.qoi_dim1_all{idx};
            obj.ydata_dim1_ui.Value = 1;
            if isequal(obj.qoiObject_y.Sidx_default_str,'-')
                obj.ydata_sidx_ui.String = '1';
                obj.ydata_sidx_ui.Enable = 'off';
            else
                obj.ydata_sidx_ui.String = '';
                obj.ydata_sidx_ui.Enable = 'on';
            end
            if isequal(obj.qoiObject_y.Tidx_default_str,'-')
                obj.ydata_tidx_ui.String = '1';
                obj.ydata_tidx_ui.Enable = 'off';
            else
                obj.ydata_tidx_ui.String = '';
                obj.ydata_tidx_ui.Enable = 'on';
            end
            obj.update_generatingCode_ui();
        end
        
        function process_QOI_append_request(obj,~,~)
            obj.QOI_Container.QOI_Master.write_QOI_values('GUIplot');
        end
        
        function gather_data_and_plot(obj,~,~)
            
            obj.QOI_Container.generate_2dplot([],...
                {obj.qoiObject_x.qoiName;obj.qoiObject_y.qoiName},...
                {obj.xdata_dim1_ui.Value;obj.ydata_dim1_ui.Value},...
                {obj.xdata_sidx_ui.String,obj.ydata_sidx_ui.String},...
                {obj.xdata_tidx_ui.String,obj.ydata_tidx_ui.String}); 
        end
        
        function update_generatingCode_ui(obj,~,~)
            if isempty(obj.xdata_sidx_ui.String) || isempty(obj.xdata_tidx_ui.String) || isempty(obj.ydata_sidx_ui.String) || isempty(obj.ydata_tidx_ui.String)
                obj.generatingCode_ui.Enable = 'off';
            else
                obj.generatingCode_ui.Visible = 'on';
                obj.generatingCode_ui.Enable = 'on';
%             obj.generatingCode_ui.String = [...
%                 'Command: >> ' obj.NBS_Master_varName '.generate_2dplot(''' , obj.QOI_Container.containerName , '''' ,...
%                 ',{''' , obj.qoiName_x , '''',...
%                 ',', num2str(obj.xdata_dim1_ui.Value) ,...
%                 ',''' , obj.xdata_sidx_ui.String , '''',...
%                 ',''' , obj.xdata_tidx_ui.String , '''}',...
%                 ',{''' , obj.qoiName_y , '''',...
%                 ',' , num2str(obj.ydata_dim1_ui.Value) ,...
%                 ',''' , obj.ydata_sidx_ui.String , '''',...
%                 ',''' , obj.ydata_tidx_ui.String , '''});'];
            obj.Command_cell = {...
                    'Command: >> ',...
                    obj.NBS_Master_varName,...
                    '.generate_2dplot(',...
                    '''' , obj.QOI_Container.containerName , '''' ,...
                    ',{''' , obj.qoiName_x , '''',...
                    ',', obj.xdata_dim1_ui.Value ,...
                    ',''' , obj.xdata_sidx_ui.String , '''',...
                    ',''' , obj.xdata_tidx_ui.String , '''}',...
                    ',{''' , obj.qoiName_y , '''',...
                    ',' , obj.ydata_dim1_ui.Value ,...
                    ',''' , obj.ydata_sidx_ui.String , '''',...
                    ',''' , obj.ydata_tidx_ui.String , '''});'};
            obj.generatingCode_ui.String = cell2char(obj.Command_cell);
            
            end
        end
        
    end
    
    
    
    
    
    methods% (Access = private)
        
        
        function QOI_append_request = get.QOI_append_request(obj)
            
            %x data request
            SpatialIndices_x = obj.xdata_sidx_ui.String;
            TemporalIndices_x = obj.xdata_tidx_ui.String;
            
            %y data request
            SpatialIndices_y = obj.ydata_sidx_ui.String;
            TemporalIndices_y = obj.ydata_tidx_ui.String;
            
            QOI_append_request = {...
                    true,obj.qoiName_x,SpatialIndices_x,TemporalIndices_x;...
                    true,obj.qoiName_y,SpatialIndices_y,TemporalIndices_y};
            
        end

    end
    
    
end

