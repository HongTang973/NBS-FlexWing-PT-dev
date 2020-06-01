classdef plotSpecialView < handle % UI_Container
    
    properties %(Access = private)
        QOI_Container
        plotType_ui
        fieldTitle_ui_cell
        fieldValue_ui_cell
        Command_cell
        generatingCode_ui
        NBS_Master_varName
    end
    
    properties
%         qoiName_x
%         qoiName_y
%         qoiObject_x
%         qoiObject_y
    end
    
    properties (Dependent)
%        QOI_append_request
    end
    
%     events
%         EVENT_PlotRequestChange
%     end
    
    methods
        function obj = plotSpecialView(UIparent,QOI_Container,NBS_Master_varName)
            
            plotRequestPanel = uipanel(UIparent,'Title','Plot Special');%,'FontSize',12,'BackgroundColor','white','Position',[20 20 440 321]);
                        
            obj.QOI_Container = QOI_Container; %<- QOI_Container
            
            obj.NBS_Master_varName = NBS_Master_varName;
            
            tVBox = uix.VBoxFlex('Parent',plotRequestPanel);
            tGrid = uix.Grid('Parent',tVBox);
            
            %tGrid = uix.Grid('Parent',plotRequestPanel);%,'Widths',[-1;-1],'Heights',[-1;-1]);
            uicontrol('Parent',tGrid,'Style','text','String','Plot Type');
            obj.plotType_ui = uicontrol('Parent',tGrid,'Style','popupmenu','String',{'-','3D Part Drawing'},'Callback',@obj.enableFields);
            
            for i_ = 1:4
                obj.fieldTitle_ui_cell{i_} = uicontrol('Parent',tGrid,'Style','text','String','');
                obj.fieldValue_ui_cell{i_} = uicontrol('Parent',tGrid,'Style','edit','String','','Enable','off','Callback',@obj.update_generatingCode_ui);
            end
            
            uicontrol('Parent',tGrid,'Style','text');
            uicontrol('Parent',tGrid,'Style','pushbutton','String','Plot','Callback',@obj.runCommand);
           %uicontrol('Parent',tGrid,'Style','pushbutton','String','button2','Callback',@obj.gather_data_and_plot);
           %uicontrol('Parent',tGrid,'Style','pushbutton','String','Generate Code','Callback',@obj.update_generatingCode_ui);
            
%            tGrid.Widths = [-0.5;-1;-1;-1;-1];
            tGrid.Heights = [22;22];
            tGrid.Padding = 10;
            tGrid.Spacing = 10;
            
            obj.generatingCode_ui = uicontrol('Parent',tVBox,'Style','edit','String','','Callback',@obj.update_generatingCode_ui);
            obj.generatingCode_ui.FontSize = 9;
            obj.generatingCode_ui.FontName = 'Monospaced';
            obj.generatingCode_ui.Visible = 'off';
            
            
            tVBox.Heights = [110;22];
            
%            obj.xSelectionChange();
%            obj.ySelectionChange();
%            QOI_Container.plotSpecialView = obj;
            
%             %create a listener
%             plotRequest_listener = event.listener(obj,'EVENT_PlotRequestChange',@obj.update_generatingCode_ui);
%             obj.Listeners = plotRequest_listener;
            
        end
    end
    
    methods
        
        function enableFields(obj,~,~)
            switch obj.plotType_ui.String{obj.plotType_ui.Value}
                case '-'
                    obj.clearFields;
                case '3D Part Drawing'
                    obj.clearFields;
                    obj.fieldTitle_ui_cell{1}.String = 'Part';
                    obj.fieldValue_ui_cell{1}.Style = 'popupmenu';
                    if strcmp(obj.QOI_Container.containerName,'aircraft')
                        obj.fieldValue_ui_cell{1}.String = {'all'};
                    else
                        obj.fieldValue_ui_cell{1}.String = {obj.QOI_Container.containerName,'all'};
                    end
                    obj.fieldValue_ui_cell{1}.Enable = 'on';
                    
                    obj.fieldTitle_ui_cell{2}.Style = 'popupmenu';
                    obj.fieldTitle_ui_cell{2}.BackgroundColor = obj.fieldTitle_ui_cell{2}.BackgroundColor + 0.01;
                    obj.fieldTitle_ui_cell{2}.String = {'Time Index','Time Value'};
                    obj.fieldValue_ui_cell{2}.Enable = 'on';
                    
                    obj.fieldTitle_ui_cell{3}.String = 'Plot Centre of Mass';
                    obj.fieldValue_ui_cell{3}.Style = 'popupmenu';
                    obj.fieldValue_ui_cell{3}.String = {'false','true'};
                    obj.fieldValue_ui_cell{3}.Enable = 'on';
            end
        end
        
        function clearFields(obj)
            for i_ = 1:numel(obj.fieldTitle_ui_cell)
                obj.fieldTitle_ui_cell{i_}.Style = 'text';
                obj.fieldTitle_ui_cell{i_}.String = '';
                obj.fieldValue_ui_cell{i_}.Style = 'edit';
                obj.fieldValue_ui_cell{i_}.String = '';
                obj.fieldValue_ui_cell{i_}.Enable = 'off';
            end
        end
        
        function update_generatingCode_ui(obj,~,~)
            switch obj.plotType_ui.String{obj.plotType_ui.Value}
                case '3D Part Drawing'
                    obj.update_generatingCode_plot3D;
            end
        end
        
        function update_generatingCode_plot3D(obj)
            
            invalidCommand = invalidCommandCondition(obj);
            if invalidCommand
                obj.generatingCode_ui.Enable = 'off';
            else
                obj.generatingCode_ui.Visible = 'on';
                obj.generatingCode_ui.Enable = 'on';
                obj.Command_cell = {...
                    'Command: >> ',...
                    obj.NBS_Master_varName,...
                    '.draw(',...
                    '''parts'',',...
                    '''' , obj.fieldValue_ui_cell{1}.String{obj.fieldValue_ui_cell{1}.Value} , ''''};
                    
                if strcmp(obj.fieldTitle_ui_cell{2}.String{obj.fieldTitle_ui_cell{2}.Value},'Time Index')
                    obj.Command_cell = [obj.Command_cell {',''Tidx'',' , obj.fieldValue_ui_cell{2}.String , ','} ];
                elseif strcmp(obj.fieldTitle_ui_cell{2}.String{obj.fieldTitle_ui_cell{2}.Value},'Time Value')
                    obj.Command_cell = [obj.Command_cell {',''t'',' , obj.fieldValue_ui_cell{2}.String , ','} ];
                end
                
                obj.Command_cell = [obj.Command_cell {'''plotCoM'',' , obj.fieldValue_ui_cell{3}.String{obj.fieldValue_ui_cell{3}.Value}}];
                obj.Command_cell = [obj.Command_cell {');'}];
                obj.generatingCode_ui.String = cell2char(obj.Command_cell);
            end
            
            function invalidCommand = invalidCommandCondition(obj)
                invalidCommand = isempty(obj.fieldValue_ui_cell{2}.String);
            end
        end
        
        function runCommand(obj,~,~)
            obj.update_generatingCode_plot3D;
            if strcmp(obj.generatingCode_ui.Enable,'on')
                f = figure;
                NBS_Master_obj = obj.QOI_Container.QOI_Master.parentObject; %#ok<NASGU>
                try
                    eval(['NBS_Master_obj' cell2char(obj.Command_cell(3:end))]);
                catch
                    msgbox('Invalid plotting command.');
                    pause(0.4);
                    close(f);
                end
            end
        end
        
    end
    
    
    
    
    
    
    
end

