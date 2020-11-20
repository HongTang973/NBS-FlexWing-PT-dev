classdef timeIndexConversion < handle % UI_Container
    
    properties
        NBS_Master
        Tidx_ui
        Tval_ui
    end
    
    methods
        function obj = timeIndexConversion(UIparent,NBS_Master)
            
            timeIndexConversionPanel = uipanel(UIparent,'Title','Time Index Conversion');%,'FontSize',12,'BackgroundColor','white','Position',[20 20 440 321]);
                        
            obj.NBS_Master = NBS_Master;
            
            tVBox = uix.VBoxFlex('Parent',timeIndexConversionPanel);
            tGrid = uix.Grid('Parent',tVBox);
            
            uicontrol('Parent',tGrid,'Style','text','String','Time Value');
            obj.Tval_ui = uicontrol('Parent',tGrid,'Style','edit','Callback',@obj.calculate_Tidx);
            
            uicontrol('Parent',tGrid,'Style','text','String','Closest Time Index');
            obj.Tidx_ui = uicontrol('Parent',tGrid,'Style','edit','Callback',@obj.calculate_t);
            
            
            tGrid.Heights = [12;22];
            tGrid.Padding = 10;
            tGrid.Spacing = 10;
            
        end
    end
    
    methods
        
        function calculate_Tidx(obj,~,~)
            O = obj.NBS_Master;
            t_request = str2num(obj.Tval_ui.String);
            [tidx,tval] = O.closestTidx(t_request);
            obj.Tidx_ui.String = num2str(tidx);
            obj.Tval_ui.String = num2str(tval);
        end
        
        function calculate_t(obj,~,~)
            O = obj.NBS_Master;
            tidx = str2num(obj.Tidx_ui.String);
            tval = O.t(tidx);
            obj.Tval_ui.String = num2str(tval);
        end
        
    end
        
end

