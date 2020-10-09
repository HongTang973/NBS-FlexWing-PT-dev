classdef QOI_GUI_Panel < handle % UI_Container
    
    properties %(Access = private)
        figureHandle
        QOI_Master
    end
    
    methods
        function obj = QOI_GUI_Panel(QOI_Master)
            obj.figureHandle = figure();
            obj.QOI_Master = QOI_Master;
        end
    end
    
    methods
        function delete(obj)
            for QOI_Container_ = obj.QOI_Master.QOIcontainers_cell
                delete(QOI_Container_{1}.TableView);
            end
            close(obj.figureHandle);
        end
    end
    
    
end

