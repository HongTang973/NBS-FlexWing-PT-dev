classdef TableView < handle % UI_Container
    
    properties %(Access = private)
        QOI_Container
        Table
        Listeners
    end
    
    properties (Dependent)
        Data
        QOI_request
    end
    
    methods
        function obj = TableView(UIparent,QOI_Container)
            
            %validate inputs
            assert(isa(QOI_Container,'QOI_objects.QOI_Container'),'Not a valid application QOI_Container');
            
            %create graphics
            table = uitable(UIparent);
            %UIparent.TabTitles{end} = QOI_Container.containerName;
            
            %create listener
            %listen to changes in QOI_Container i.e. obj.QOI_Container
            %callback function will trigger obj.update()
            QOI_Container_listener = event.listener(QOI_Container,'EVENT_updateTable',@obj.onDataChanged);
            obj.Listeners = QOI_Container_listener;
            
            %store the QOI_Container
            obj.QOI_Container = QOI_Container; %<- QOI_Container
            obj.Table = table;
            QOI_Container.TableView = obj;
            %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>obj.Listeners = lh;
            
            %update graphics
            obj.update();
        end
    end
    
    methods
        function v = get.Data(obj)
            v = obj.Table.Data;
        end
        function v = get.QOI_request(obj)
            v = obj.Table.Data(:,1:4);
        end
    end
    
    methods% (Access = private)
        function onDataChanged(obj,source,event)
            obj.update();
        end
        
        function update(obj)
            
            QOI_Container = obj.QOI_Container;
            REQUESTED = ~QOI_Container.qoiFlat.is_default;    REQUESTEDedit = true;
            QOINAME = QOI_Container.qoiFlat.qoiName;          QOINAMEedit = false;
            SIDX = QOI_Container.qoiFlat.Sidx_str;            SIDXedit = true;
            TIDX = QOI_Container.qoiFlat.Tidx_str;            TIDXedit = true;
            DIMENSION = QOI_Container.qoiFlat.dimension_str;  DIMENSIONedit = false;
            
            obj.Table.Data =  [num2cell(REQUESTED)   ,QOINAME    ,SIDX    ,TIDX    ,DIMENSION];
            obj.Table.ColumnName = {'Requested','QOI Name','Spatial Indices','Temporal Indices','Dimension of QOI'};
            obj.Table.ColumnEditable = [REQUESTEDedit,QOINAMEedit,SIDXedit,TIDXedit,DIMENSIONedit];

        end
        
    end 
    
end

