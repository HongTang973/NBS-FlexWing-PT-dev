classdef QOI_Master < handle
    
    properties
        parentObject
        QOIcontainers_struct %a structure of qoi objects
        QOIcontainers_flat %flattened QOIcontainers_struct structure of qoi objects
        QOIcontainers_cell
        QOIcontainers_array
        qoi_request_systemLevel
        Tidx_request
    end
    
    properties
        nQOI_Containers
    end
    
    methods
        function obj = QOI_Master(parentObject)
            obj.parentObject = parentObject;
        end
        
        function add_QOI_Container(obj,address,partObj)
            QOI_Container = QOI_objects.QOI_Container(obj,address,partObj);
            set_field(obj,address,QOI_Container);
            
            address_flat = address;
            address_flat = strrep(address_flat,'.QOIcontainers_struct.','');
            address_flat = ['QOIcontainers_flat.' strrep(address_flat,'.','___')];
            set_field(obj,address_flat,QOI_Container);
            
            obj.QOIcontainers_cell = [obj.QOIcontainers_cell , {QOI_Container}];
            obj.QOIcontainers_array = [obj.QOIcontainers_array , QOI_Container];
        end
        
        function write_QOI_values(obj,qoiRequestSource,varargin)
            
            Display = get_option(varargin,'display',true);
            
            addpath('./utility_functions');
            addpath('./static_method_groups');
            addpath('./aerodynamic_codes');
            
            if Display, disp('Pushing request to qoi objects'); end
            obj.preProcess_qoi_request(qoiRequestSource);
                        
            if Display
                disp('Calculating requested values');
                waitBar = waitbar(0,'Processing Request');
            end
            
            Tidx_request_ = obj.Tidx_request;
            nTidx_request_ = numel(Tidx_request_);
            for i_ = 1:nTidx_request_
                tidx = Tidx_request_(i_);
                f(obj.parentObject.t(tidx),obj.parentObject.Q(:,tidx),obj.parentObject,'qoi',tidx);
                if Display, waitbar(i_/nTidx_request_,waitBar); end
            end
            
            obj.postProcess_qoi_request();
            if Display
                disp('done');            
                waitbar(1,waitBar);pause(0.1);close(waitBar);
            end
        end
        
        function preProcess_qoi_request(obj,qoiRequestSource)
            Tidx_eval = [];
            
            for i_ = 1:numel(obj.QOIcontainers_cell)
                QOI_Container_ = obj.QOIcontainers_cell{i_};
                QOI_Container_.pass_request_to_qois(qoiRequestSource);
                Tidx_eval = [Tidx_eval,QOI_Container_.Tidx_append_all];
            end
            
            obj.Tidx_request = unique(Tidx_eval);
        end
        
        function postProcess_qoi_request(obj)
            for i_ = 1:numel(obj.QOIcontainers_cell)
                QOI_Container_ = obj.QOIcontainers_cell{i_};
                QOI_Container_.postProcess_qoi_request();
            end
        end
    end
    
    
    methods
        function val = get.nQOI_Containers(obj)
            val = numel(obj.QOIcontainers_cell);
        end
    end
    
end
