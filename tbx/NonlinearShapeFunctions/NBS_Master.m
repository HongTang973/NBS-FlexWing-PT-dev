
classdef NBS_Master < handle

    properties%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
     %...properties..... .........................
      %..................derived_properties.......
       %|---------------|-------------------------

        %Global Properties <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        uVec_freeStream_G = [1;0;0]                                        %unit vector in free stream velocity direction, can provide function handle instead of value
        V = 0                                                              %[m/s] free steam velocity
        rho                                                                %[kg/m^3] air density
        grav_acc                                                           %[m/s^2] magnitude of gravitational acceleration
        gravVec_G = [0;0;-1]                                               %[-] global orientation of gravitational acceleration (unit vector)
        aerodynamics                                                       %[-] aerodynamics switch
        prescribedMotion_fnc                                               %handle to a function that prescribes an enforced motion of the aircraft reference point
        CUSTOM_free_states                                                 %handle to a function that prescribes a mapping from a set of kinematic states to position/rotational quantities and their variations
        ff_h
        R_G_A_0                                                            %initial orientation
        StateInfo
                        StateMap_cell
                        StateMap
                        IC                                                 %[-] initial condition vector
        qg2nd_idx                                                          %[-] indices of the second order states in the global 1st order state vector Q
        dqg2nd_idx                                                         %[-] indices of the second order state derivatives in the global 1st order state vector Q
        qg1st_idx                                                          %[-] indices of the fisrt order states in the global 1st order state vector Q
        qRigidT
        qRigidR
        dqRigidT = 0
        dqRigidR = 0
                        nqg2nd
        %               nqa                                                %total attitude shape degrees of freedom
        %               nqs                                                %total shear shape degrees of freedom
        %               nqf                                                %total flexible shape degrees of freedom
        %               nqr                                                %total rigid degrees of freedom
        %               nqAero                                             %total aerodynamic degrees of freedom
        %nmo                                                               %[-] [No._theta_functions, No._psi_functions, No._phi_functions, No._tau_z_functions, No._tau_x_functions, No._tau_y_functions]
        Mu = 0                                                             %[kg] point mass representing the aircraft rigid part
        I_Mu_A = zeros(3)                                                  %[kg m^2] point mass representing the aircraft rigid part
        rMu_A = [0;0;0]                                                    %[m] vector from the aircraft reference point to the rigid body centre of mass Mu

                        FLAG_free_free                                     %[true/false] true if rigid body states included in problem

        int_fnc = @utility_functions.integrate2_delxConst                                    %[function handle] integration function
        Gamma_int_fnc = @utility_functions.int_midpointShooting                              %[function handle] integration function for Gamma quantities
        plotBounds                                                         %[m] chosen plotBounds for 3D deflection plot
        %TODO recast fuselage, HTP and VTP as objects
        fuselage                                                           %structure containing information about the aircraft fuselage
        HTP                                                                %structure containing information about the horizontal tail plane
        VTP                                                                %structure containing information about the vertical tail plane

        %simulation results
        analysis_notes                                                     %any additional analysis notes
        t, Q                                                               %output time steps (t [s]) and corresponding state vector values (Q [-])
        Q_control
        %                nt                                                 %total number of time evaluation points
        QOI                                                                %additional quantities of interest
        QOI_index                                                          %a listing of all QOIs contained in the object and sub objects; object address and dimensions given for each QOI
        QOI_Master %experimental object-based QOI storage
        QOI_Container

        %temporary parameters of any data type
        %can be cleared at the end of the analysis
        FLAG_static
        runTime
        profileStructure
        object_creation_date
        fileArchive
        Camlight = [45,-90]
        temp_properties
        mainDir
        %flexParts_nonlinear_cell
        %               nflexParts_nonlinear
        flexParts_nonlinear
                        flexParts_nonlinear_fields
        rigidParts
                        nrigidParts
        rigidParts_cell
                        rigidParts_fields

        allParts_struct
        allParts_cell
                        nParts

        aeroPartNames = cell(0)

            rT_idx
            rR_idx
            drT_idx
            drR_idx

    end

    properties (SetAccess = {?NBS_Master,?NBS_flexPart_nonlinear,?NBS_rigidPart})
        flexParts_nonlinear_cell
                        nflexParts_nonlinear

                        nt
    end

    properties
        ModelName = 'NBS_Model_Version_8'
        partName = 'aircraft'
        mult3d_mex = (exist('mtimesx','file') == 3)                        %true if there exists a 'mtimesx' mex version of the AtimesB matrix operation on the matlab search path
    end

    %#ok<*MCSUP>

%==========================================================================
%//////////////////////////////////////////////////////////////////////////
%==========================================================================
    methods %constuctor and parameter methods

    function obj = NBS_Master(varargin)
%
%         package_dir = fileparts(mfilename('fullpath'));
%         addpath(fullfile(package_dir, 'utility_functions'));
%         addpath(fullfile(package_dir, 'static_method_groups'));
%         addpath(fullfile(package_dir, 'aerodynamic_codes'));
%         addpath(fullfile(package_dir, 'user_functions'));

        obj.object_creation_date = datetime;
        obj.ModelName = utility_functions.get_option(varargin,'Global',obj.ModelName);
        obj.partName = utility_functions.get_option(varargin,'partName',obj.partName);
%                                                                          obj.archive('f.m');
%                                                                          obj.archive('NBS_Master.m');
%                                                                          obj.archive('NBS_flexPart_nonlinear.m');
        obj.archive('ShapeFunctionObject.m');
        obj.archive('BeamPropertiesObject.m');
        obj.archive('+QOI_objects');
        obj.archive(['.' filesep 'utility_functions' filesep 'multiprod.m']);
%                                                                          obj.archive('utility_functions');
%                                                                          obj.archive('static_method_groups');

        obj.allParts_cell{end+1} = obj;
        obj.allParts_struct.(obj.partName) = obj;
    end
    
    function set_flex_properties(obj)
        obj.aeroPartNames = [];
        for partObj_ = obj.allParts_cell
            partObj = partObj_{1};
            if ~isa(partObj,'NBS_Master')
                partObj.set_dependent_properties();
            end
        end
    end
    
    function set_dependent_properties(obj)

        %call the set_dependent_properties methods of all aircraft sub parts
        obj.aeroPartNames = [];
        for partObj_ = obj.allParts_cell
            partObj = partObj_{1};
            if ~isa(partObj,'NBS_Master')
                partObj.set_dependent_properties();               
            end
        end

        %populate additional dependent aircraft parameters
        %if numel(O.nmo)==3, O.nmo(6) = 0; end %pad zeros for shear if required

        %--------------------
        obj.nflexParts_nonlinear = numel(obj.flexParts_nonlinear_cell);
        try
            obj.rigidParts_fields = fields(obj.rigidParts);
        catch
            obj.rigidParts_fields = {};
        end
        obj.nrigidParts = numel(obj.rigidParts_cell);
        obj.nParts = numel(obj.allParts_cell);
        %--------------------

        %------------------------------------------------------------------
        %Initialise StateInfo table structure and specify 1st order and 2nd order state groups
        obj.StateInfo = cell2table(cell(5,0));
        obj.StateInfo.Properties.RowNames = {'part','stateSet','nStates','Index','Initial Condition'};
        [~,stateSet_rowNum] = ismember('stateSet',obj.StateInfo.Properties.RowNames);
        [~,nStates_rowNum] = ismember('nStates',obj.StateInfo.Properties.RowNames);
        [~,Index_rowNum] = ismember('Index',obj.StateInfo.Properties.RowNames);
        StateSets2ndOrder = {'qth','qsi','qph','qSx','qSy','qSz','rigidBody','qCustom'};
        StateSets1stOrder = {'qAero'};
        StateSets = [StateSets2ndOrder , StateSets1stOrder];
        %----------------------------------
        %Populate the first row of StateInfo according to the 'n' and 'group' state information written to each system part
        flex_counter = 0;
        flexNames = cell(1,numel(obj.flexParts_nonlinear_cell));
        for ii_ = 1:obj.nParts
            PartObj = obj.allParts_cell{ii_};           
            if isa(PartObj,'NBS_flexPart_nonlinear')
                flex_counter = flex_counter  + 1;
                flexNames{flex_counter} = PartObj.partName;                
            end  
            
            for j_ = StateSets
                StateSet = j_{:};
                if isprop(PartObj,StateSet) && ~isempty(PartObj.(StateSet))
                    if isfield(obj.StateInfo,PartObj.(StateSet).group)
                        obj.StateInfo.(PartObj.(StateSet).group) = {...
                            PartObj.partName;...
                            StateSet;...
                            max(obj.StateInfo.(PartObj.(StateSet).group) , PartObj.(StateSet).n );...
                            inf;0};
                    else
                        obj.StateInfo.(PartObj.(StateSet).group) = {...
                            PartObj.partName;...
                            StateSet;...
                            PartObj.(StateSet).n;...
                            inf;0};
                    end

                    if ismember(StateSet,StateSets2ndOrder)
                        column = obj.StateInfo.(PartObj.(StateSet).group);
                        column{stateSet_rowNum} = ['d' column{stateSet_rowNum}];
                        obj.StateInfo.(['d' PartObj.(StateSet).group]) = column;
                    end
                end
            end
        end
        %----------------------------------
        %Populate the rigid state StateInfo row 1 entries
        obj.StateInfo.qRigidT = {'global';'rigidBody';numel(obj.qRigidT);inf;obj.qRigidT};
        obj.StateInfo.qRigidR = {'global';'rigidBody';0;inf;obj.qRigidR};
        obj.StateInfo.dqRigidT = {'global';'drigidBody';numel(obj.qRigidT);inf;obj.dqRigidT};
        obj.StateInfo.dqRigidR = {'global';'drigidBody';numel(obj.qRigidR);inf;obj.dqRigidR};
        %----------------------------------
        %Populate the row 2 StateInfo entries and get the state vector indices qg2nd_idx, dqg2nd_idx and qg1st_idx
        StateGroups = obj.StateInfo.Properties.VariableNames;
        idx_counter = 0;
        flexInd = ~strcmp(StateSets2ndOrder(1),StateGroups(1)); 
        flex_counter = reshape(repmat(1:numel(flexNames), 12, 1), 12*numel(flexNames), []);
        
        obj.qg1st_idx = []; obj.qg2nd_idx = []; obj.dqg2nd_idx = [];

        for j_ = 1:numel(StateGroups)
            StateGroup = StateGroups{j_};
            StateSet = obj.StateInfo.(StateGroup){stateSet_rowNum};
            nStatesInGroup = obj.StateInfo.(StateGroup){nStates_rowNum};
            
            if nStatesInGroup==0
                idx = [];
            else
                idx = (1:nStatesInGroup) + idx_counter;
                idx_counter = idx(end);
            end
            obj.StateInfo.(StateGroup){Index_rowNum} = idx;
            
            if any(strcmp(StateGroup,{'qRigidT','qRigidR'}))
                StateGroup_ = StateGroup(end);
                obj.(['r',StateGroup_,'_idx']) = idx;
            elseif any(strcmp(StateGroup,{'dqRigidT','dqRigidR'}))
                StateGroup_ = StateGroup(end);
                obj.(['dr',StateGroup_,'_idx']) = idx;
            end
            
            if ismember(StateSet,StateSets1stOrder)
                obj.qg1st_idx = [obj.qg1st_idx(:);obj.StateInfo.(StateGroup){Index_rowNum}(:)];
            elseif ismember(StateSet,StateSets2ndOrder)
                obj.qg2nd_idx = [obj.qg2nd_idx(:);obj.StateInfo.(StateGroup){Index_rowNum}(:)];
            elseif strcmp(StateSet(1),'d')
                obj.dqg2nd_idx = [obj.dqg2nd_idx(:);obj.StateInfo.(StateGroup){Index_rowNum}(:)];
            else
                error('This Should Not Be Reached')
            end
            obj.IC(idx) = obj.StateInfo{'Initial Condition',StateGroup}{:};
        end
        
        idx_counter = 0;
        for j_ = 1:numel(StateGroups)
            StateGroup = StateGroups{j_};
            StateSet = obj.StateInfo.(StateGroup){stateSet_rowNum};
            nStatesInGroup = obj.StateInfo.(StateGroup){nStates_rowNum};
            
            if nStatesInGroup==0
                idx = [];
                if any(strcmp(StateSet,{'qth','qsi','qph','qSx','qSy','qSz'}))
                    StateGroup_ = StateSet(2:end);
                    if flexInd
                        obj.flexParts_nonlinear.(flexNames{flex_counter(j_)}).([StateGroup_,'_idx']) = idx;
                    else
                        for i_ = 1:numel(flexNames)
                            obj.flexParts_nonlinear.(flexNames{i_}).([StateGroup_,'_idx']) = idx;
                        end
                    end
                elseif any(strcmp(StateSet,{'dqth','dqsi','dqph','dqSx','dqSy','dqSz'}))
                    StateGroup_ = [StateSet(1), StateSet(3:end)];
                    if flexInd
                        obj.flexParts_nonlinear.(flexNames{flex_counter(j_)}).([StateGroup_,'_idx']) = idx;
                    else
                        for i_ = 1:numel(flexNames)
                            obj.flexParts_nonlinear.(flexNames{i_}).([StateGroup_,'_idx']) = idx;
                        end
                    end
                end
            else
                idx = (1:nStatesInGroup) + idx_counter;
                if any(strcmp(StateSet,{'qth','qsi','qph','qSx','qSy','qSz'}))
                    StateGroup_ = StateSet(2:end);
                    if flexInd
                        obj.flexParts_nonlinear.(flexNames{flex_counter(j_)}).([StateGroup_,'_idx']) = idx;
                    else
                        for i_ = 1:numel(flexNames)
                            obj.flexParts_nonlinear.(flexNames{i_}).([StateGroup_,'_idx']) = idx;
                        end
                    end
                elseif any(strcmp(StateSet,{'dqth','dqsi','dqph','dqSx','dqSy','dqSz'}))
                    StateGroup_ = [StateSet(1), StateSet(3:end)];
                    if flexInd                    
                    obj.flexParts_nonlinear.(flexNames{flex_counter(j_)}).([StateGroup_,'_idx']) = idx;
                     else
                        for i_ = 1:numel(flexNames)
                            obj.flexParts_nonlinear.(flexNames{i_}).([StateGroup_,'_idx']) = idx;
                        end
                    end                   
                end
                idx_counter = idx(end);
            end
        end

        obj.IC = reshape(obj.IC,[],1);
        obj.nqg2nd = numel(obj.qg2nd_idx);

        %--------------------

        if obj.StateInfo.qRigidT{nStates_rowNum}==3, obj.FLAG_free_free = true; else, obj.FLAG_free_free = false; end
        
        obj.get_StateMap;
    end

    function initialise_QOI_Master(obj)

%         O.nflexParts_nonlinear = numel(O.flexParts_nonlinear_cell);
%         try
%             O.rigidParts_fields = fields(O.rigidParts);
%         catch
%             O.rigidParts_fields = {};
%         end
%         O.nrigidParts = numel(O.rigidParts_fields);

        QOI_Master_object = QOI_objects.QOI_Master(obj);

        QOI_Master_object.add_QOI_Container(['.QOIcontainers_struct.' obj.partName],obj);
        for i_flexPart = 1:obj.nflexParts_nonlinear
            flexPart_obj = obj.flexParts_nonlinear_cell{i_flexPart};
            QOI_Master_object.add_QOI_Container(['.QOIcontainers_struct.flexParts_nonlinear.' flexPart_obj.partName],flexPart_obj);
        end
        for i_rigidPart = 1:obj.nrigidParts
            rigidPart_obj = obj.rigidParts_cell{i_rigidPart};
            QOI_Master_object.add_QOI_Container(['.QOIcontainers_struct.rigidParts.' rigidPart_obj.partName],rigidPart_obj);
        end
        obj.QOI_Master = QOI_Master_object;

    end

    function initialise_qois(obj)
        f(0,reshape(obj.IC,[],1,1),obj,'qoi',1);
    end

    function populate_QOIs(obj,qoiRequest)

        if nargin == 2

            %generate quantities of interest from qoiRequest command if
            %supplied rather than generate GUI

            update_QOIs(obj,qoiRequest);

        else

            %construct GUI window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            s = ver('layout'); %does current matlab version include the GUI layout toolbox
            if isempty(s)
                error('Requires GUI Layout Toolbox');
            end


            %create a new figure to mount the qoi request panels
            GUIPANEL = QOI_objects.QOI_GUI_Panel(obj.QOI_Master);
            f1 = GUIPANEL.figureHandle;
            screenSize = get(0,'screensize');
            figureWidth = min(screenSize(3),760);
            figureHeight = min(screenSize(4),860);
            f1.OuterPosition = [...
                0.4*screenSize(3)-figureWidth/2,...
                40+(screenSize(4)-40)/2-figureHeight/2,...
                figureWidth,...
                figureHeight];

            %create vertical panel container
            fVbox = uix.VBoxFlex('Parent',f1,'Spacing',5);
            %create tap panel container
            tab = uix.TabPanel('Parent',fVbox);

            %get the names of all of the QOI containers
            QOIcontainer_names = fields(obj.QOI_Master.QOIcontainers_flat);
            %loop over the QOI containers
            nQOIcontainers = numel(QOIcontainer_names);

            for i_ = 1:nQOIcontainers
                QOIcontainer_name = QOIcontainer_names{i_};
                QOI_Container_ = obj.QOI_Master.QOIcontainers_flat.(QOIcontainer_name);

                %for each QOI container create a new tab panel
                tabVbox = uix.VBoxFlex('Parent',tab,'Spacing',5);
                tab.TabTitles{end} = QOI_Container_.containerName;

                %attach a table to the tab panel detailing the qois
                %available for request from the current QOI container
                Table = QOI_objects.TableView(tabVbox,QOI_Container_); %#ok<NASGU>

                PlotRequest = QOI_objects.plotRequestView(tabVbox,QOI_Container_,inputname(1)); %#ok<NASGU>

                PlotSpecial = QOI_objects.plotSpecialView(tabVbox,QOI_Container_,inputname(1)); %#ok<NASGU>

                tabVbox.Heights = [-2;-1;-1];

            end

            TimeIndexConversion = QOI_objects.timeIndexConversion(fVbox,obj); %#ok<NASGU>

            buttonHeight = 25;
            buttonWidth  = 100;
            fHButtonBox = uix.HButtonBox('Parent',fVbox,'ButtonSize',[buttonWidth,buttonHeight]);

            uicontrol('parent',fHButtonBox,'Style','pushbutton','String','Push Request','Callback',{@pushRequest,obj});
            uicontrol('parent',fHButtonBox,'Style','pushbutton','String','Close','Callback',{@closeGUI,GUIPANEL});

            fVbox.Heights(end-1) = 90;
            fVbox.Heights(end) = buttonHeight+10;

        end


        function pushRequest(~,~,O)
            O.QOI_Master.write_QOI_values('GUI');
        end
        function [] = closeGUI(~,~,GUIPANEL)
            delete(GUIPANEL);
        end

    end

    function generate_2dplot(obj,partName,qoiStructure_x,qoiStructure_y,varargin)

        if nargin == 1
            disp(' ')
            disp('function call of the form ''obj.generate_2dplot(partName,qoiStructure_x,qoiStructure_y,options)''')
            disp('where')
            disp('    qoiStructure: {qoiName,componentNumber,sidx,tidx}')
            disp('        qoiName [string]')
            disp('        componentNumber [+ve integer]')
            disp('        sidx [string] (can use the variable ''ns'' in this argument)')
            disp('        tidx [string] (can use the variable ''nt'' in this argument)')
            disp(' ')
            disp('additional keyword/value pairs specified in ''options'' argument')
            disp(' ')
            disp('optional keywords (default values marked with ''*''):')
            disp('        axisHandle: [@axisHandle/[]*] - specify axis handle to which plot objects will be drawn')
            disp('        plotOptions: [cell array] additional keyword value pairs accepted by Matlab ''plot()'' function')
            disp(' ')
            disp('example:')
            disp('O.generate_2dplot(''halfWing'',{''t'',1,''1'',''1:nt''},{''Gamma_G'',3,''ns'',''1:nt''},''axisHandle'',gca,''plotOptions'',{''LineWidth'',2});')
            disp(' ')
            return
        end

        axisHandle = utility_functions.get_option(varargin,'axisHandle',[]);
        plotOptions = utility_functions.get_option(varargin,'plotOptions',{});
        figureName = utility_functions.get_option(varargin,'figureName','');
        %example: O.generate_2dplot('halfWing',{'t',1,'1','1:nt'},{'Gamma_G',3,'ns','1:nt'})
        %qoiStructure: {qoiName,componentNumber,sidx,tidx}
        [qoiName_x,componentNumber_x,SpatialIndices_x,TemporalIndices_x] = qoiStructure_x{:};
        [qoiName_y,componentNumber_y,SpatialIndices_y,TemporalIndices_y] = qoiStructure_y{:};

        %generatedDataTreatment = get_option(varargin,'data','append');

        qoiRequest = {...
            true, qoiName_x,SpatialIndices_x,TemporalIndices_x;...
            true, qoiName_y,SpatialIndices_y,TemporalIndices_y};

        requestedPartObj = get_partObj(obj,partName);
        QOI_Container_requestedPart = requestedPartObj.QOI_Container;
        QOI_Container_requestedPart.qoi_request_partLevel = qoiRequest;


%         QOI_Container_requestedPart = [];
%
%         %find the QOI_Container for partName
%         %and write the qoiRequest to QOI_Container.qoi_request_partLevel
%         for qoiContainer = O.QOI_Master.QOIcontainers_array
%             if strcmp(qoiContainer.containerName,partName)
%                 QOI_Container_requestedPart = qoiContainer;
%                 QOI_Container_requestedPart.qoi_request_partLevel = qoiRequest;
%             end
%         end

        %run QOI_Master.write_QOI_values methods with partLevel request
        obj.QOI_Master.write_QOI_values('partLevel','display',false);

        %run the 2d plotting method contained in qoiContainer
        QOI_Container_requestedPart.generate_2dplot(axisHandle,...
            {qoiName_x;qoiName_y},...
            {componentNumber_x;componentNumber_y},...
            {SpatialIndices_x;SpatialIndices_y},...
            {TemporalIndices_x;TemporalIndices_y},figureName,plotOptions{:});

    end

    function archive(obj,item)

        itemType = exist(item); %#ok<EXIST>
        switch itemType
            case 2 %item is a .m file in the working directory
                fileName = item;
                archive_file(obj,fileName);
            case 6 %item is a .p file in the working directory
                fileName = item;
                archive_file(obj,fileName);
            case 7 %item is a sub-directory
                folderName = item;
                dirList = dir(folderName);
                for i_ = 1:length(dirList)
                    fileName = dirList(i_).name;
                    address = ['.' filesep folderName filesep];
                    if ~ismember(fileName,{'.','..'}) %skip '.' and '..' special entries
                        archive_file(obj,strcat(address,fileName));
                    end
                end
        end

        function archive_file(O,fileAddress)
            file_contents = {};
            [address_,fname_] = utility_functions.splitFilePath(fileAddress);
            fid = fopen(fileAddress,'r');
            while feof(fid) ~= 1
                fileLine = fgetl(fid);
                file_contents = [file_contents;{fileLine}]; %#ok<AGROW>
            end
            fclose(fid);
            fieldname = strrep(fname_,'.','_');
            O.fileArchive.(fieldname).fileContents = file_contents;
            O.fileArchive.(fieldname).archiveDate = datetime;
            O.fileArchive.(fieldname).address = [address_ filesep];
        end
    end

    function PartObject = get_partObj(obj,PartName)
        if ~isfield(obj.allParts_struct,PartName)
            disp('Input argument must be a valid aircraft part.');
            disp('Choose from:');
            disp(fields(obj.allParts_struct));
            error('Invalid Aircraft Part Name');
        end
        PartObject = obj.allParts_struct.(PartName);
    end

    function setParameterVal(obj,PartObject,PropertyName,Value)
        if ischar(PartObject), PartObject = obj.get_partObj(PartObject); end
        PartObject.(PropertyName) = Value;
    end

    function Value = getParameterVal(obj,PartObject,PropertyName)
        if ischar(PartObject), PartObject = obj.get_partObj(PartObject); end
        Value = PartObject.(PropertyName);
    end

    function val = get_qoiValue(obj,PartObject,qoiName,varargin)

        %returns the specified qoi value for the current flexPart
        %user may specify Sidx and Tidx indices to return

        Sidx = utility_functions.get_option(varargin,'Sidx',':');
        Tidx = utility_functions.get_option(varargin,'Tidx',':');
        generate_QOIs = utility_functions.get_option(varargin,'generate_QOIs',false);
        Display = utility_functions.get_option(varargin,'display',true);
        
        if ischar(PartObject), PartObject = obj.get_partObj(PartObject); end

        if generate_QOIs
            qoiRequest = {true, qoiName, Sidx, Tidx};
            PartObject.QOI_Container.qoi_request_partLevel = qoiRequest;
            obj.QOI_Master.write_QOI_values('partLevel', 'display', Display);
        end

        qoi = PartObject.QOI_Container.qoiStruct.(qoiName);

        if strcmp(Sidx,':')
            sIndex_data = ':';
        elseif strcmp(Sidx,'end')
            sIndex_data = numel(qoi.Sidx);
        else
            Sidx_mat = PartObject.QOI_Container.StringToMat(Sidx);
            [boolArray,sIndex_data] = ismember(Sidx_mat,qoi.Sidx);
            assert(all(boolArray),'Sidx argument must reference populated data entries.');
        end

        if strcmp(Tidx,':')
            tIndex_data = ':';
        else
            Tidx_mat = PartObject.QOI_Container.StringToMat(Tidx);
            [boolArray,tIndex_data] = ismember(Tidx_mat,qoi.Tidx);
            assert(all(boolArray),'Tidx argument must reference populated data entries.');
        end

        val = qoi.value(:,sIndex_data,tIndex_data);

    end

    function draw(obj,varargin)

        if strcmpi(varargin,'info')
            disp(' ')
            disp('Produce a 3d drawing of an aircraft part or parts')
            disp('method call of the form ''NBS_Master_Object.draw(options)''')
            disp(' ')
            disp('additional keyword/value pairs specified in ''options'' argument')
            disp(' ')
            disp('optional keywords (default values marked with ''*''):')
            disp('        parts: [''all*''/''partName''] - specification of aircraft parts to be plotted')
            disp('        plotCoM: [false*/true] - will plot the centre of mass of the selected parts')
            disp('        request_qoi_write: [true*/false] - if true the function will generate the required qoi data to fulfil the plotting request.')
            disp('        ...')
            return
        end
        
        request_qoi_write = utility_functions.get_option(varargin,'qoiRequest',true);
        plot_CoM = utility_functions.get_option(varargin,'plotCoM',false);
        parts = utility_functions.get_option(varargin,'parts','all');
        azimuth = utility_functions.get_option(varargin,'az',120);
        elevation = utility_functions.get_option(varargin,'el',20);
        AeroForce = utility_functions.get_option(varargin,'AeroForce',false);
        newFig = utility_functions.get_option(varargin,'newFig',true);
        figureName = utility_functions.get_option(varargin,'figureName','');
        if newFig
            figure('WindowStyle','docked','DockControls','on','Name',figureName);
        end
        axis equal;
        %AppliedForce = get_option(varargin,'AppliedForce',false);

        %view(ax,[120 20]);

        %---------------------------------------
        Tidxs = utility_functions.get_option(varargin,'Tidx',[]);
        t_requests = utility_functions.get_option(varargin,'t',[]);

        if isempty(Tidxs) && isempty(t_requests), Tidxs = obj.nt; end

        assert(isempty(Tidxs) || isempty(t_requests),'Must supply either one ''Tidx'' OR one ''t'' argument')
        %------------
        hold on;
        for i_ = 1:max(numel(Tidxs),numel(t_requests))
            if ~isempty(Tidxs), Tidx = Tidxs(i_); else, Tidx = []; end
            if ~isempty(t_requests), t_request = t_requests(i_); else, t_request = []; end

        %%Special Inputs
        %if strcmpi(Tidx,'end'), Tidx = numel(obj.t); end
        %if strcmpi(t_request,'end'), t_request = obj.t(end); end
        %------------
        tVec = obj.t;
        if ~isempty(Tidx) && isempty(t_request)
            t_ = tVec(Tidx);
        elseif ~isempty(Tidx) && ~isempty(t_request)
            warning(['''t'' argument to ' class(obj) ' ignored if ''Tidx'' specified']);
            t_ = tVec(Tidx);
        elseif isempty(Tidx) && isempty(t_request)
            Tidx = numel(tVec);
            t_ = tVec(Tidx);
        elseif isempty(Tidx) && ~isempty(t_request)
            [~,Tidx] = min(abs(tVec-t_request));
            t_ = tVec(Tidx);
        end
        %---------------------------------------

        %If requested, write the required qoi data for plotting.
        if request_qoi_write

            qoiRequest = [];

            if plot_CoM
                qoiRequest = [qoiRequest;{true,'CoM_G','1',num2str(Tidx)}]; %#ok<AGROW>
            end

            if AeroForce
                qoiRequest = [qoiRequest;{true,'Aero_Forces_G',':',num2str(Tidx)}]; %#ok<AGROW>
                qoiRequest = [qoiRequest;{true,'Aero_Forces_Gamma_G',':',num2str(Tidx)}]; %#ok<AGROW>
            end

            %specify a qoi Request for the required plotting fields
            qoiRequest = [ qoiRequest ; {...
                true,'ex_G','1:ns',num2str(Tidx);...
                true,'ey_G','1:ns',num2str(Tidx);...
                true,'ez_G','1:ns',num2str(Tidx);...
                true,'Gamma_G','1:ns',num2str(Tidx)} ]; %#ok<AGROW>

            %write qoi Request to the QOI_Container
            obj.QOI_Master.qoi_request_systemLevel = qoiRequest;

            %populate any extra info required based on above request
            obj.QOI_Master.write_QOI_values('systemLevel','display',false);
        end

        XLIM = obj.plotBounds(1,:); YLIM = obj.plotBounds(2,:); ZLIM = obj.plotBounds(3,:);

        %invoke the draw_part method for each aircraft part
        for partObj = obj.allParts_cell
            if ~isa(partObj,'NBS_Master') && (strcmp(parts,'all') || ismember({partObj{:}.partName},parts))
                partObj{:}.draw_part('Tidx',Tidx,...
                    'qoiRequest',false,...
                    'XLIM',XLIM,...
                    'YLIM',YLIM,...
                    'ZLIM',ZLIM,varargin{:});
            end
        end

        %Plot the centre of mass of the system if requested
        if plot_CoM
            CoM = obj.get_qoiValue(obj,'CoM_G','Tidx',Tidx,'Sidx',1,'generate_QOIs',false);
            CoM_plotStyle = {'marker','o','markerEdgeColor','none','markerFaceColor',[1 0 0]};
            plot3(CoM(1),CoM(2),CoM(3),CoM_plotStyle{:});
            plot3(XLIM(2),CoM(2),CoM(3),CoM_plotStyle{:}); plot3(CoM(1),YLIM(1),CoM(3),CoM_plotStyle{:}); plot3(CoM(1),CoM(2),ZLIM(1),CoM_plotStyle{:});
        end

        if AeroForce
            FAero = obj.get_qoiValue(obj,'Aero_Forces_G','Tidx',Tidx,'Sidx',':','generate_QOIs',false);
            FAero_Gamma = obj.get_qoiValue(obj,'Aero_Forces_Gamma_G','Tidx',Tidx,'Sidx',':','generate_QOIs',false);
            quiver3(FAero_Gamma(1,:),FAero_Gamma(2,:),FAero_Gamma(3,:),FAero(1,:),FAero(2,:),FAero(3,:));
        end

        end
        hold off;

        %set extra figure options
        ax = gca;
        view(ax,[azimuth,elevation]);
        box on; grid on;
%         title(['System Parts: ' utility_functions.cell2char(parts,'parseChars',', ') ' , t = ' num2str(t_)]);

    end

    function generate_video(obj,varargin)        
        %------------------------------------------------------------------
        framesPerSecond = utility_functions.get_option(varargin,'framerate',20);
        playSpeed = utility_functions.get_option(varargin,'playSpeed',1);
        tBounds = utility_functions.get_option(varargin,'tBounds',[0,obj.nt]);
        fileName = utility_functions.get_option(varargin,'fileName','video.avi');
        
        [tidx_start,~] = obj.closestTidx(tBounds(1));
        [tidx_end,~] = obj.closestTidx(tBounds(2));
        
        t_ = obj.t; nt_ = obj.nt;
        t0 = t_(1);
        for i_ = 2:nt_
            deltaT = (t_(i_)-t0);
            if deltaT > 1/framesPerSecond
                framesPerSecond_closestDiscreteFit = 1/deltaT;
                break
            end
        end
        Tidx = tidx_start:(i_-1):tidx_end;
        Tidx_str = ['[' num2str(Tidx) ']'];
        %------------------------------------------------------------------
        plot_CoM = utility_functions.get_option(varargin,'plotCoM',false);
        if plot_CoM
            qoiRequest = {true,'CoM_G','1',Tidx_str};
        else
            qoiRequest = [];
        end   
        
        %------------------------------------------------------------------
        
        %specify a qoi Request for the required plotting fields
        qoiRequest = [qoiRequest;{...
            true,'ex_G','1:ns',Tidx_str;...
            true,'ey_G','1:ns',Tidx_str;...
            true,'ez_G','1:ns',Tidx_str;...
            true,'R_C_W','1',Tidx_str;...
            true,'Gamma_G','1:ns',Tidx_str}];
        
        AeroForce = utility_functions.get_option(varargin,'AeroForce',false);
        if AeroForce
            qoiRequest = [qoiRequest;{true,'Aero_Forces_G',':',Tidx_str}];
            qoiRequest = [qoiRequest;{true,'Aero_Forces_Gamma_G',':',Tidx_str}]; 
        end    
        %write qoi Request to the QOI_Container
        obj.QOI_Master.qoi_request_systemLevel = qoiRequest;
        
        %populate any extra info required based on above request
        obj.QOI_Master.write_QOI_values('systemLevel','display',false);
        
        writerObj = VideoWriter(fileName);
        writerObj.FrameRate = framesPerSecond_closestDiscreteFit*playSpeed;
        open(writerObj);
        
        figure('WindowStyle','normal');
        set(gcf,'outerPosition',[100 100 800 800],'color',[1 1 1]);
        pause(0.1);
        
        for i_ = 1:numel(Tidx)
            tidx = Tidx(i_);
            cla
            obj.draw('parts','all','Tidx',tidx,'qoiRequest',false,'newFig',false,varargin{:})
            set(gcf,'Renderer','zbuffer');
            ax = gca;
            set(ax,'Units','pixels');
            view([90,-70])
%             set(ax, 'Xdir', 'reverse')
%             set(ax, 'Zdir', 'reverse')
            pos = get(ax,'Position');
            marg = 30;
%             rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
            title(['t = ' num2str(t_(tidx))]);
            F = getframe(gca);
            writeVideo(writerObj,F);
        end
        close(gcf);
    end

    function [Tidx,t] = closestTidx(obj,t_request)
        tColVec = reshape(obj.t,[],1);
        tRowVec_request = reshape(t_request,1,[]);
        [~,Tidx] = min(abs(tColVec-tRowVec_request));
        t = obj.t(Tidx);
    end

    function [] = plotShapes(obj,varargin)
        output_detail = utility_functions.get_option(varargin,'output_detail','final_1');
        flexParts_nonlinear_names = fields(obj.flexParts_nonlinear);
        for i_ = 1:obj.nflexParts_nonlinear
            flexPartName = flexParts_nonlinear_names{i_};
            shapeObject_bend = utility_functions.get_field(obj,['.flexParts_nonlinear.' flexPartName '.shapeObject_bend']);
            shapeObject_bend.plotShapes('output_detail',output_detail);
            shapeObject_twist = utility_functions.get_field(obj,['.flexParts_nonlinear.' flexPartName '.shapeObject_twist']);
            shapeObject_twist.plotShapes('output_detail',output_detail);
        end
    end

    function get_StateMap(obj)
        stateGroups_cell = obj.StateInfo.Properties.VariableNames;
        Part_cell = obj.StateInfo{'part',:};
        Index_cell = obj.StateInfo{'Index',:};
        nStateGroups = width(obj.StateInfo);
        val = cell(nStateGroups,5);
        for i_ = 1:nStateGroups
            stateGroup = stateGroups_cell{i_};
            part = Part_cell{i_};
            groupIndices = Index_cell{i_};
            for j_ = 1:numel(groupIndices)
                idx = groupIndices(j_);
                val(idx,:) = {['Component ',num2str(idx),' | ',part,' | ',stateGroup,' | DoF ',num2str(j_)] , num2str(idx) , part , stateGroup , ['DoF ' num2str(j_)]};
            end
        end
        obj.StateMap_cell = val(:,1);
        obj.StateMap = table;
        obj.StateMap.Component = val(:,2);
        obj.StateMap.Part = val(:,3);
        obj.StateMap.StateGroup = val(:,4);
        obj.StateMap.DoF = val(:,5);
    end

    function val = annoteState(obj,Q,varargin)

        idx = utility_functions.get_option(varargin,'indices',':');

        val = table(Q);
        val.Row = obj.StateMap_cell(idx,1);
    end

    function staticTrim(obj,partName_partProperty,partName_QOIName_QOIValue)

        nParameters = size(partName_partProperty,1);
        nQOIs = size(partName_QOIName_QOIValue,1);

        assert(nParameters == nQOIs, 'Number of free parameters must match the number of QOI conditions');
        nCon = nQOIs;

        Lambda_0 = zeros(nCon,1);
        Y_Target = zeros(nQOIs,1);

        for i_ = 1:nCon
            Lambda_0(i_) = obj.getParameterVal(partName_partProperty{i_,1},partName_partProperty{i_,2});
            Y_Target(i_) = partName_QOIName_QOIValue{i_,3};
        end

        warning('off','optim:fsolve:NonSquareSystem');
        Lambda = fsolve(@staticEval,Lambda_0,[],obj,Y_Target,nCon,partName_partProperty,partName_QOIName_QOIValue);
        warning('on','optim:fsolve:NonSquareSystem');

        %write the trimmed Lambda values to the object
        for i_ = 1:nCon
            obj.setParameterVal(partName_partProperty{i_,1},partName_partProperty{i_,2},Lambda(i_));
            obj.set_dependent_properties();
        end

        function Cost = staticEval(Lambda,obj_,Y_Target,nCon,partName_partProperty,partName_QOIName_QOIValue)
            %Lambda is the vector of free parameters

            %write these free parameters to ObjCopy
            for j_ = 1:nCon
                obj_.setParameterVal(partName_partProperty{j_,1},partName_partProperty{j_,2},Lambda(j_));
            end

            %run a static simulation of objCopy
            runSim(0,0,'analysisType','static','fromObject',obj_,'waitbar',false,'suppressIter',true);

            %use the output of the static simulation to write the output vector 'Y'
            Y = zeros(nCon,1);
            for k_ = 1:nCon
                Y(k_) = obj_.get_qoiValue(partName_QOIName_QOIValue{k_,1},partName_QOIName_QOIValue{k_,2},'Tidx','nt');
            end

            disp(['Property Value: ' num2str(Lambda)]);
            disp(['Output Value: ' num2str(Y)]);

            %assign a Cost based on the difference between 'Y' and 'Y_Target'
            Cost = Y - Y_Target;
        end

    end

    function simulate(obj,t1,t2,varargin)

        if nargin == 0
            disp(' ')
            disp('function call of the form ''obj = runSim(t1,t2,options)''')
            disp('where t1 and t2 are the start and end times (ignored for static obj)')
            disp(' ')
            disp('additional keyword/value pairs specified in ''options'' argument')
            disp(' ')
            disp('optional keywords (default values marked with ''*''):')
            disp('        analysisType: [''dynamic''*/''static''] - specify dynamic or static analysis to be run')
            disp('        solver:       [''ode15s''*/''NewmarkBeta''] - choose between ode15s and Newmark Beta solvers (ignored for static obj)')
            disp('        intFnc:       [@functionName] - specify integration function to be used; Expected form, intVal = functionName(x,y)')
            disp('                                        where: size(x) = [1 nx], size(y) = [nfunc nx], size(intVal) = [nfunc nx]')
            disp('        profileCode:  [false*/true] - run code profiler if true')
            disp('        fromObject:   [object of type NBS_MasterObject] - inherit model and parameters from existing object (e.g. useful for re-starts)')
            return
        end

        %mtimesx('LOOPS');
        %addpath(genpath('.')) <- add all subfolders to search path

        analysisType = utility_functions.get_option(varargin,'analysisType','dynamic');
        solver = utility_functions.get_option(varargin,'solver','ode15s');
        intFnc = utility_functions.get_option(varargin,'intFnc',[]);
        profileCode = utility_functions.get_option(varargin,'profileCode',false);
        displayWaitBar = utility_functions.get_option(varargin,'waitBar',true);
        fHandle = utility_functions.get_option(varargin,'fHandle',@f);

        if strcmp(analysisType,'dynamic')
            delta_t = utility_functions.get_option(varargin,'delta_t',min((t2-t1)/400,1/20));
            tsteps = t1:delta_t:t2; %time steps at which output is requested
            t_end = tsteps(end); t_temp = tsteps(1);
        else
            displayWaitBar = false;
        end

        if profileCode
            profile on;
        end

        if displayWaitBar, waitBar = waitbar(0,'Initialising','CreateCancelBtn',@cancelFnc); end

        if ~isempty(intFnc), obj.int_fnc = intFnc; end

        if isequal(analysisType,'none')
            return;
        end

        if isequal(analysisType,'dynamic') %//dynamic analysis/////////////////////

            switch solver
                case 'ode15s'

                    %f(0, obj.IC.', obj,'qoi');
                    tic

                    %JPattern = sparse([zeros(nqf) eye(nqf);ones(nqf,2*nqf)]);
                    %options = odeset('OutputFcn',@outputFunction,'JPattern',JPattern,'BDF','off','relTol',1e-5);

                    options = odeset('OutputFcn',@outputFunction,'BDF','off','relTol',1e-5);

                    % %uncomment this code to return an estimate of the system jacobian
                    % [J_est,err] = jacobianest(fHandle,obj.IC(:),{t1,obj.IC(:),obj,'dQ'});disp('done');J_est(abs(J_est)<1e-9)=0;
                    % sort(imag(eig(J_est))/(2*pi)),

                    %options = odeset('OutputFcn',@outputFunction,'BDF','off','relTol',1e-5,'Jacobian',@jacobian);
                    [t_,u] = ode15s(fHandle, tsteps, obj.IC, options, obj,'dQ');

                    runTime_ = toc;
                    obj.runTime = runTime_;
                    obj.t = t_.';
                    obj.Q = u.';
                    obj.nt = length(t_);
            end

        elseif isequal(analysisType,'static') %/////////////////////////////////static analysis/////////////////////

            suppressIter = utility_functions.get_option(varargin,'suppressIter',false);
            %maxEvals = 200;
            options = optimoptions('fsolve','Display','iter','MaxIter',1e3,'MaxFunctionEvaluations',10000,'OutputFcn',@getQ_iter);

            if suppressIter, options.Display = 'none'; end

            %if suppressIter, options = optimoptions('fsolve','MaxFunEvals',maxEvals,'OutputFcn',@getQ_iter);
            %else, options = optimoptions('fsolve','Display','iter','MaxIter',1e3,'MaxFunctionEvaluations',10000,'OutputFcn',@getQ_iter); end

            Q_iter = [];

            tic,

            %[x,fval,~,output] = fsolve(fHandle,obj.IC',options, obj);

            State_idx_static = [obj.qg2nd_idx(:) ; obj.qg1st_idx(:)];
            x0 = obj.IC(State_idx_static);

            [x,fval,~,output] = fsolve(fHandle, x0, options, State_idx_static,obj,analysisType);

            obj.Q = obj.IC*0;
            obj.Q(State_idx_static) = x;

            obj.t = 0;
            obj.temp_properties.Q_iter = Q_iter;
            runTime_ = toc;
            obj.runTime = runTime_;

        end

        if profileCode
            profile off
            obj.profileStructure = profile('info');
            profile viewer
        end




        %% //////////////////////////////////////////////////////// Post Processing
        if displayWaitBar, waitbar(1,waitBar,'Post Processing'); end

        profile on;
        obj.initialise_QOI_Master();
        obj.initialise_qois();
        obj.QOI_Master.write_QOI_values('default');
        profile off;
        p = profile('info');

        rootAddress = pwd;
        for i_ = 1:numel(p.FunctionTable)
            fileAddress = p.FunctionTable(i_).FileName;
            if contains(fileAddress,rootAddress)
                fileAddress = strrep(fileAddress,pwd,'.');
                obj.archive(fileAddress);
            end
        end

        if displayWaitBar, close(waitBar); end
        %//////////////////////////////////////////////////////////////////////////

        %>>>>>>>>>>>>>>>>>>>>>>nested functions
        function out = getQ_iter(Q,optimValues,state,state_idx,obj,outputFormat)

            if isequal(state,'iter')
                Q_iter = [Q_iter Q];
            end
            out = 0;

        end

        function out = outputFunction(t,Q,flag,~,~)
            waitBarSteps = 0.2; t_temp = -waitBarSteps;
            if displayWaitBar && ~isempty(t) && (t(1)-t_temp) >= waitBarSteps
                %disp(['t: ' num2str(t(1))]);
                waitbar(t(1)/t_end,waitBar,['Simulating System.  t: ' num2str(t(1))]);
                t_temp = t_temp + waitBarSteps;
            end

            out = 0;
        end

        function cancelFnc(~,~)
            obj_ = gcbo;
            if isprop(obj_,'Style') %true if callback initiated from Cancel button
                delete(obj_.Parent);
            else                   %else callback initiated from X button
                delete(obj_);
            end
        end

    end

    function draw_part(~,varargin)
    end

    end

%==========================================================================
%//////////////////////////////////////////////////////////////////////////
%==========================================================================

    methods (Static) %Utility functions performing useful operations

    function [] = dumpArchive(obj,fname_cell)
        if ischar(fname_cell)
            fname_cell = {fname_cell};
        end
        if nargin < 2
            fname_cell = fields(obj.fileArchive);
        end

        for ii = 1:length(fname_cell)
            fname = fname_cell{ii};
            append_str = '_archive_copy';
            dot_idx = strfind(fname,'_');
            fname_pt1 = fname(1:dot_idx(end)-1);
            fname_pt2 = fname(dot_idx(end)+1:end);
            fid = fopen([fname_pt1 append_str '.' fname_pt2],'w');
            fprintf(fid,'%s\n',obj.fileArchive.(fname).fileContents{:});
            fclose(fid);
        end
    end

    function prescribed_motion = prescribedMotion(obj,t)
        prescribed_motion = obj.prescribedMotion_fnc(obj,t);
    end

    function rmat = get_R_A_W(obj,wing)
        %update only the R_A_W rotational transform
        %useful if for example running multiple angle of attack cases for the same wing
        if strcmp(wing,'right'); siR = obj.sweep*pi/180; thR = atan(cos(siR)*tan(obj.dihedral*pi/180));
        elseif strcmp(wing,'left'); siR = -obj.sweep*pi/180; thR = -atan(cos(siR)*tan(obj.dihedral*pi/180));
        else, error('error in O.get_R_A_W');
        end
       %thR = atan(cos(siR)*tan(O.dihedral*pi/180));
       %phR = atan(cos(thR)/cos(siR)*(tan(O.alpha_0*pi/180)+tan(thR)*sin(siR)));
        phR = atan(1/cos(siR)*(tan(obj.alpha_0*pi/180)*cos(thR)+sin(siR)*sin(thR)))-abs(obj.alpha_0>90)*pi;
        eyRmat = [cos(thR).*sin(siR) ; cos(thR).*cos(siR) ; sin(thR)];
        exRmat = [cos(siR).*cos(phR) + sin(thR).*sin(siR).*sin(phR);- sin(siR).*cos(phR) + sin(thR).*cos(siR).*sin(phR);-cos(thR).*sin(phR)];
        ezRmat = [cos(siR).*sin(phR) - sin(thR).*sin(siR).*cos(phR);- sin(siR).*sin(phR) - sin(thR).*cos(siR).*cos(phR); cos(thR).*cos(phR)];

        rmat = [exRmat eyRmat ezRmat];
    end

    function Asample = sample(A,eval_idx,dim)
        %linearly sample matrix A at the indices eval_idx along the dimension 'dim'
        %example:
        %A = [1 2 3;2 4 6];
        %>> sample(A,[1.1 2.5 3],2)
        %    1.1    2.5    3.0
        %    2.2    5.0    6.0

        %---- extract the integer part of the evaluation indices
        eval_idx_floor  = floor(eval_idx);
        indices_fl  = {':',':',':'}; indices_fl{dim}  = eval_idx_floor;
        %---- modified integer part to handle special case eval_idx(i_) = size(A,dim)
        eval_idx_floor_ = min(eval_idx_floor,size(A,dim)-1);
        indices_fl_ = indices_fl;    indices_fl_{dim} = eval_idx_floor_;
        %---- extract the decimal part of the evaluation indices
        indices_dc = eval_idx - eval_idx_floor;
        dimensions_indices_dc = [1,1,1]; dimensions_indices_dc(dim) = length(eval_idx);

        indices_dc_rs = reshape(indices_dc,dimensions_indices_dc);

        Andim = length(size(A));
        Adiff = diff(A,1,dim);

        %linear sampling of A matrix
        Asample = A(indices_fl{1:Andim}) + bsxfun(@times,Adiff(indices_fl_{1:Andim}),indices_dc_rs);
    end

    function AB = mult3d(A,B)
        %quick performance for matrices where dim1/dim2 << dim3
        dimL = size(A,1); dim_inner = size(A,2); dimR = size(B,2);
        AB = zeros(dimL,dimR,max(size(A,3),size(B,3)));
        for ii_ = 1:dimL
            for jj_ = 1:dimR
                AB_ = 0;
                for kk_ = 1:dim_inner
                    AB_ = AB_ + A(ii_,kk_,:).*B(kk_,jj_,:);
                end
                AB(ii_,jj_,:) = AB_;
            end
        end
    end

    function AB = mult_Anm1_Bmpz(A,B)
        %multiplication of the 2D matrix A by each 2D slice (:,:,i) of matrix B
        szA = size(A);
        szB = size(B);
        AB = reshape(A*reshape(B,szB(1),szB(2)*szB(3),1),[szA(1) szB(2) szB(3)]);
    end

    function AB = mult_A11z_Bnmz(A,B)
        %element-wise multiplication of the 1x1xz vector A by each 1D strip (i,j,:) of matrix B
        AB = bsxfun(@times,A,B);
    end

    function AB = mult_Anmz_Bmp1(A,B)
        %multiplication of each 2D slice (:,:,i) of matrix A by the 2D matrix B
        B_ = permute(A,[2 1 3]);
        A_ = B.';
        szA_ = size(A_);
        szB_ = size(B_);
        if length(szB_)==3, szB_3=szB_(3); else, szB_3=1; end
        AB_ = reshape(A_*reshape(B_,szB_(1),szB_(2)*szB_3,1),[szA_(1) szB_(2) szB_3]);
        AB = permute(AB_,[2 1 3]);
    end

    function rotationVector = get_rotationVector(obj,E_1,E_2)
        E1x = E_1(:,1,:);
        E1y = E_1(:,2,:);
        E1z = E_1(:,3,:);
        E2x = E_2(:,1,:);
        E2y = E_2(:,2,:);
        E2z = E_2(:,3,:);

        Vector = (...
            obj.MultiProd(getSkewMat(E1x),E2x)+...
            obj.MultiProd(getSkewMat(E1y),E2y)+...
            obj.MultiProd(getSkewMat(E1z),E2z));
        Magnitude = sum(Vector.^2,1).^0.5;
        Magnitude = max(Magnitude,1e-15);
        NormalVector = Vector./Magnitude;
        rotationVector = NormalVector.*asin(Magnitude/2);

%         rotationVector = asin(1/2.*(...
%             O.MultiProd(getSkewMat(E1x),E2x)+...
%             O.MultiProd(getSkewMat(E1y),E2y)+...
%             O.MultiProd(getSkewMat(E1z),E2z)));
    end

    function int_val = intVal(x,y,int_fnc,dim)
        if nargin < 4
            dim = 3;
        end
        fullValOnly = true;
        int_val = int_fnc(x,y,fullValOnly,dim);
    end

    end

    methods %special copy method
        function obj_copy = copy(obj) 
            %bit of a workaround!
            %produce a deep copy of the object 'obj' preserving all handle structures within the object but severing all links between 'obj_copy' and 'obj'
            %accomplished via save and reload operation on obj

            unique_identifyer = num2str(floor(rand()*1e12));
            temp_fileName = ['temp_NBS_object_copy_' unique_identifyer '.mat'];
            save(temp_fileName,'obj');
            struct = load(temp_fileName);
            varName = fields(struct);
            obj_copy = struct.(varName{1});
            delete(temp_fileName);
        end
    end

    methods (Static = true, Hidden = true)

        function [] = draw_genericPart(s_draw,E,r,w,h,beam_cntr,varargin)

        %A generic beam drawing function that may be accessed and utilized
        %by the .draw methods of each aircraft part

        ns = numel(s_draw);
        r = squeeze(r);
        if numel(beam_cntr)==1, beam_cntr(1:ns) = beam_cntr; end
        ns_draw = length(s_draw); L_ = s_draw(end);
        CrossSectionProfiles_ = utility_functions.get_option(varargin,'CrossSectionProfiles','box');
        CamLight = utility_functions.get_option(varargin,'camlight',[45,-90]);
        idx_ribs_ = utility_functions.get_option(varargin,'idx_ribs',[]);
        idx_intrinsic = utility_functions.get_option(varargin,'idx_intrinsic',[]);
        XLIM = utility_functions.get_option(varargin,'XLIM',[]);
        YLIM = utility_functions.get_option(varargin,'YLIM',[]);
        ZLIM = utility_functions.get_option(varargin,'ZLIM',[]);
        plotReferenceLine = utility_functions.get_option(varargin,'plotReferenceLine',true);
        plot2DProjections = utility_functions.get_option(varargin,'plot2DProjections',[1,1,1]);

        if length(w) == 1, w = ones(1,ns_draw)*w; end

        if strcmp(CrossSectionProfiles_,'box')
            CrossSectionProfiles_ = {L_,linspace(0,1,15),ones(1,15)/2,[linspace(1,0,15),0],[-ones(1,15),1]/2};
        elseif strcmp(CrossSectionProfiles_,'NACA0012')
            x = linspace(0,1,51);
            y = (0.2969*(x).^0.5-0.1260*(x)-0.3516*(x).^2+0.2843*(x).^3-0.1015*(x).^4);
            y = y/max(y)/2;
            CrossSectionProfiles_ = {L_,x(1:end-1),y(1:end-1),fliplr(x),fliplr(-y)};
        elseif strcmp(CrossSectionProfiles_,'circ')
            a_ = linspace(0,2*pi,51);
            x = (cos(a_)+1)/2;
            y = sin(a_)/2;
            CrossSectionProfiles_ = {L_,x(1:25),y(1:25),x(26:51),y(26:51)};
        end
        if isempty(idx_ribs_)
            no_ribs = 20;
            idx_ribs_ = ns_draw:-floor((ns_draw-1)/(no_ribs-1)):1; idx_ribs_(end) = 1;
        end

        %//////////////////////////////////////////////////////////////////////////

%         campos([1 mean(ylim) mean(zlim)]);
%         camtarget([0 mean(ylim) mean(zlim)]);

        draw_section(s_draw,E,r,w,h,beam_cntr,CrossSectionProfiles_,idx_ribs_,plotReferenceLine,CamLight);

        if ~isempty(idx_intrinsic)
            for i_ = idx_intrinsic
                plot3_pairs([r(:,i_) r(:,i_)+w(end)*L_/8*E(1:3,i_)],[1 2],{'r','LineWidth',2});%ex
                plot3_pairs([r(:,i_) r(:,i_)+w(end)*L_/8*E(4:6,i_)],[1 2],{'g','LineWidth',2});%ey
                plot3_pairs([r(:,i_) r(:,i_)+w(end)*L_/8*E(7:9,i_)],[1 2],{'b','LineWidth',2});%ez
            end
        end

        hold on;
        xlim(XLIM); ylim(YLIM); zlim(ZLIM);
        if plot2DProjections
            plot_2D_projections(r,XLIM,YLIM,ZLIM,plot2DProjections);
        end

        function [] = draw_section(s_draw,E,r,w,h,beam_cntr,CrossSectionProfiles_,idx_ribs,plotReferenceLine,CamLight)
            % draw airfoil section profile over plotting elements
            if plotReferenceLine
                plot3(r(1,:),r(2,:),r(3,:),'color',[0 0.5 0],'lineWidth',2);
            end

            nx = length(CrossSectionProfiles_{1,2}) + length(CrossSectionProfiles_{1,4});
            no_stringers = 2;

            P_airfoil_prev(3,2*nx) = 0;

            AP_counter = 1;
            CrossSectionProfiles_ = [CrossSectionProfiles_(1,:) ; CrossSectionProfiles_]; CrossSectionProfiles_{1,1} = 0;
            s_aeroProfiles = CrossSectionProfiles_(:,1);
            P_Airfoil = zeros(3,nx,numel(s_draw));

            for j_ = 1:size(r,2)

                if s_draw(j_)>s_aeroProfiles{AP_counter+1}
                    AP_counter = AP_counter + 1;
                end
                profileScalingFactor = (s_draw(j_) - s_aeroProfiles{AP_counter})/(s_aeroProfiles{AP_counter+1} - s_aeroProfiles{AP_counter});

                idx_stringers = 1:floor((nx-1)/(no_stringers)):nx;
                idx_stringers(end) = [];
                idx_stringers = [idx_stringers ; idx_stringers+nx].';

                ex = E(1:3,j_);
                %ey = E(4:6,j_);
                ez = E(7:9,j_);
                stripVec = ex;

                xu = CrossSectionProfiles_{AP_counter,2}*(1-profileScalingFactor) + CrossSectionProfiles_{AP_counter+1,2}*(profileScalingFactor);
                yu = CrossSectionProfiles_{AP_counter,3}*(1-profileScalingFactor) + CrossSectionProfiles_{AP_counter+1,3}*(profileScalingFactor);
                P_airfoil_up = bsxfun(@plus,r(:,j_),bsxfun(@times,(xu-beam_cntr(j_))*w(j_),stripVec)+bsxfun(@times,yu*h(j_),ez));

                xd = CrossSectionProfiles_{AP_counter,4}*(1-profileScalingFactor) + CrossSectionProfiles_{AP_counter+1,4}*(profileScalingFactor);
                yd = CrossSectionProfiles_{AP_counter,5}*(1-profileScalingFactor) + CrossSectionProfiles_{AP_counter+1,5}*(profileScalingFactor);
                P_airfoil_dn = bsxfun(@plus,r(:,j_),bsxfun(@times,(xd-beam_cntr(j_))*w(j_),stripVec)+bsxfun(@times,yd*h(j_),ez));

                P_airfoil = [P_airfoil_up P_airfoil_dn];
                Paf = [P_airfoil_prev P_airfoil];
                if j_>1
                    plot3_pairs(Paf,idx_stringers,{'color',[0.4,0.4,0.4]});
                end
                if ismember(j_,idx_ribs)
                    %fill3(P_airfoil(1,:),P_airfoil(2,:),P_airfoil(3,:),'r','edgeColor','r','faceAlpha',0.4);
                    fill3(P_airfoil(1,:),P_airfoil(2,:),P_airfoil(3,:),'k','edgeColor',[0.4,0.4,0.4],'faceAlpha',0.0);
                    %fill3(P_airfoil(1,:),P_airfoil(2,:),P_airfoil(3,:),'r','edgeColor','r','faceAlpha',1.0);
                    %surf(P_airfoil(1,:),P_airfoil(2,:),P_airfoil(3,:));
                    %plot3(r(1,j_),r(2,j_),r(3,j_),'marker','+','markerEdgeColor','b');
                end
                P_airfoil_prev = P_airfoil;
                P_Airfoil(:,:,j_) = P_airfoil;
            end
            Xsurf = squeeze(P_Airfoil(1,:,:));
            Ysurf = squeeze(P_Airfoil(2,:,:));
            Zsurf = squeeze(P_Airfoil(3,:,:));
            surfaceOptions = {...
                'FaceAlpha',0.7+0.3,...
                'FaceColor',[0.6,0.7,0.8]+0.2,...
                'MeshStyle','row',...
                'edgeColor','none',...
                'FaceLighting','gouraud'};
            surf(Xsurf,Ysurf,Zsurf,surfaceOptions{:});
            camlight(CamLight(1),CamLight(2));
        end

        function [] = plot3_pairs(z,idc,args)
            %z: 3xn, vector of all points
            %idc: nx2, list of index pairs referencing z
            for j_ = 1:size(idc,1)
                idc_ = idc(j_,:);
                plot3(z(1,idc_),z(2,idc_),z(3,idc_),args{:});
            end
        end

        function [] = plot_2D_projections(Gamma,XLIM,YLIM,ZLIM,projectionFacesXYZ)
            GammaX = Gamma(1,:); GammaY = Gamma(2,:); GammaZ = Gamma(3,:);

            XLIM_projection = XLIM(projectionFacesXYZ(1));
            YLIM_projection = YLIM(projectionFacesXYZ(2));
            ZLIM_projection = ZLIM(projectionFacesXYZ(3));

            hold on;
            zrs = GammaX*0;
            plot3(zrs+XLIM_projection,GammaY,GammaZ,'color',[0.75,0.75,0.75],'lineWidth',1.5);
            plot3(GammaX,GammaY,zrs,'color',[0,0.25,0.75],'lineWidth',1.5);
            
%             plot3(GammaX,zrs+YLIM_projection,GammaZ,'color',[0,0.75,0.75],'lineWidth',1.5);
            
            plot3(GammaX,GammaY,zrs+ZLIM_projection,'color',[0.75,0.75,0.75],'lineWidth',1.5);
            plot3(zrs+GammaX(end),zrs+GammaY(end),GammaZ,'color',[0,0.25,0.75],'lineWidth',1.5);
            
            plot3([0,XLIM_projection],[0,0],[0,0],'color',[0.75,0.75,0.75],'lineStyle',':','lineWidth',1);
            plot3([0,0],[0,0],[0,ZLIM_projection],'color',[0.75,0.75,0.75],'lineStyle',':','lineWidth',1);
            plot3([GammaX(end),XLIM_projection],[0,0]+GammaY(end),[0,0]+GammaZ(end),'color',[0.75,0.75,0.75],'lineStyle',':','lineWidth',1);
            plot3([0,0]+GammaX(end),[0,0]+GammaY(end),[GammaZ(end),ZLIM_projection],'color',[0.75,0.75,0.75],'lineStyle',':','lineWidth',1);

        end

        end

    end

    methods
        function set.uVec_freeStream_G(obj,val_or_func)
            if isa(val_or_func,'function_handle')
                val = val_or_func(0);
            else
                val = val_or_func;
            end
            assert(numel(val)==3,'uVec_freeStream must return a vector of length 3');
            if ~utility_functions.eq_tol(norm(val),1,'relTol',1e-6)
                warning('uVec_freeStream did not return a unit vector')
            end

            obj.uVec_freeStream_G = val_or_func;
        end

        function set.flexParts_nonlinear_cell(obj,val)
            obj.flexParts_nonlinear_cell = val;
            obj.nflexParts_nonlinear = numel(val);
        end

        function set.t(obj,val)
            obj.t = val(:).';
            obj.nt = numel(val);
        end
    end

end
