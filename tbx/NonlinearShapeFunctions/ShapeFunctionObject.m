classdef ShapeFunctionObject

    properties
        shapeSetTemplate                                                   %[string] template set type upon which object is based (e.g. Chebyshev 1st)
        BCs                                                                %[-] requested boundary conditions of the shape set. BCs = [displacement_left disp_right ; 1st_derivative_left 1stDrv_right ; 2ndDrv_left 2ndDrv_right] where entries are either '1' non-zero, or, '0' zero. E.g. for an unforced cantilever beam BCs = [0 1;1 0;1 0]
        halfShape
        s
        x, dx_ds
        E, dE_ds
        customShapes, dcustomShapes_ds
        setHistory
        isOrthogonal
        setName                                                            %optional title string to describe entire object
    end

    properties (Dependent)
        y, y_all
        dy_dx
        dy_ds, dy_ds_all
        weightFunction, weightFunction_all
        nShapes
    end

    methods

        function obj = ShapeFunctionObject(shapeSetTemplate,varargin)

            %addpath('./utility_functions');

            obj.s = utility_functions.get_option(varargin,'s',linspace(0,1,101)); obj.s = reshape(obj.s,1,[],1);
            obj.BCs = utility_functions.get_option(varargin,'BCs',[0 1;1 1;1 1]);
            nShapes = utility_functions.get_option(varargin,'nShapes',10);
            obj.halfShape = utility_functions.get_option(varargin,'halfShape',false);
            obj.setName = utility_functions.get_option(varargin,'setName',[]);
            obj.setHistory{1,3} = zeros(nShapes,1);

            ShapeSetTemplates = {...
                'polynomial',...
                'reverse_poly',...
                'cantilever_mode_shapes',...
                'trigonometric',...
                'chebyshev_1st',...
                'chebyshev_2nd',...
                'legendre'};

            obj.shapeSetTemplate = lower(shapeSetTemplate);
            assert(ismember(obj.shapeSetTemplate,ShapeSetTemplates),'Argument must specify an existing shape set template from ''ShapeSetTemplates''');

            obj = generateScalingFunctions(obj);
            obj = get_templateSet(obj);
            obj = applyScaling(obj);

        end

    end

    methods

        function obj = get_templateSet(obj)
            switch obj.shapeSetTemplate
                case 'chebyshev_1st'
                    [y_,dy_dx_,weightFunction_] = obj.chebyshev(obj.x,obj.nShapes,'1st');
                case 'chebyshev_2nd'
                    [y_,dy_dx_,weightFunction_] = obj.chebyshev(obj.x,obj.nShapes,'2nd');
                case 'legendre'
                    [y_,dy_dx_,weightFunction_] = obj.legendre(obj.x,obj.nShapes);
                case 'polynomial'
                    [y_,dy_dx_,weightFunction_] = obj.polynomial(obj.x,obj.nShapes);
            end
            obj.setHistory = {[obj.shapeSetTemplate ' template set']...
                , 'y(x)' , y_ , 'dy_dx(x)' , dy_dx_ , 'W(x)' , weightFunction_};

            obj.isOrthogonal = ~isempty(obj.weightFunction);
        end

        function obj = applyScaling(obj)
            y_ = obj.y.*obj.E;
            dy_ds_ = (obj.dy_dx.*obj.dx_ds.*obj.E) + obj.y.*obj.dE_ds;
            weightFunction_ = (obj.weightFunction.*obj.dx_ds)./max(obj.E,1e-3).^2;

            obj.setHistory = [obj.setHistory ; {'shape set scaled using envelope E(s) and mapping function x(s)'...
                , 'y(s) = E y(x(s))' , y_ , 'dy_ds(s)' , dy_ds_ , 'W(s) = ( W(x(s)) / E(s)^2 ) * dx_ds' , weightFunction_}];
        end

        function obj = generateScalingFunctions(obj)
            %BC_str = num2str(reshape(BCs.',1,[]));
            d0 = obj.BCs(1,:);
            d1 = obj.BCs(2,:);
            d2 = obj.BCs(3,:);

            L = (max(obj.s) - min(obj.s));
            s01 = (obj.s - min(obj.s)) / L;

            switch obj.shapeSetTemplate
                case {'chebyshev_1st','chebyshev_2nd','legendre'}
                    %allowable_BCs = {'','','',''}
                    %ismember check
                    d12_d22 = [d1(2) d2(2)];
                    if isequal(d0,[0 1]) && d1(1)==1 && d2(1)==1
                        if isequal(d12_d22,[1 1])
                            E_ = ((s01).^3-3*(s01).^2+3*(s01)); dE_ds_ = (3*(s01).^2-6.*(s01)+3)*(1/L);
                            xcheby = -2*(s01)+1; dxcheby_ds = -2/L;
                        elseif isequal(d12_d22,[0 1])
                            E_ = sin(s01*pi/2); dE_ds_ = pi/(2*L)*cos(s01*pi/2);
                            xcheby =2*(sin((-s01)*pi/2))+1; dxcheby_ds = 2*-pi/(2*L)*cos((-s01)*pi/2);
                        elseif isequal(d12_d22,[0 0])
                            E_ = ((s01).^3-3*(s01).^2+3*(s01)); dE_ds_ = (3*(s01).^2-6.*(s01)+3)*(1/L);
                            S = -2*(s01)+1;
                            xcheby = (1/5.*S.^5-2/3.*S.^3+S)*15/8; dxcheby_ds = (S.^4-2.*S.^2+1)*(-2/L*15/8);
                        else
                            error('chebyshev boundary condition not available')
                        end
                    elseif isequal(d0,[1 1]) && d1(1)==1 && d2(1)==1
                        E_ = 1; dE_ds_ = 0;
                        if isequal(d12_d22,[1 1])
                            xcheby = -2*(s01)+1; dxcheby_ds = -2/L;
                        elseif isequal(d12_d22,[0 1])
                            xcheby = 2*(sin((-s01)*pi/2))+1; dxcheby_ds = 2*-pi/(2*L)*cos((-s01)*pi/2);
                        elseif isequal(d12_d22,[0 0])
                            S = -2*(s01)+1;
                            xcheby = (S.^3+3*S.^2+3*S-3)*1/4; dxcheby_ds = (3*S.^2+6.*S+3)*(-2/L*1/4);
                        else
                            error('chebyshev boundary condition not available')
                        end
                    else
                        error('chebyshev boundary condition not available')
                    end
                    obj.E = E_;
                    obj.dE_ds = dE_ds_;
                    obj.x = (1-xcheby)/2;
                    obj.dx_ds = dxcheby_ds*(-1/2);

                case 'polynomial'
                    obj.x = s01;
                    obj.dx_ds = 1 / L;
                    obj.E = 1;
                    obj.dE_ds = 0;
            end

            if obj.halfShape ~= false
                if isequal(obj.halfShape,'left')
                    obj.x = obj.x/2;
                    obj.dx_ds = obj.dx_ds/2;
                elseif isequal(obj.halfShape,'right')
                    obj.x = (obj.x+1)/2;
                    obj.dx_ds = obj.dx_ds/2;
                end
            end

        end

        function obj = addCustomFunctions(obj,s_custom,y_custom,dy_ds_custom,varargin)

            s_custom = reshape(s_custom,1,[]); ns = numel(s_custom);
            y_custom = reshape(y_custom,[],ns,1);

            if ~isempty(dy_ds_custom)
                dy_ds_custom = reshape(dy_ds_custom,[],ns,1);
            else
                del_s = diff(s_custom);
                dy_ds_custom = zeros(size(y_custom));
                for i_ = 1:size(y_custom,1)
                    del_y = diff(y_custom(i_,:));
                    dy_ds_custom(i_,:) = ([0 del_y] + [del_y 0]) ./ ([0 del_s] + [del_s 0]);
                end
            end
            sMat = utility_functions.sampleMat(s_custom,obj.s);
            assert(max(abs(s_custom*sMat-obj.s)) < (max(obj.s)-min(obj.s))*1e-8 , 'error with custom function sampling matrix ''sMat''');

            yCustom = y_custom*sMat;
            dydsCustom = dy_ds_custom*sMat;
            obj.customShapes = yCustom;
            obj.dcustomShapes_ds = dydsCustom;

            delete_ = get_option(varargin,'delete',[]);
            y_ = obj.y; y_(delete_,:) = []; y_ = [yCustom;y_];
            dy_ds_ = obj.dy_ds; dy_ds_(delete_,:) = []; dy_ds_ = [dydsCustom;dy_ds_];

            obj.setHistory = [obj.setHistory ; {'append custom functions to base set'...
                , 'y(s) subset of union( E y(x(s)) , y_custom(s) )' , y_ , 'dy_ds(s)' , dy_ds_ , 'W(s)' , obj.weightFunction}];
            obj.isOrthogonal = [];
        end

        function obj = orthNorm(obj)
            %==========================================================================
            %normalise the orthogonal shape set obj.y(s)
            %==========================================================================
            y_ = obj.y; dy_ds_ = obj.dy_ds;
            s_ = obj.s; weightFunction_ = obj.weightFunction;

            for i_ = 1:obj.nShapes
                orthMagnitude = obj.orth_integral(s_,y_(i_,:),y_(i_,:),weightFunction_)^0.5;

                y_(i_,:) = y_(i_,:)/orthMagnitude;
                dy_ds_(i_,:) = dy_ds_(i_,:)/orthMagnitude;
            end
            obj.setHistory = [obj.setHistory ; {'orthonormalise the shape set using supplied weightFunction'...
                , 'y(s) norm (y(s))' , y_ , 'dy_ds(s)' , dy_ds_ , 'W(s)' , obj.weightFunction}];
            obj.isOrthogonal = true;
        end

        function obj = orth(obj)
            %==========================================================================
            %performs Gram Schmidt orthogonalisation on the shape set
            %obj.y(s) given the weighting function obj.weightFunction(s)
            %==========================================================================
            y_ = obj.y; dy_ds_ = obj.dy_ds;
            s_ = obj.s; weightFunction_ = obj.weightFunction;

            for i_ = 2:obj.nShapes
                y_corrections = zeros(i_-1,numel(obj.s));
                dy_corrections = zeros(i_-1,numel(obj.s));
                for j_ = 1:(i_-1)
                    correctionFactor = obj.orth_integral(s_,y_(i_,:),y_(j_,:),weightFunction_);
                    y_corrections(j_,:) = correctionFactor * y_(j_,:);
                    dy_corrections(j_,:) = correctionFactor * dy_ds_(j_,:);
                end
                ith_y_correction = sum(y_corrections,1);
                ith_dy_correction = sum(dy_corrections,1);

                y_(i_,:) = y_(i_,:) - ith_y_correction;
                dy_ds_(i_,:) = dy_ds_(i_,:) - ith_dy_correction;

                orthMagnitude = obj.orth_integral(s_,y_(i_,:),y_(i_,:),weightFunction_)^0.5;
                y_(i_,:) = y_(i_,:)/orthMagnitude;
                dy_ds_(i_,:) = dy_ds_(i_,:)/orthMagnitude;
            end
            obj.setHistory = [obj.setHistory ; {'reorthogonalise shape set using supplied weightFunction'...
                , 'y(s) orthogonalise(y(s))' , y_ , 'dy_ds(s)' , dy_ds_ , 'W(s)' , obj.weightFunction}];
            obj.isOrthogonal = true;
        end

        function [] = plotShapes(obj,varargin)

            output_detail = get_option(varargin,'output_detail','final');

            switch output_detail

                case 'final_1' %plot only the final shape set

                    figure('windowStyle','docked');
                    plot(obj.s,obj.y_all{end}); xlabel('s'); ylabel(obj.setHistory{end,2}); title(obj.setName,'interpreter','none');

                case 'final_2' %plot only the final shape set plus a few supporting panels

                    figure('windowStyle','docked');
                    subplot(6,3,[1 2 4 5 7 8]), plot(obj.s,obj.y_all{end}); xlabel('s'); ylabel(obj.setHistory{end,2}); title(obj.setName,'interpreter','none');
                    subplot(6,3,[10 11 13 14 16 17]), plot(obj.s,obj.dy_ds_all{end}); xlabel('s'); ylabel(obj.setHistory{end,4});
                    subplot(6,3,[3 6]), plot(obj.s,obj.x); xlabel('s'); ylabel('x');
                    subplot(6,3,[9 12]), plot(obj.s,obj.E); xlabel('s'); ylabel('E');
                    subplot(6,3,[15 18]), plot(obj.s,obj.weightFunction_all{end}); xlabel('s'); ylabel(obj.setHistory{end,6});

                case 'full_panel' %plot final set plus all intermediate steps in the same figure window

                    nRow = size(obj.setHistory,1)+1;

                    figure;
                    subplot(nRow,3,1), plot(obj.x,obj.y_all{1}); xlabel('x'); ylabel(obj.setHistory{1,2}); title(obj.setHistory{1,1},'interpreter','none');
                    subplot(nRow,3,2), plot(obj.x,obj.dy_dx); xlabel('x'); ylabel(obj.setHistory{1,4});
                    hold on, x = obj.x; y = obj.y_all{1}; xDiff = diff(x); yDiff = diff(y.').'; plot((x(1:end-1)+x(2:end))/2,bsxfun(@times,yDiff,1./xDiff),'r:');
                    subplot(nRow,3,3), plot(obj.x,obj.weightFunction_all{1}); xlabel('x'); ylabel(obj.setHistory{1,6});

                    subplot(nRow,3,1+3), plot(obj.s,obj.y_all{1}); xlabel('s'); ylabel('y(x(s))');
                    subplot(nRow,3,2+3), plot(obj.s,bsxfun(@times,obj.dy_dx,obj.dx_ds)); xlabel('s'); ylabel('dy_ds');
                    hold on, s = obj.s; y = obj.y_all{1}; sDiff = diff(s); yDiff = diff(y.').'; plot((s(1:end-1)+s(2:end))/2,bsxfun(@times,yDiff,1./sDiff),'r:');
                    subplot(nRow,3,3+3), plot(obj.s,obj.x); xlabel('s'); ylabel('x');

                    subplot(nRow,3,1+6), plot(obj.s,obj.y_all{2}); xlabel('s'); ylabel(obj.setHistory{2,2}); title(obj.setHistory{2,1},'interpreter','none');
                    subplot(nRow,3,2+6), plot(obj.s,obj.dy_ds_all{2}); xlabel('s'); ylabel(obj.setHistory{2,4});
                    hold on, s = obj.s; y = obj.y_all{2}; sDiff = diff(s); yDiff = diff(y.').'; plot((s(1:end-1)+s(2:end))/2,bsxfun(@times,yDiff,1./sDiff),'r:');
                    subplot(nRow,3,3+6), plot(obj.s,obj.E); xlabel('s'); ylabel('E');

                    %             if ~isempty(obj.customShapes)
                    %             figure('windowStyle','docked');
                    %             subplot(1,3,1), plot(obj.s,obj.customShapes,'lineWidth',2,'color','r'); hold on;
                    %                             plot(obj.s,obj.y_all{3}(size(obj.customShapes,1)+1:end,:)); xlabel('s'); ylabel(obj.setHistory{3,2}); title(obj.setHistory{3,1});
                    %             subplot(1,3,2), plot(obj.s,obj.dcustomShapes_ds,'lineWidth',2,'color','r'); hold on;
                    %                             plot(obj.s,obj.dy_ds_all{3}(size(obj.customShapes,1)+1:end,:)); xlabel('s'); ylabel(obj.setHistory{3,4});
                    %             subplot(1,3,3),
                    %             end

                    for i_ = 3:size(obj.setHistory,1)
                        subplot(nRow,3,1+i_*3), plot(obj.s,obj.y_all{i_}); xlabel('s'); ylabel(obj.setHistory{i_,2}); title(obj.setHistory{i_,1},'interpreter','none');
                        subplot(nRow,3,2+i_*3), plot(obj.s,obj.dy_ds_all{i_}); xlabel('s'); ylabel(obj.setHistory{i_,4});
                        subplot(nRow,3,3+i_*3), plot(obj.s,obj.weightFunction_all{i_}); xlabel('s'); ylabel(obj.setHistory{i_,6});
                    end

                case 'full_docked' %plot final set plus all intermediate steps in separate docked tabs

                    figure('windowStyle','docked');
                    subplot(1,3,1), plot(obj.x,obj.y_all{1}); xlabel('x'); ylabel(obj.setHistory{1,2}); title(obj.setHistory{1,1},'interpreter','none');
                    subplot(1,3,2), plot(obj.x,obj.dy_dx); xlabel('x'); ylabel(obj.setHistory{1,4});
                    checkGradient(obj.x,obj.y_all{1});
                    subplot(1,3,3), plot(obj.x,obj.weightFunction_all{1}); xlabel('x'); ylabel(obj.setHistory{1,6});

                    figure('windowStyle','docked');
                    subplot(1,3,1), plot(obj.s,obj.y_all{1}); xlabel('s'); ylabel('y(x(s))');
                    subplot(1,3,2), plot(obj.s,bsxfun(@times,obj.dy_dx,obj.dx_ds)); xlabel('s'); ylabel('dy_ds');
                    checkGradient(obj.s,obj.y_all{1});
                    subplot(1,3,3), plot(obj.s,obj.x); xlabel('s'); ylabel('x');

                    figure('windowStyle','docked');
                    subplot(1,3,1), plot(obj.s,obj.y_all{2}); xlabel('s'); ylabel(obj.setHistory{2,2}); title(obj.setHistory{2,1},'interpreter','none');
                    subplot(1,3,2), plot(obj.s,obj.dy_ds_all{2}); xlabel('s'); ylabel(obj.setHistory{2,4});
                    checkGradient(obj.s,obj.y_all{2});
                    subplot(1,3,3), plot(obj.s,obj.E); xlabel('s'); ylabel('E');

                    %             if ~isempty(obj.customShapes)
                    %             figure('windowStyle','docked');
                    %             subplot(1,3,1), plot(obj.s,obj.customShapes,'lineWidth',2,'color','r'); hold on;
                    %                             plot(obj.s,obj.y_all{3}(size(obj.customShapes,1)+1:end,:)); xlabel('s'); ylabel(obj.setHistory{3,2}); title(obj.setHistory{3,1});
                    %             subplot(1,3,2), plot(obj.s,obj.dcustomShapes_ds,'lineWidth',2,'color','r'); hold on;
                    %                             plot(obj.s,obj.dy_ds_all{3}(size(obj.customShapes,1)+1:end,:)); xlabel('s'); ylabel(obj.setHistory{3,4});
                    %             subplot(1,3,3),
                    %             end

                    for i_ = 3:size(obj.setHistory,1)
                        figure('windowStyle','docked');
                        subplot(1,3,1), plot(obj.s,obj.y_all{i_}); xlabel('s'); ylabel(obj.setHistory{i_,2}); title(obj.setHistory{i_,1},'interpreter','none');
                        subplot(1,3,2), plot(obj.s,obj.dy_ds_all{i_}); xlabel('s'); ylabel(obj.setHistory{i_,4});
                        checkGradient(obj.s,obj.y_all{i_});
                        subplot(1,3,3), plot(obj.s,obj.weightFunction_all{i_}); xlabel('s'); ylabel(obj.setHistory{i_,6});
                    end

            end

            function [] = checkGradient(x,y)
                hold on
                xDiff = diff(x);
                yDiff = diff(y.').';
                plot((x(1:end-1)+x(2:end))/2,bsxfun(@times,yDiff,1./xDiff),'r:');
            end

        end

        function isOrth = orthTest(obj,varargin)

            PLOT = get_option(varargin,'PLOT',false);

            nShapes = obj.nShapes;
            orthMatrix = zeros(nShapes);
            %--------------------------------------------------------------
            % calculate orthogonal integrals
            for i_ = 1:nShapes
                for j_ = 1:nShapes
                    orthMatrix(i_,j_) = obj.orth_integral(obj.s,obj.y(i_,:),obj.y(j_,:),obj.weightFunction);
                end
            end
            maxOrthVal = max(abs(orthMatrix(:)));
            %--------------------------------------------------------------
            if PLOT
                % plot orthogonal integrals
                figure;
                for i_ = 1:nShapes
                    for j_ = 1:nShapes
                        subplot(nShapes,nShapes,i_+(j_-1)*nShapes), area(obj.s,obj.y(i_,:).*obj.y(j_,:).*obj.weightFunction), hold on, plot(obj.s,0);
                        set(gca,'color',[1-abs(orthMatrix(i_,j_)/maxOrthVal) 1 1]);
                    end
                end
                %-----------------------------------
                % plot orthMatrix
                figure; hold on;
                rectangle('position',[0.5 -nShapes-0.5 nShapes nShapes]);
                for i_ = 1:nShapes,
                    for j_ = 1:nShapes,
                        plot(i_,-j_,'ko','markerFaceColor','blue','markerSize',abs(orthMatrix(i_,j_)/maxOrthVal)*20+0.01);
                    end
                end; hold off
            end
            %--------------------------------------------------------------
            isOrth = isdiag((abs(orthMatrix) > maxOrthVal/50)+0);
        end

    end

    methods %get, set, methods for dependent properties
        function y_all = get.y_all(obj)
            y_all = obj.setHistory(:,3);
        end
        function y = get.y(obj)
            y = obj.y_all{end};
        end
        function dy_dx = get.dy_dx(obj)
            dy_dx = obj.setHistory{1,5};
        end
        function dy_ds_all = get.dy_ds_all(obj)
            dy_ds_all = obj.setHistory(:,5);
            dy_ds_all{1} = bsxfun(@times,dy_ds_all{1},obj.dx_ds);
        end
        function dy_ds = get.dy_ds(obj)
            dy_ds = obj.dy_ds_all{end};
        end
        function weightFunction_all = get.weightFunction_all(obj)
            weightFunction_all = obj.setHistory(:,7);
        end
        function weightFunction = get.weightFunction(obj)
            weightFunction = obj.weightFunction_all{end};
        end
        function nShapes = get.nShapes(obj)
            nShapes = size(obj.y,1);
        end
        function obj = set.halfShape(obj,str)
            assert(isequal(str,false) || ismember(str,{'left','right'}),'optional keyword ''halfShape'' must take a value in the set {false,''left'',''right''}');
            obj.halfShape = str;
        end
    end

    methods (Static) %methods that return default shape sets

        function [y,dy_dx,weightFnc] = chebyshev(x,N,type)

            xx = -2*x+1;
            dxx_dx = -2;

            cheby_1st = zeros(N,length(xx));
            cheby_2nd = zeros(N,length(xx));
            dcheby_dxx_1st = zeros(N,length(xx));

            %1st kind recursive definitions -------------------------------
            if N>0
                cheby_1st(1,:) = 1;
                dcheby_dxx_1st(1,:) = 0;
            end

            if N>1
                cheby_1st(2,:) = xx;
                dcheby_dxx_1st(2,:) = 1;
            end

            for n = 3:N
                % y[n] = 2.x.y[n-1] - y[n-2]
                cheby_1st(n,:) = 2*xx.*cheby_1st(n-1,:) - cheby_1st(n-2,:);
                % dy[n] = (n-1).(2.y[n-1] + 1/(n-3).dy[n-2])
                dcheby_dxx_1st(n,:) = (n-1)*(2*cheby_1st(n-1,:)+1/(max(n-3,1)).*dcheby_dxx_1st(n-2,:));
            end

            %2nd kind definitions, calculated from 1st kind ---------------

            cheby_2nd(1:2:N,:) = 2*cumsum(cheby_1st(1:2:N,:),1)-1;
            cheby_2nd(2:2:N,:) = 2*cumsum(cheby_1st(2:2:N,:),1);
            dcheby_dxx_2nd(1:2:N,:) = 2*cumsum(dcheby_dxx_1st(1:2:N,:),1);
            dcheby_dxx_2nd(2:2:N,:) = 2*cumsum(dcheby_dxx_1st(2:2:N,:),1);

            switch type
                case '1st'
                    y = cheby_1st;
                    dy_dx = dcheby_dxx_1st*dxx_dx;
                    x_ = xx; x_(1) = x_(1)-1e-3; x_(end) = x_(end)+1e-3;    %perturb x limits away from singularities
                    weightFnc = (1-x_.^2).^(-0.5);
                case '2nd'
                    y = cheby_2nd;
                    dy_dx = dcheby_dxx_2nd*dxx_dx;
                    weightFnc = (1-xx.^2).^(0.5);
            end

        end

        function [y,dy_dx,weightFnc] = legendre(x,N)

            xx = -1 + x*2;
            dxx_dx = 2;

            legendre = zeros(N,length(xx));
            dlegendre_dxx = zeros(N,length(xx));

            if N>0
                legendre(1,:) = 1;
                dlegendre_dxx(1,:) = 0;
            end

            if N>1
                legendre(2,:) = xx;
                dlegendre_dxx(2,:) = 1;
            end

            for n = 3:N
                J = n-2;
              % (J+1) P_{J+1}(x) = (2J+1) x P_n(x) - J P_{J-1}(x): J = 1 for 3rd shape
                legendre(n,:) = ( (2*J+1)*x.*legendre(n-1,:) - J*legendre(n-2,:) )/(J+1);
				dlegendre_dxx(n,:) = ( (2*J+1)*legendre(n-1,:) + (2*J+1)*x.*dlegendre_dxx(n-1,:) - J*dlegendre_dxx(n-2,:) )/(J+1);
            end

            y = legendre;
            dy_dx = dlegendre_dxx*dxx_dx;
            weightFnc = x*0+1;

        end

        function [y,dy_dx,weightFnc] = polynomial(x,N)

            y = zeros(N,length(x));
            dy_dx = zeros(N,length(x));

            for n = 1:N
                y(n,:) = x.^n;
                dy_dx(n,:) = n.*x.^(n-1);
                weightFnc = [];
            end
            disp('BCs not available for polynomial shape set')
        end

        function orthInt = orth_integral(x,y1,y2,weightFnc)
            int = integrate2(x,y1.*y2.*weightFnc);
            orthInt = int(end);
        end

    end
end
