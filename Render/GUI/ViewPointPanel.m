classdef ViewPointPanel < ListLayoutComponent
    
    
    properties
        PosColor = [0 .3 0]
        ApColor = [.1 .7 .1]
        ApAlpha = .5
    end
    
    properties ( SetAccess = private )
        u0=0
        v0=0
        radius=0
        Ap=[]
        panel
        plotAxes
        ApAxes
        ApImg
        ApPos
        
        uPosEditor
        uText
        vPosEditor
        vText
        
        U = []
        V = []
        UVLim = [-5 5 -5 5]
        uPlotCenter = 0
        vPlotCenter = 0
        hasUVcoords = false
        
        callbackFcn = []
        
    end
    
    properties ( Access = private )
        uDragCurrent
        vDragCurrent
        middleButtonClicked = false
        buttonClicked = false
        
        ValboxDecimalDigits=2
    end

    
    methods
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%     Constructor     %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = ViewPointPanel(Parent, u0, v0, radius, Ap, U, V, UVLim)
            obj@ListLayoutComponent(Parent);
            if(isa(Parent,'ListLayoutPanel'))
                ParentUI = Parent.panel;
            else
                ParentUI = Parent;
            end
            
            
            
            if(exist('u0','var') && exist('v0','var'))
                obj.u0 = u0;
                obj.v0 = v0;
            end
            if(exist('radius','var')), obj.radius = radius; end
            if(exist('Ap','var')), obj.Ap = Ap; end
            
            %Set optional U,V coordinates of intitial viewpoints
            if(exist('U','var') && exist('V','var'))
                obj.U = U; obj.V = V;
            end
            obj.hasUVcoords = ~isempty(obj.U) && ~isempty(obj.V);
            if(numel(obj.U) ~= numel(obj.V))
                warning('U and V lists of initial view positions have different sizes. They will be ignored.');
                obj.hasUVcoords = false;
                obj.U = []; obj.V = [];
            end
            
            %Set limits U,V coordinates of the panel
            if(exist('UVLim','var'))
                obj.UVLim = UVLim;
                obj.uPlotCenter = mean(obj.UVLim(1:2));
                obj.vPlotCenter = mean(obj.UVLim(3:4));
            elseif(obj.hasUVcoords)
                scale = 3;
                obj.uPlotCenter = (max(obj.U(:))+min(obj.U(:)))/2;%mean(U(:));%
                obj.vPlotCenter = (max(obj.V(:))+min(obj.V(:)))/2;%mean(V(:));%
                obj.UVLim(1:2) = [obj.uPlotCenter - scale*(obj.uPlotCenter-min(obj.U(:))), obj.uPlotCenter + scale*(max(obj.U(:))-obj.uPlotCenter)];
                obj.UVLim(3:4) = [obj.vPlotCenter - scale*(obj.vPlotCenter-min(obj.V(:))), obj.vPlotCenter + scale*(max(obj.V(:))-obj.vPlotCenter)];
                
                %dU = ( max(obj.U(:)) - min(obj.U(:)) )/2;
                %obj.UVLim(1:2) = [min(obj.U(:))-scale*dU, max(obj.U(:))+scale*dU];
                %dV = ( max(obj.V(:)) - min(obj.V(:)) )/2;
                %obj.UVLim(3:4) = [min(obj.V(:))-scale*dV, max(obj.V(:))+scale*dV];
                
            end

            %Main viewpoint panel
            obj.panel = uipanel('Parent', ParentUI, 'Units', 'Pixels');
            obj.panel.set('ResizeFcn',@(src,evnt)ViewPointPanel.panelResize(obj,src,evnt));

            %Axes of the viewpoints
            obj.plotAxes = axes('pos',[0 0 1 1], 'Parent', obj.panel);
            scatter(obj.U(:), obj.V(:), [], [.5,.5,.5], '+', 'Parent', obj.plotAxes);
            line([0,0],[intmin intmax],'Color','red','LineStyle',':','Parent', obj.plotAxes);
            line([intmin intmax],[0 0],'Color','red','LineStyle',':','Parent', obj.plotAxes);
            obj.plotAxes.set('XLim', obj.UVLim(1:2), 'YLim', obj.UVLim(3:4));
            obj.plotAxes.set('XTickMode','manual','YTickMode','manual', 'XColor',[.5 .5 .5], 'YColor',[.5 .5 .5], 'XGrid','on', 'YGrid','on', 'TickLength', [0 0]);
            obj.plotAxes.set('XTick', round(obj.UVLim(1)):round(obj.UVLim(2)) );
            obj.plotAxes.set('YTick', round(obj.UVLim(3)):round(obj.UVLim(4)) );
            
            %Axes of the aperture image
            obj.ApAxes = axes('Parent', obj.panel);hold on;

            %Initialize Aperture image
            obj.ApImg = imshow(bsxfun(@times,ones([size(obj.Ap),3]),permute(obj.ApColor,[3 1 2])), 'Parent', obj.ApAxes);
            obj.ApImg.set('AlphaData',max(0,(obj.Ap*obj.ApAlpha)));

            %Marker of the view position
            obj.ApPos = line(0,0,'Marker','+','MarkerSize',10,'Color',obj.PosColor,'LineWidth',1.5, 'Parent', obj.ApAxes);hold off

            %Invisible axes on top of other elements only used to receive mouse clicks.
            clickAxes = axes('Parent', obj.panel, 'Color', 'none', 'XColor','none', 'YColor','none', 'Position',[0 0 1 1], 'Visible','on');
            clickAxes.set('ButtonDownFcn',@(src,evnt)ViewPointPanel.mouseClick(obj,src,evnt));
            
            %Text box editors for u and v positions
            obj.uText = uicontrol('Parent', obj.panel, 'Style','text','String','u=', 'Units','Pixels', 'Position',[5,22,20,14], 'BackgroundColor','white');
            obj.uPosEditor = uicontrol('Parent',obj.panel, 'Style','edit', 'Units','Pixels', 'Position',[21,22,33,14]);
            set(obj.uPosEditor, 'Callback', @(src,evnt) ViewPointPanel.editUVPos(obj,src,evnt));
            obj.vText = uicontrol('Parent', obj.panel, 'Style','text','String','v=', 'Units','Pixels', 'Position',[5,5,20,14], 'BackgroundColor','white');
            obj.vPosEditor = uicontrol('Parent',obj.panel, 'Style','edit', 'Units','Pixels', 'Position',[21,5,33,14]);
            set(obj.vPosEditor, 'Callback', @(src,evnt) ViewPointPanel.editUVPos(obj,src,evnt));
            
            %Callbacks for mouse events
            addlistener(gcf,'WindowMouseMotion',@(src,evnt)ViewPointPanel.mouseMotion(obj,src,evnt));
            addlistener(gcf,'WindowMouseRelease',@(src,evnt)ViewPointPanel.mouseRelease(obj,src,evnt));
            addlistener(gcf,'WindowScrollWheel',@(src,evnt)ViewPointPanel.mouseScrollWheel(obj,src,evnt));
            
            obj.setViewPosition(obj.u0,obj.v0);
            
        end
        
%%%%%%%%%%%%%%%     Methods updating the panel display     %%%%%%%%%%%%%%%%

        function updateAperturePosScale(obj)
            xlim = obj.plotAxes.get('XLim');
            ylim = obj.plotAxes.get('YLim');
            obj.ApAxes.set('Position',...
                [(obj.u0-obj.radius-xlim(1))/(xlim(2)-xlim(1)),...
                (obj.v0-obj.radius-ylim(1))/(ylim(2)-ylim(1)),...
                2*obj.radius/(xlim(2)-xlim(1)), 2*obj.radius/(ylim(2)-ylim(1))]);
        end

        function updateApertureShape(obj)
            if(isempty(obj.Ap))
                xLim = [0 1];
                yLim = [0 1];
            else
                xLim = [0 size(obj.Ap,2)];
                yLim = [0 size(obj.Ap,1)];
            end
            
            obj.ApAxes.set('XLim',xLim);
            obj.ApAxes.set('YLim',yLim);
            obj.ApPos.set('XData',xLim(2)/2,'YData',yLim(2)/2);
            obj.ApImg.set('Cdata',bsxfun(@times,ones([size(obj.Ap),3]),permute(obj.ApColor,[3 1 2])));
            obj.ApImg.set('AlphaData',max(0,(obj.Ap*obj.ApAlpha)));
        end
        
        function centerView(obj)
            if(obj.hasUVcoords)
                scale = 3;
                obj.uPlotCenter = (max(obj.U(:))+min(obj.U(:)))/2;%mean(U(:));%
                obj.vPlotCenter = (max(obj.V(:))+min(obj.V(:)))/2;%mean(V(:));%
                obj.UVLim(1:2) = [obj.uPlotCenter - scale*(obj.uPlotCenter-min(obj.U(:))), obj.uPlotCenter + scale*(max(obj.U(:))-obj.uPlotCenter)];
                obj.UVLim(3:4) = [obj.vPlotCenter - scale*(obj.vPlotCenter-min(obj.V(:))), obj.vPlotCenter + scale*(max(obj.V(:))-obj.vPlotCenter)];
            else
                obj.uPlotCenter = 0;
                obj.vPlotCenter = 0;
            end
            obj.UpdatePlotLimits(1);
        end
        
        
        function UpdatePlotLimits(obj,zoomScale)
            %Change limits of the plot
            % - Force horizontal sale to remain equal to the vertical scale.
            % - Force the plot to remain centered on the point (uPlotCenter, vPlotCenter).
            pposPanel = getpixelposition(obj.panel);
            dV = (obj.UVLim(4)-obj.UVLim(3)) * zoomScale;
            dU = dV * pposPanel(3) / pposPanel(4);
            obj.UVLim(1) = obj.uPlotCenter - dU/2;
            obj.UVLim(2) = obj.uPlotCenter + dU/2;
            obj.UVLim(3) = obj.vPlotCenter - dV/2;
            obj.UVLim(4) = obj.vPlotCenter + dV/2;
            
            obj.plotAxes.set('XLim', obj.UVLim(1:2));
            obj.plotAxes.set('YLim', obj.UVLim(3:4));
            
            obj.plotAxes.set('XTick', round(obj.UVLim(1)):round(obj.UVLim(2)) );
            obj.plotAxes.set('YTick', round(obj.UVLim(3)):round(obj.UVLim(4)) );
            
            obj.updateAperturePosScale();
        end
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%     Access methods     %%%%%%%%%%%%%%%%%%%%%%%%%%

        function setViewPosition(obj,u0,v0)
            obj.u0 = u0;
            obj.v0 = v0;
            obj.updateAperturePosScale();
            obj.uPosEditor.set('String',num2str(round(u0,obj.ValboxDecimalDigits)));
            obj.vPosEditor.set('String',num2str(round(v0,obj.ValboxDecimalDigits)));
            pause(.00001);%drawnow%
        end
         
        function setRadius(obj, radius)
            obj.radius = radius;
            obj.updateAperturePosScale();
        end
        
        function setApertureShape(obj, Ap)
            obj.Ap = Ap;
            obj.updateApertureShape();
        end
        
        function setCallbackFcn(obj, callbackFcn)
            if(isa(callbackFcn, 'function_handle'))
                obj.callbackFcn = callbackFcn;
            else
                obj.callbackFcn=[];
            end
        end
        
% Implement abstract method from ListLayoutComponent class
        function setPixelPosition(obj,pos)
            setpixelposition(obj.panel,pos);
            obj.updateApertureShape();
            obj.updateAperturePosScale();
        end
        
        
    end % methods
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%     Callback functions     %%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Static )
        
        function panelResize(obj,src,event)
            obj.UpdatePlotLimits(1);
        end
        
        function mouseRelease(obj,src,event)
            obj.buttonClicked = false;
            obj.middleButtonClicked = false;
        end
        
        function mouseMotion(obj,src,event)
            if(obj.buttonClicked)
                P = get(obj.plotAxes, 'CurrentPoint');
                PSnap = round([P(1,1:2)]);
                if( sum((PSnap-P(1,1:2)).^2) < .02 )
                    P = PSnap;
                end
                obj.setViewPosition(P(1,1), P(1,2));
                obj.updateAperturePosScale();
                if(~isempty(obj.callbackFcn)), obj.callbackFcn(obj,'MovePosition');end
                
            elseif(obj.middleButtonClicked)
                P = get(obj.plotAxes, 'CurrentPoint');
                obj.uPlotCenter = obj.uPlotCenter + obj.uDragCurrent - P(1,1);
                obj.vPlotCenter = obj.vPlotCenter + obj.vDragCurrent - P(1,2);
                
                obj.UpdatePlotLimits(1);
                
            end
        end
        
        function mouseClick(obj,src,event)
            if(event.Button==2)
                obj.middleButtonClicked = true;
                P = get(obj.plotAxes, 'CurrentPoint');
                obj.uDragCurrent = P(1,1);
                obj.vDragCurrent = P(1,2);
            else
                P = get(obj.plotAxes, 'CurrentPoint');
                
                switch get(gcf,'SelectionType')
                  case 'open' %double-click
                    if(obj.hasUVcoords)
                        %search closest point to the clicked position in the input U,V coordinates.
                        [~,idClosest] = min((P(1,1)-obj.U(:)).^2 + (P(1,2)-obj.V(:)).^2);
                        P(1,1) = obj.U(idClosest);
                        P(1,2) = obj.V(idClosest);
                        if(P(1,1)<obj.UVLim(1)||P(1,1)>obj.UVLim(2)||P(1,2)<obj.UVLim(3)||P(1,2)>obj.UVLim(4)), obj.centerView();end
                    else
                        obj.buttonClicked = true;
                    end
                  otherwise
                    obj.buttonClicked = true;
                end
                
                obj.setViewPosition(P(1,1), P(1,2));
                obj.updateAperturePosScale();
                if(~isempty(obj.callbackFcn)), obj.callbackFcn(obj,'MovePosition');end
                
            end
        end
        
        function mouseScrollWheel(obj,src,event)
            P = get(obj.plotAxes, 'CurrentPoint');
            if( P(1,1) > obj.UVLim(1) && P(1,1) <= obj.UVLim(2) && P(1,2) > obj.UVLim(3) &&  P(1,2) <= obj.UVLim(4) )
                obj.UpdatePlotLimits(1+.1*event.VerticalScrollCount);
            end
        end
        
        function editUVPos(obj,src,event)
            val = str2double(src.get('String'));
            if(isnan(val))
                if(src == obj.uPosEditor)
                    src.set('String',num2str(round(obj.u0,obj.ValboxDecimalDigits)));
                else
                    src.set('String',num2str(round(obj.v0,obj.ValboxDecimalDigits)));
                end
            else
                if(src == obj.uPosEditor)
                    obj.setViewPosition( val, obj.v0 );
                    uicontrol(obj.vPosEditor);%set the focus to the other box.
                else
                    obj.setViewPosition( obj.u0, val );
                    uicontrol(obj.uPosEditor);%set the focus to the other box.
                end
                if(~isempty(obj.callbackFcn)), obj.callbackFcn(obj,'MovePosition');end
            end
        end
        
    end % methods ( Static )
    
end