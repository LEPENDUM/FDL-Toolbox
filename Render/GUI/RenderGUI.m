%Graphical user interface for the light field rendering application
classdef RenderGUI < handle
    
    properties
        model
        
        U = []
        V = []
        
        Fig
        UIPanel
        viewpointPanel
        focusPanel
        radiusPanel
        numBladesPanel
        thicknessPanel
        anglePanel
        
        imAxes
        imHandle
        
        ImgClicked = false;
        ImgMotion = false;
        
        displayMode=0 %0=show rendered image / 1=show disparity map.
        
        defaultSaveFolder = ''
        defaultSaveImage = ''
        
    end
    
    methods
        function obj = RenderGUI(model, D, U, V)
            
            if(~exist('model','var') || isempty(model))
                emptyModel = RenderModel();
                obj.model = emptyModel;
                obj.createGUI([],[],[-1,1]);
            else
                
                obj.model = model;
                obj.model.setApShape('polygon');

                addlistener(obj.model,'ChangeApShape', @(src,evnt)RenderGUI.changedApertureShape(obj, src, evnt));
                addlistener(obj.model,'ChangeNumBlades', @(src,evnt)RenderGUI.changedApertureShape(obj, src, evnt));
                addlistener(obj.model,'ChangeApThickness', @(src,evnt)RenderGUI.changedApertureShape(obj, src, evnt));
                addlistener(obj.model,'ChangeApAngle', @(src,evnt)RenderGUI.changedApertureShape(obj, src, evnt));

                addlistener(obj.model,'ChangeRadius', @(src,evnt)RenderGUI.changedRenderParam(obj, src, evnt));
                addlistener(obj.model,'ChangeFocus', @(src,evnt)RenderGUI.changedRenderParam(obj, src, evnt));
                addlistener(obj.model,'ChangePosition', @(src,evnt)RenderGUI.changedRenderParam(obj, src, evnt));
                addlistener(obj.model,'ToggleSkipInvTr', @(src,evnt)RenderGUI.changedRenderParam(obj, src, evnt));


                if(~exist('U','var') || ~exist('V','var'))
                    U=[]; V=[];
                end

                obj.U=U;
                obj.V=V;
                obj.createGUI(U,V,D);
            end
            
        end
        
        
        function createGUI(obj, U, V, D)
            obj.Fig = figure('Color','k');
            obj.Fig.set('WindowButtonMotionFcn',@RenderGUI.voidFcn); %tip: force update of the currentPoint property of all UI elements before triggering WindowMouseMotion event.
            
            %Setup save toolbar buttons
            hSave = findall(obj.Fig, 'tag', 'figMenuFileSave');
            hSavebutton = findall(obj.Fig, 'tag','Standard.SaveFigure');
            hSaveFig = findall(obj.Fig, 'tag', 'figMenuFileSaveAs');
            if(~isempty(hSave))
                set(hSave,'Label','&Save Image...', 'Callback',@(src,evnt)RenderGUI.saveImageCallback(obj,src,evnt));
            end
            if(~isempty(hSavebutton))
                hSavebutton.TooltipString = 'Save Image';
                set(hSavebutton,'ClickedCallback',@(src,evnt)RenderGUI.saveImageCallback(obj,src,evnt));
            end
            if(~isempty(hSaveFig))
                set(hSaveFig,'Label','&Save FDL...', 'Callback',@(src,evnt)RenderGUI.saveAsCallback(obj,src,evnt));
            end
            
            %Callbacks for the main GUI
            obj.Fig.set('KeyPressFcn',@(src,evnt)RenderGUI.KeyPress(obj, src, evnt));
            obj.Fig.set('KeyReleaseFcn',@(src,evnt)RenderGUI.KeyRelease(obj, src, evnt));
            addlistener(obj.Fig,'WindowMouseMotion', @(src,evnt)RenderGUI.mouseMotion(obj,src,evnt));
            addlistener(obj.Fig,'WindowMouseRelease', @(src,evnt)RenderGUI.mouseRelease(obj,src,evnt));
            addlistener(obj.Fig,'WindowScrollWheel', @(src,evnt)RenderGUI.scrollWheel(obj,src,evnt));
            
            
            %Main control panel
            UIPanelSize = .3;
            UIPanelColor = [.98 .98 .98];
            obj.UIPanel = ListLayoutPanel('Position',[.0 .0 UIPanelSize 1], 'units','normalized');
            obj.UIPanel.panel.set('BackgroundColor',UIPanelColor);
            
            %Image area
            obj.imAxes = axes('pos',[UIPanelSize+.01 .01 .98-UIPanelSize .98]);
            obj.imHandle = imshow(obj.model.Image);
            obj.imHandle.set('ButtonDownFcn',@(src,evnt)RenderGUI.ImgButtonDown(obj, src, evnt));

            
            vertMinMargin = 5;
            vertMaxMargin = 15;
            horizMinMargin = 10;
            horizMaxSize = 500;
            horizMinSize = 100;
            
            %Add empty Component to leave space at the top of the panel.
            EmptyListLayoutComponent(obj.UIPanel, 0, 0, 5, 50, 0, 0, 0, 0);
            
            %Add panel controling the viewpoint
            obj.viewpointPanel = ViewPointPanel(obj.UIPanel, obj.model.u0, obj.model.v0, obj.model.radius, obj.model.Ap, U, V);
            obj.viewpointPanel.setCallbackFcn(@(src,evnt)RenderGUI.moveViewPoint(obj, src, evnt));
            obj.viewpointPanel.minMargin = [horizMinMargin-5 0];
            obj.viewpointPanel.maxMargin(2) = 0;
            obj.viewpointPanel.minSize(1) = horizMinSize;
            obj.viewpointPanel.maxSize = [horizMaxSize+10, 400];
            
            %Setup focus slider range and tick spacing
            minFocus = floor(min(D(:))*(1+1.5)/2 + max(D(:))*(1-1.5)/2);
            maxFocus = ceil(max(D(:))*(1+1.5)/2 + min(D(:))*(1-1.5)/2);
            FocusStep = (maxFocus-minFocus)/4;
            radiusMinorStep = FocusStep/5;
            
            %Add slider controling the focus.
            obj.focusPanel = SliderPanel(obj.UIPanel, @(src,evnt)RenderGUI.slideFocus(obj, src, evnt), 'Refocus', minFocus, maxFocus, obj.model.s, radiusMinorStep,FocusStep);
            obj.focusPanel.minMargin = [horizMinMargin vertMinMargin+25];
            obj.focusPanel.maxMargin(2) = vertMaxMargin + 30;
            obj.focusPanel.minSize(1) = horizMinSize;
            obj.focusPanel.maxSize(1) = horizMaxSize;
            
            %Setup radius slider range and tick spacing.
            if(numel(U) > 1 && numel(V) == numel(U))
                radiusSliderMax = max( sqrt( ( U(:)-mean(U(:)) ).^2 + ( V(:) - mean(V(:)) ).^2 ) );
                radiusSliderMax = ceil(3*radiusSliderMax);
                radiusMinorStep = radiusSliderMax / 20;
                radiusStep = radiusSliderMax / 4;
            else
                radiusSliderMax = 20;
                radiusMinorStep = 1;
                radiusStep = 5;
            end
            
            %Add slider controling the aperture radius.
            obj.radiusPanel = SliderPanel(obj.UIPanel, @(src,evnt)RenderGUI.slideRadius(obj, src, evnt), 'Aperture radius', 0, radiusSliderMax, obj.model.radius, radiusMinorStep,radiusStep);
            obj.radiusPanel.minMargin = [horizMinMargin vertMinMargin];
            obj.radiusPanel.maxMargin(2) = vertMaxMargin;
            obj.radiusPanel.minSize(1) = horizMinSize;
            obj.radiusPanel.maxSize(1) = horizMaxSize;
            
            %Add slider controling the number of blades of the aperture.
            obj.numBladesPanel = SliderPanel(obj.UIPanel, @(src,evnt)RenderGUI.slideNumBlades(obj, src, evnt), 'Number of blades', 3, 15, obj.model.numBlades, 1,3);
            obj.numBladesPanel.minMargin = [horizMinMargin vertMinMargin + 25];
            obj.numBladesPanel.maxMargin(2) = vertMaxMargin + 40;
            obj.numBladesPanel.minSize(1) = horizMinSize;
            obj.numBladesPanel.maxSize(1) = horizMaxSize;
            
            %Add slider controling the aperture angle.
            obj.anglePanel = SliderPanel(obj.UIPanel, @(src,evnt)RenderGUI.slideAngle(obj, src, evnt), 'Aperture angle', 0, 180, obj.model.apAngle, 30,60);
            obj.anglePanel.minMargin = [horizMinMargin vertMinMargin];
            obj.anglePanel.maxMargin(2) = vertMaxMargin;
            obj.anglePanel.minSize(1) = horizMinSize;
            obj.anglePanel.maxSize(1) = horizMaxSize;
            
            %Add slider controling the aperture thickness.
            obj.thicknessPanel = SliderPanel(obj.UIPanel, @(src,evnt)RenderGUI.slideThickness(obj, src, evnt), 'Aperture thickness', 0, 1, obj.model.apThickness, .2,.2);
            obj.thicknessPanel.minMargin = [horizMinMargin vertMinMargin];
            obj.thicknessPanel.maxMargin(2) = vertMaxMargin;
            obj.thicknessPanel.minSize(1) = horizMinSize;
            obj.thicknessPanel.maxSize(1) = horizMaxSize;            
            
            %Add empty Component to leave space at the bottom of the panel.
            EmptyListLayoutComponent(obj.UIPanel, 0, 0, 0, 50, 0, 0, 0, 0);
            
            %initialize the ListLayoutPanel (compute sizes of all the added components).
            obj.UIPanel.initLayout();
            
        end % createGUI function
        
        
        function refreshApertureShape(obj)
            %recompute Aperture shape from modified model properties.
            obj.model.computeAperture();
            %display new aperture shape on screen.
            obj.viewpointPanel.setApertureShape(obj.model.Ap);
        end
        
        function refreshRender(obj)
            if(obj.displayMode==0)
                %recompute image from modified modified model properties.
                obj.model.renderImage();
                %display image on screen.
                set(obj.imHandle,'Cdata', obj.model.Image);
            elseif(obj.displayMode==1)
                set(obj.imHandle,'Cdata', uint8(repmat( 255*(obj.model.DispMap-min(obj.model.DispMap(:))) / (max(obj.model.DispMap(:))-min(obj.model.DispMap(:))), 1,1,3)));
            end
            
            pause(.02);
        end
        
    end %methods
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Callbacks functions                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Static )


%%%%%%%%%%%%%%%%%%%     Callbacks for model change      %%%%%%%%%%%%%%%%%%%
function changedApertureShape(obj, src, evnt)
    switch evnt.EventName
        case 'ChangeApShape'
            obj.viewpointPanel.setRadius(src.trueRadius);
        case 'ChangeNumBlades'
            obj.numBladesPanel.setValue(src.numBlades);
        case 'ChangeApThickness'
            obj.thicknessPanel.setValue(src.apThickness);
        case 'ChangeApAngle'
            obj.anglePanel.setValue( src.apAngle * 180/pi);
    end
    
    obj.refreshApertureShape();
    obj.refreshRender();
end

function changedRenderParam(obj, src, evnt)
    
    switch evnt.EventName
        
        case 'ChangeRadius'
            obj.radiusPanel.setValue(src.radius);
            obj.viewpointPanel.setRadius(src.trueRadius);
        case 'ChangeFocus'
            obj.focusPanel.setValue(src.s);
        case 'ChangePosition'
            obj.viewpointPanel.setViewPosition(src.u0, src.v0);
        case 'ToggleSkipInvTr'
    end
    obj.refreshRender();
end


%%%%%%%%%%%%%%%%%%%%%%%     Mouse and Keyboard     %%%%%%%%%%%%%%%%%%%%%%%%

function KeyPress(obj,src,event)
    switch event.Key
        case 'rightarrow'
            obj.model.setFocus(obj.model.s + .1);
        case 'leftarrow'
            obj.model.setFocus(obj.model.s - .1);
        case 'uparrow'
            obj.model.setRadius(obj.model.radius + .5);
        case 'downarrow'
            obj.model.setRadius(obj.model.radius - .5);
        case 'space'
             if(strcmp(event.Modifier,'control'))
                figure,imshow(obj.model.Image); %ctrl+space -> shows the image in a new figure
            else
                obj.model.setSkipInvTr(true);   %space -> skip inverse Fourier transform (show directly the result in the Fourier domain).
            end
        case 'd' %show disparity map
            if(obj.displayMode~=1)
                obj.displayMode = 1;
                obj.refreshRender();
            end
    end
    
    %Not necessary anymore : the ring, disk and square shapes can all be approximated with the polygon shape.
    %{
    if(strfind(event.Key,'numpad'))
        obj.model.setApShapeId(str2double(event.Character));
    end
    %}
end

function KeyRelease(obj,src,event)
    switch event.Key
        case 'space'
            obj.model.setSkipInvTr(false);
        case 'd'
            obj.displayMode = 0;
            obj.refreshRender();
    end
end

function ImgButtonDown(obj,src,event)
    obj.ImgClicked=true;
end


function scrollWheel(obj,src,event)
    
    P = get(obj.imAxes, 'CurrentPoint');
    if( P(1,1) > obj.imAxes.XLim(1) && P(1,1) <= obj.imAxes.XLim(2) && P(1,2) > obj.imAxes.YLim(1) &&  P(1,2) <= obj.imAxes.YLim(2) )
        if(strcmp(RenderModel.ApShapes{obj.model.ApShapeId},'ring'))
            obj.model.setApThickness(obj.model.apThickness + 5*event.VerticalScrollCount);
        elseif(strcmp(RenderModel.ApShapes{obj.model.ApShapeId},'polygon'))
            obj.model.setNumBlades(min(20, obj.model.numBlades - event.VerticalScrollCount));
        end
    end
end

function mouseRelease(obj,src,event)
    if(obj.ImgClicked)
        obj.ImgClicked = false;
        if(~obj.ImgMotion)
            P = get(obj.imAxes, 'CurrentPoint');
            x = round(P(1,1));
            y = round(P(1,2));
            if(x>0 && x<=size(obj.model.DispMap,2) && y>0 && y<=size(obj.model.DispMap,1))
                %x=min(max(round(P(1,1)),1),size(DispMap,1));
                %y=min(max(round(P(1,2)),1),size(DispMap,2));
                sAnimList = linspace(obj.model.s, obj.model.DispMap(y,x), 5);
                for s=sAnimList(2:end)
                    obj.model.setFocus(s);
                end
            end
        else
            obj.ImgMotion=false;
        end        
    end
end

function mouseMotion(obj,src,event)
    if(obj.ImgClicked)
        obj.ImgMotion=true;
        P = get(obj.imAxes, 'CurrentPoint');
        xLim = get(obj.imAxes,'XLim');
        yLim = get(obj.imAxes,'YLim');
        u0 = obj.viewpointPanel.uPlotCenter + ((P(1,1) - xLim(1))/(xLim(2)-xLim(1)) -.5)*8;
        v0 = obj.viewpointPanel.vPlotCenter - ((P(1,2) - yLim(1))/(yLim(2)-yLim(1)) -.5)*8; %minus sign to compensate for the reverse Ydir property of the image Axes (imposed by matlab for showing the image).
        obj.model.setPosition(u0,v0);
    end
    
end



function saveImageCallback(obj,src,event)
    if(isempty(obj.defaultSaveImage))
        [file,path] = uiputfile([obj.defaultSaveFolder '/*.png'], 'Save Image');
    else
        [file,path] = uiputfile(fullfile(obj.defaultSaveFolder,'/', obj.defaultSaveImage), 'Save Image');
    end
    if(path)
        obj.defaultSaveFolder = path;
        obj.defaultSaveImage = file;
        filename = fullfile(path,'/', file);
        imwrite(obj.model.Image,filename);
    end
end

function saveAsCallback(obj,src,event)
    [file,path] = uiputfile([obj.defaultSaveFolder '/*.fdl'], 'Save FDL');
    if(path)
        obj.defaultSaveFolder = path;
        filename = fullfile(path,'/', file);
        obj.model.saveFDL(filename,obj.U,obj.V);
    end
end

function voidFcn(obj,src)
end


%%%%%%%%%%%%%%%%%%%      Callbacks for UI elements      %%%%%%%%%%%%%%%%%%%

function moveViewPoint(obj, src, event)
    obj.model.setPosition( src.u0, src.v0 );
end

function slideRadius(obj, src, event)
    obj.model.setRadius( src.value );
end

function slideFocus(obj,src,event)
    obj.model.setFocus( src.value );
end

function slideNumBlades(obj,src,event)
    obj.model.setNumBlades( src.value );
end

function slideThickness(obj,src,event)
    obj.model.setApThickness( src.value );
end

function slideAngle(obj,src,event)
    obj.model.setApAngle( src.value*pi/180 );
end

    end %method ( Static )
    
end %classdef
