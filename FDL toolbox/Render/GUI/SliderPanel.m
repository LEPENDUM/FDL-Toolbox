classdef SliderPanel < ListLayoutComponent
    
    properties ( SetAccess = private )
        callbackFcn=[]
        value
    end
    
    properties ( Access = private )
        
        slideScale = 100
        ValboxDecimalDigits = 2;
        txtHeight = 20;
        txtMargin = 0;
        minSliderHeight = 10;
        maxSliderHeight = 30;
        sliderHeightSwitch
        totMin
        totMax
        
        
        panel
        text
        slider
        sliderHG
        ValEditor
        
        skipJSliderCallback = false
        
    end
    
    
    methods
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%     Constructor     %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = SliderPanel(Parent, callbackFcn, name, min, max, val, stepMinor, stepMajor, Position)
            obj@ListLayoutComponent(Parent);
            if(isa(Parent,'ListLayoutPanel'))
                ParentUI = Parent.panel;
            else
                ParentUI = Parent;
            end            
            
            
            obj.minSize(2) = obj.minSliderHeight + obj.txtHeight + obj.txtMargin;
            obj.maxSize(2) = obj.maxSliderHeight + obj.txtHeight + obj.txtMargin;
            obj.sliderHeightSwitch = [obj.maxSliderHeight,  17 ];

            ParentColor = get(ParentUI,'BackgroundColor');
            
            %Main matlab uipanel
            obj.panel = uipanel('Parent',ParentUI, 'Units','Pixels', 'BackgroundColor', ParentColor,'BorderType','none');
            obj.panel.set('ResizeFcn',@(src,evnt) SliderPanel.panelResize(obj,src,evnt));
            if(exist('Position','var') && ~isempty(Position)), obj.panel.set('Position',Position); end

            %Text component (slider name)
            obj.text = uicontrol('Parent', obj.panel, 'Style','text','String',name, 'BackgroundColor',ParentColor);
            obj.text.set('Fontsize',8.5);
            try
                obj.text.set('FontName', 'Microsoft JhengHei');
            catch
            end

            %Create slider component
            [obj.slider,obj.sliderHG] = javacomponent(javax.swing.JSlider,[], obj.panel);
            obj.slider.setPaintTicks(true);
            set(obj.slider, 'Minimum', min*obj.slideScale, 'Maximum', max*obj.slideScale);
            set(obj.slider, 'MinorTickSpacing',stepMinor*obj.slideScale, 'MajorTickSpacing',stepMajor*obj.slideScale, 'PaintLabels',true, 'Background',java.awt.Color(ParentColor(1),ParentColor(2),ParentColor(3)));


            %Create tick labels of the sllider
            labTable = java.util.Hashtable();
            for i = min:stepMajor:max
                lbl = javax.swing.JLabel(num2str(i));
                lbl.setFont(lbl.getFont().deriveFont(16));
                labTable.put(java.lang.Integer(int32(i*obj.slideScale)), lbl);
            end
            obj.slider.setLabelTable(labTable);
            
            %Create value editable box
            obj.ValEditor = uicontrol('Parent', obj.panel, 'Style','edit', 'BackgroundColor',ParentColor);
            
            %Initialize slider's value
            obj.setValue(val);
            
            %Set the callback function of the SliderPanel object
            obj.setCallbackFcn(callbackFcn);
        end
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%     Access methods     %%%%%%%%%%%%%%%%%%%%%%%%%%
        function setValue(obj, value)
            if(value ~= obj.slider.value)
                %Set the slider's value without performing the actions in callbackFromJSlider (these actions are only required when the change of value is operated directly from slider's action).
                obj.skipJSliderCallback = true;
                set(obj.slider, 'Value', max(min(value*obj.slideScale, get(obj.slider,'Maximum')),get(obj.slider,'Minimum')) );
            end
            
            set(obj.ValEditor,'String',num2str(round(value,obj.ValboxDecimalDigits)));
            obj.value = value;
        end

        
        function setCallbackFcn(obj, callbackFcn)
            if(isa(callbackFcn, 'function_handle'))
                obj.callbackFcn = callbackFcn;
                set(obj.slider, 'StateChangedCallback', @(src,evnt) SliderPanel.callbackFromJSlider(obj,src,evnt));
                set(obj.ValEditor, 'Callback', @(src,evnt) SliderPanel.callbackFromEditBox(obj,src,evnt));
            else
                obj.callbackFcn=[];
                set(obj.slider, 'StateChangedCallback', []);
                set(obj.ValEditor, 'Callback', []);
            end
        end
        
        %{
        function setPrecisionScale(obj,scale)
            
            min = get(obj.slider, 'Minimum') / obj.slideScale;
            max = get(obj.slider, 'Maximum') / obj.slideScale;
            minorTickSpacing = get(obj.slider, 'MinorTickSpacing') / obj.slideScale;
            majorTickSpacing = get(obj.slider, 'MajorTickSpacing') / obj.slideScale;
            set(obj.slider, 'MinorTickSpacing',minorTickSpacing*scale, 'MajorTickSpacing',majorTickSpacing*scale);
            set(obj.slider, 'Minimum', min*scale, 'Maximum', max*scale);
            
            obj.slideScale = scale;
            
            labTable = java.util.Hashtable();
            for i = min:majorTickSpacing:max
                lbl = javax.swing.JLabel(num2str(i));
                lbl.setFont(lbl.getFont().deriveFont(16));
                labTable.put(java.lang.Integer(int32(i*obj.slideScale)), lbl);
            end
            obj.slider.setLabelTable(labTable);
        end
        %}

% Implement abstract method from ListLayoutComponent class
        function setPixelPosition(obj,pos)
            setpixelposition(obj.panel,pos);
        end
        

    end % methods
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%     Callback functions     %%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Static )
        
        %Resizing function
        function panelResize(obj,src,event)
            pposPanel = getpixelposition(src);
            
            pposSlider = [1, 1, pposPanel(3)*.83, pposPanel(4)-obj.txtHeight-obj.txtMargin];
            pposText = [1, 1+pposSlider(4)+obj.txtMargin, pposPanel(3), obj.txtHeight];
            pposValEditor = [1+pposPanel(3)*.85, 1+pposPanel(4)*.2, pposPanel(3)*.15, pposPanel(4)*.5];

            if(pposSlider(4)<obj.sliderHeightSwitch(1)), set(obj.slider,'PaintLabels',false);else, set(obj.slider,'PaintLabels',true);end
            if(pposSlider(4)<obj.sliderHeightSwitch(2)), set(obj.slider,'PaintTicks',false);else, set(obj.slider,'PaintTicks',true);end

            setpixelposition(obj.text, pposText);
            setpixelposition(obj.ValEditor, pposValEditor);
            setpixelposition(obj.sliderHG, pposSlider);

        end
        
        function callbackFromJSlider(obj,src,event)
            if(~obj.skipJSliderCallback) %only true when the action comes from the slider (to allow values beyond the displayable slider range without interfering with other components).
                obj.value = obj.slider.value/obj.slideScale;
                obj.ValEditor.set('String',num2str(obj.value));
                if(~isempty(obj.callbackFcn)), obj.callbackFcn(obj,event); end
            end
            obj.skipJSliderCallback = false;
        end
        
        function callbackFromEditBox(obj,src,event)
            val = str2double(obj.ValEditor.get('String'));
            if(isnan(val))
                obj.ValEditor.set('String',num2str(obj.value));
            else
                obj.setValue(val);
                if(~isempty(obj.callbackFcn)), obj.callbackFcn(obj,event); end
            end
        end

    end % methods ( Static )

end