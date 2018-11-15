classdef ListLayoutPanel < handle
    
    properties ( SetAccess = private )
        panel
        componentList = {}   
    end
    
    
    methods
        
        function obj = ListLayoutPanel(varargin)
            obj.panel = uipanel(varargin{:});
            obj.panel.set('ResizeFcn',@(src,evnt)ListLayoutPanel.panelResize(obj,src,evnt));
        end
        
        function addComponent(obj,component)
            if(isa(component,'ListLayoutComponent'))
                alreadyContained = false;
                for i=1:length(obj.componentList)
                    if(component == obj.componentList{i})
                        alreadyContained = true;
                        break;
                    end
                end
                if(~alreadyContained)
                    obj.componentList{end+1} = component;
                end
                
            end
        end
        
        function initLayout(obj)
            ListLayoutPanel.panelResize(obj,obj.panel,[]);
        end
        
    end
    
    methods( Static )
        
        function panelResize(obj,src,event)
            
            pposPanel = getpixelposition(src).*get(src,'InnerPosition')./get(src,'Position');
            
            
            %Compute remaining space to allocate after summing the minimum sizes and margins of all component.
            minTotSizeVertical = 0;
            for i=1:length(obj.componentList)
                minTotSizeVertical = minTotSizeVertical + obj.componentList{i}.minSize(2) + 1 * obj.componentList{i}.minMargin(2);
            end
            vFreeSpaceTotal = max( 0, ( pposPanel(4) - minTotSizeVertical ) );
            
            
            %Find how much of the remaining space to allocate to each component (avoids allocating a component more than its maximum size+margin).
            vFreeSpace = ones(1,length(obj.componentList)) * vFreeSpaceTotal / length(obj.componentList); %start with even allocation.
            reachedMaxSize = false(1,length(obj.componentList));
            vFreeSpaceUnused = 1;
            while(vFreeSpaceUnused > 0)
                vFreeSpaceUnused = 0;
                for i=1:length(obj.componentList)
                    minSize = obj.componentList{i}.minSize(2) + obj.componentList{i}.minMargin(2);
                    maxSize = obj.componentList{i}.maxSize(2) + obj.componentList{i}.maxMargin(2);
                    if(~reachedMaxSize(i) && vFreeSpace(i) > maxSize - minSize )
                        vFreeSpaceUnused = vFreeSpaceUnused + vFreeSpace(i) + minSize - maxSize;
                        vFreeSpace(i) = maxSize - minSize;
                        reachedMaxSize(i) = true;
                    end
                end
                for i=1:length(obj.componentList)
                    if(~reachedMaxSize(i))
                        vFreeSpace(i) = vFreeSpace(i) + vFreeSpaceUnused / sum(~reachedMaxSize);
                    end
                end
            end
            
            
            %Compute sizes and margins of all component.
            vPosCumulate = 0;
            for i=1:length(obj.componentList)
                
                hFreeSpace = max( 0, pposPanel(3) - ( obj.componentList{i}.minSize(1) + 2 * obj.componentList{i}.minMargin(1) ) );
                
                szX = min(obj.componentList{i}.maxSize(1), obj.componentList{i}.minSize(1) + hFreeSpace );
                marginFreeSpace = hFreeSpace - (szX - obj.componentList{i}.minSize(1));
                marginX = min(obj.componentList{i}.maxMargin(1), obj.componentList{i}.minMargin(1) + marginFreeSpace/2 );
                
                
                szY = min(obj.componentList{i}.maxSize(2), obj.componentList{i}.minSize(2) + vFreeSpace(i) );
                marginFreeSpace = vFreeSpace(i) - (szY - obj.componentList{i}.minSize(2));
                marginY = min(obj.componentList{i}.maxMargin(2), obj.componentList{i}.minMargin(2) + marginFreeSpace/1 );
                
                
                vPosCumulate = vPosCumulate + marginY;
                pos = [1+marginX, vPosCumulate, szX, szY];
                pos(2) = 1 + pposPanel(4) - pos(2) - szY; %invert top-bottom
                vPosCumulate = vPosCumulate + szY;
                
                pos = round(pos);
                obj.componentList{i}.setPixelPosition(pos);
            end
            
        end
        
    end
    
    
end