classdef ListLayoutComponent < handle
    
    properties
        
        minSize = [50,50]     %vector of length 2 (horizontal, vertical).
        maxSize = [inf,inf]   %vector of length 2 (horizontal, vertical).
        minMargin = [0,0]     %vector of length 2 (horizontal, vertical).
        maxMargin = [inf,inf] %vector of length 2 (horizontal, vertical).
        
    end
    
    methods
        
        %function obj = ListLayoutComponent( listLayoutPanel, minSizeHorizontal, maxSizeHorizontal, minSizeVertical, maxSizeVertical, minMarginLeft, maxMarginLeft, minMarginTop, maxMarginTop, minMarginRight, maxMarginRight, minMarginBottom, maxMarginBottom )
        function obj = ListLayoutComponent( listLayoutPanel, minSizeHorizontal, maxSizeHorizontal, minSizeVertical, maxSizeVertical, minMarginHorizontal, maxMarginHorizontal, minMarginVertical, maxMarginVertical )
            if(exist('minSizeHorizontal','var')   && isscalar(minSizeHorizontal)   && isnumeric(minSizeHorizontal)),   obj.minSize(1) = minSizeHorizontal;end
            if(exist('maxSizeHorizontal','var')   && isscalar(maxSizeHorizontal)   && isnumeric(maxSizeHorizontal)),   obj.maxSize(1) = maxSizeHorizontal;end
            
            if(exist('minSizeVertical','var')     && isscalar(minSizeVertical)     && isnumeric(minSizeVertical)),     obj.minSize(2) = minSizeVertical;end
            if(exist('maxSizeVertical','var')     && isscalar(maxSizeVertical)     && isnumeric(maxSizeVertical)),     obj.maxSize(2) = maxSizeVertical;end
            
            if(exist('minMarginHorizontal','var') && isscalar(minMarginHorizontal) && isnumeric(minMarginHorizontal)), obj.minMargin(1) = minMarginHorizontal;end
            if(exist('maxMarginHorizontal','var') && isscalar(maxMarginHorizontal) && isnumeric(maxMarginHorizontal)), obj.maxMargin(1) = maxMarginHorizontal;end
            
            if(exist('minMarginVertical','var')   && isscalar(minMarginVertical)   && isnumeric(minMarginVertical)),   obj.minMargin(2) = minMarginVertical;end
            if(exist('maxMarginVertical','var')   && isscalar(maxMarginVertical)   && isnumeric(maxMarginVertical)),   obj.maxMargin(2) = maxMarginVertical;end
            
            %if(exist('minMarginRight','var')   && isscalar(minMarginRight)   && isnumeric(minMarginRight)),   obj.minMargin(3) = minMarginRight;end
            %if(exist('maxMarginRight','var')   && isscalar(maxMarginRight)   && isnumeric(maxMarginRight)),   obj.maxMargin(3) = maxMarginRight;end
            
            %if(exist('minMarginBottom','var')   && isscalar(minMarginBottom)   && isnumeric(minMarginBottom)),   obj.minMargin(4) = minMarginBottom;end
            %if(exist('maxMarginBottom','var')   && isscalar(maxMarginBottom)   && isnumeric(maxMarginBottom)),   obj.maxMargin(4) = maxMarginBottom;end
                        
            if(exist('listLayoutPanel','var')), obj.addToListLayout(listLayoutPanel); end
            
        end
        
        function addToListLayout(obj, listLayoutPanel)
            if(isa(listLayoutPanel,'ListLayoutPanel'))
                listLayoutPanel.addComponent(obj);
            end
        end
        
    end
    
    methods( Abstract )
        setPixelPosition(obj,pos)
    end
    
end