classdef EmptyListLayoutComponent < ListLayoutComponent
    methods
        %constructor inherited from ListLayoutComponent with parameters:
        % -listLayoutPanel
        % -minSizeHorizontal 
        % -maxSizeHorizontal
        % -minSizeVertical
        % -maxSizeVertical
        % -minMarginLeft
        % -maxMarginLeft
        % -minMarginTop
        % -maxMarginTop
        % -minMarginRight
        % -maxMarginRight
        % -minMarginBottom
        % -maxMarginBottom
        function obj=EmptyListLayoutComponent(varargin)
            obj@ListLayoutComponent(varargin{:});
        end
        
% Implement abstract method from ListLayoutComponent class
        function setPixelPosition(obj,pos)
        end
    end
end