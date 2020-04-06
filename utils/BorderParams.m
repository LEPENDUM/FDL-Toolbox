% Class of parameters for padding/cropping and windowing of padded borders. The fields are:
%   - L: padding/crop on the left side (default=0).
%   - R: padding/crop on the right side (default=0).
%   - T: padding/crop on the top (default=0).
%   - B: padding/crop on the bottom (default=0).
%   - window:
%       - 'none' -> no windowing of padded borders.
%       - 'hann' (default) -> use hann window type for padded borders.
%       - 'linear' -> use linear window type for padded borders.
%   - padVal: padding value: either numeric value or 'replicate', 'symmetric' or 'circular' (default='symmetric').
%
% A BorderParam object is constructed with the following optional parameters using either a (name, value) pair input or an input structure with field names:
%   - L, R, T, B, window, padVal.
%   - X : sets L and R, if they are not specified individually.
%   - Y : sets T and B, if they are not specified individually.
%   - All : sets L, R, T, and B parameters, if they are not specidied individually or with X, Y.
%
% Construction examples:
%   %Construct a BorderParams object with 10 pixels on the left and right, 15 pixels on top and 12 pixels on bottom
%   bp = BorderParams('X',10,'T',15,'B',12);
%   %Same result is obtained using a structure as input of the constructor
%   bpStruct = [];
%   bpStruct.X=10;
%   bpStruct.T=15;
%   bpStruct.B=12;
%   bp = BorderParams(bpStruct);

classdef BorderParams
    properties
        L=0
        R=0
        T=0
        B=0
        window
        padVal
    end
    
    methods
        function obj = BorderParams(varargin)
            p = inputParser;
            validWindow = {'none','hann','linear'};
            checkWindow = @(x) any(validatestring(x,validWindow));
            validPadType = {'symmetric','circular','replicate'};
            checkPadVal = @(x) (isscalar(x) && isnumeric(x)) || any(validatestring(x,validPadType));
            
            addParameter (p,'All',0,@isnumeric);
            addParameter (p,'X',  0,@isnumeric);
            addParameter (p,'Y',  0,@isnumeric);
            addParameter (p,'L',  0,@isnumeric);
            addParameter (p,'R',  0,@isnumeric);
            addParameter (p,'T',  0,@isnumeric);
            addParameter (p,'B',  0,@isnumeric);
            addParameter (p,'window','hann',checkWindow);
            addParameter (p,'padVal','symmetric',checkPadVal);
            
            parse(p,varargin{:});
            
            if(~ismember('All',p.UsingDefaults))
                obj.L = p.Results.All;
                obj.R = p.Results.All;
                obj.T = p.Results.All;
                obj.B = p.Results.All;
            end
            if(~ismember('X',p.UsingDefaults))
                obj.L = p.Results.X;
                obj.R = p.Results.X;
            end
            if(~ismember('Y',p.UsingDefaults))
                obj.T = p.Results.Y;
                obj.B = p.Results.Y;
            end
            
            if(~ismember('L',p.UsingDefaults)), obj.L = p.Results.L; end
            if(~ismember('R',p.UsingDefaults)), obj.R = p.Results.R; end
            if(~ismember('T',p.UsingDefaults)), obj.T = p.Results.T; end
            if(~ismember('B',p.UsingDefaults)), obj.B = p.Results.B; end
            
            obj.window = p.Results.window;
            obj.padVal = p.Results.padVal;
            
        end
    end
    
end