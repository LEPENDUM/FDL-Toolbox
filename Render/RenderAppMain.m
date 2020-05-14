% Main FDL rendering application :
% Create the graphical user interface.
%--------------------------------------------------------------------------
%
% Function Call:
%*****************************************************************
% 1-From the FDL model with parameters:
%*****************************************************************
%RenderAppMain(FDL, fullSize, crop, Disps, U,V, DispMap, gammaOptions);
% with input arguments:
% -FDL:
%	Fourier Disparity Layers (in the Fourier Domain) given as a complex array with dimensions 1.Vertical axis (Y), 2.Horizontal axis (X), 3.Color channels, 4.Layers.
%	Only left half of the spectrum (i.e. negative or null horizontal frequencies) may be given.
%
% -fullSize:
%	Full image size in the format [verticalResolution, horizontalResolution], including the padded borders.
%	The horizontal resolution may be different from the size of the 2nd dimension of the input FDL array if it only contains half of the spectrum.
%	
% -crop:
%	Number of pixels to crop (to remove the padded borders) on each side of the image : [left, right, top, bottom].
%	Alternatively, a BorderParams object can be given with the fields L,R,T,B respectively used for left, right, top and bottom crop.
%
% -Disps:
%	Vector of disparity values associated with the layers.
%	The number of elements must be the same as the number of layers (size of the 4th dimension of the input FDL array).
%
% -U,V (Optional, empty by default):
%	Lists of angular coordinates of the views from which the layer model was constructed (displayed in the viewpoint panel and needed to define useful range for some render parameters).
%	U is the list of horizontal angular coordinates, and V is the list of vertical angular coordinates.
%
% -DispMap (Optional, empty by default):
%	Disparity map (used by GUI application for refocusing automatically when the user clicks on the image).
%   If empty or not specified, a disparity map will be estimated.
%   Can be set to any scalar value to skip disparity estimation (automatic refocus on click won't work in this case).
%
% -gammaOptions (Optional, empty by default)
%	Cell array of the form {isLinear, gammaOffset} with:
%		isLinear (default=false): Set to true to apply gamma correction after rendering.
%		gammaOffset (default=0): Offset parameter to apply after gamma correction (only used if isLinear is true).
%
% -useGPU (Optional, true by default):
%	true-> use GPU acceleration if a CUDA Device is available and matlab's parallel processing toolbox is active.
%	false-> CPU only.
%
%*****************************************************************
% 2-From a file:
%*****************************************************************
%RenderAppMain(filename);
%
%*****************************************************************
% 3-Without input argument (opens a dialog box to load a fdl file)
%*****************************************************************
%RenderAppMain();
%
%
%--------------------------------------------------------------------------
%
% GUI Features:
% - Change viewpoint (--> from the viewpoint panel, or by dragging the image).
% - Change focus( --> click on the image or use the slider).
% - Change aperture radius, polygonal aperture shape, angle and thickness (--> sliders).
% - Show disparity map (--> 'd' key).
% - Show render in the Fourier domain (i.e. skip inverse transform) (--> space key).
% - Display rendered image in a separate figure (--> ctrl+space).
% - Save image or FDL model (--> menu bar).
%--------------------------------------------------------------------------
%
% See also RenderModel

function gui = RenderAppMain(FDL, fullSize, crop, Disps, U,V, DispMap, gammaOptions, useGPU)


    if(nargin==0)
        [file, path] = uigetfile('*.fdl');
        disp(['loading FDL : ''' fullfile(path,'/', file) '''']);
        load(fullfile(path,'/', file),'-mat');
    elseif(nargin==1)
        path=fileparts(FDL);
        load(FDL,'-mat');
    end

    if(~exist('U','var') || ~exist('V','var'))
        U=[]; V=[];
    end
    
    if(~exist('DispMap','var'))
        DispMap = [];
    end
    
    if(~exist('isLinear','var')), isLinear=[];end
    if(~exist('gammaOffset','var')), gammaOffset=[];end
    if(~exist('gammaOptions','var')), gammaOptions={};end
    if(~isempty(gammaOptions))
        isLinear = gammaOptions{1};
        if(length(gammaOptions)>1)
            gammaOffset = gammaOptions{2};
        end
    end
    if(~exist('useGPU','var')), useGPU=true;end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%
    rMod = RenderModel(FDL, fullSize, crop, Disps, DispMap, isLinear, gammaOffset, useGPU);
    gui = RenderGUI(rMod, Disps, U, V);
    if(exist('path','var')), gui.defaultSaveFolder = path;end
    
    set(gui.Fig,'Units','Pixel','OuterPosition',[-5 244 1182.5 837]);
    
    rMod.setPosition(gui.viewpointPanel.uPlotCenter, gui.viewpointPanel.vPlotCenter);
    
end