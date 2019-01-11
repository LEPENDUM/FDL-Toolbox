% Main FDL rendering application :
% Create the graphical user interface.
%
% Usage:
% 1-From the FDL model with parameters:
%RenderAppMain(FDL, fullSize, crop, Disps, U,V, DispMap, gammaOptions);
%   -FDL   : Fourier Disparity layers (only left half of the spectrum (i.e. negative or null horizontal frequencies) may be given).
%   -fullSize : Full image size in the format [verticalResolution, horizontalResolution], including the padded borders.
%   -crop  : Number of pixels to crop (to remove the padded borders) on each side of the image : [top, bottom, left, right].
%   -Disps : List of disparity values associated with each layer.
%   -U,V   : (Optional) Lists of angular coordinates of the views from which the layer model was constructed (displayed in the viewpoint panel and needed to define useful range for some render parameters).
%   -DispMap: (Optional) Disparity map used for refocusing automatically when the user clicks on the image (if not given, a disparity map will be estimated).
%   -gammaOptions: (Optional) Cell array of the form {isLinear, gammaOffset}. Use isLinear=true if gammaCorrection is needed / gammaOffset: offset parameter to apply after gamma correction (only used if isLinear is true).
%   -useGPU : (Optional) true(default): use GPU acceleration if a CUDA Device is available and matlab's parallel processing toolbox is active. / false: CPU only.
%
% 2-From a file:
%RenderAppMain('path/filename.fdl');
%
% 3-Without input argument (opens a dialog box to load a fdl file).
%RenderAppMain();
%
%
% Features:
% - Change viewpoint (--> from the viewpoint panel, or by dragging the image).
% - Change focus( --> click on the image or use the slider).
% - Change aperture radius, polygonal aperture shape, angle and thickness (--> sliders).
% - Show disparity map (--> 'd' key).
% - Show render in the Fourier domain (i.e. skip inverse transform) (--> space key).
% - Display rendered image in a separate figure (--> ctrl+space).
% - Save image or FDL model (--> menu bar).

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