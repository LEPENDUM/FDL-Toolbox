%Main class for light field rendering from the FDL model.
%--------------------------------------------------------------------------
%
% The RenderModel object constructor has the following arguments:
%
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
% -Disps :
%	Vector of disparity values associated with the layers.
%	The number of elements must be the same as the number of layers (size of the 4th dimension of the input FDL array).
%
% -DispMap (Optional, empty by default):
%	Disparity map (used by GUI application for refocusing automatically when the user clicks on the image).
%   If empty or not specified, a disparity map will be estimated.
%   Can be set to any scalar value to skip disparity estimation.
%
% -isLinear (Optional, false by default):
%	Set to true to apply gamma correction after rendering.
%
% -gammaOffset (Optional, 0 by default
%	Offset parameter to apply after gamma correction (only used if isLinear is true).
%
% -useGPU (Optional, true by default):
%	true-> use GPU acceleration if a CUDA Device is available and matlab's parallel processing toolbox is active.
%	false-> CPU only.
%
%--------------------------------------------------------------------------
%
% The methods available in the RenderModel class are:
%
% - obj.renderImage()
%		Renders the image with the current RenderModel's parameters (e.g. focus, aperture parameters).
%		The rendered image can be accessed either with the property obj.Image (uint8 format with cropped borders), or with the method obj.getInternalImage (see method documentation).
%
% - obj.computeAperture()
%		Updates the internal representation of the aperture with the current aperture parameters (e.g. shape type, thickness, number of blades for polygon shape).
%		This method must be called before rendering the image (i.e. call to obj.renderImage) if aperture parameters are changed (i.e. call to the methods with names starting with setAp) to take the changes into account in the final render.
%
% - obj.setSkipInvTr(skipInvTr)
%		Set the skipInvTr property to true or false:
%			-true-> Skip the inverse Fourier transform when rendering the image (the rendered image in obj.Image gives a visual representation of the Fourier magnitude spectrum).
%			-false-> Apply inverse Fourier Transform when rendering the image.
%		skipInvTr is set to false by default when constructing a RenderModel object.
%
% - obj.getInternalImage(applyCrop)
%		Outputs the rendered image used internally (in floating point format).
%		The input argument applyCrop (false by default) can be set to true to apply the cropping operation on the output image (defined by the crop property of the RenderModel object).
%		If the skipInvTr property is true, the ouput image is the complex Fourier domain representation of the rendered image (and no cropping is applied even if applyCrop is set to true).
%
% - obj.saveFDL(filename,U,V)
%		Save the FDL model (layers with disparity values and dimension variables) in the file given in filename argument (uses the matlab mat file system).
%		
%		Optionally, a list of angular coordinates defined by the input vectors U (horizontal coordinates) and V (vertical coordinates) can be saved in the same file (e.g. coordinates of the views used for the FDL construction).
%
% - obj.setRadius(radius)
%		Set the aperture radius of the image to render (relative to the scale of the camera plane defined by the angular coordinates used as parameters for FDL construction).
%		A value of 0 can be used to render images with full depth of field (simulate pinhole aperture).
%		The meaning of the radius depends on the aperture type (see below for a description of aperture types).
%		The focus parameter can be accessed for reading with the property obj.radius (set to 0 for a default RenderModel object).
%
% - obj.setFocus(focus)
%		Set the focus parameter of the image to render (relative to the disparity values of the FDL model).
%		The focus parameter can be accessed for reading with the property obj.s (set to 0 for a default RenderModel object).
%
% - obj.setPosition(u,v)
%		Set the angular coordinates of the image to render (relative to the scale of the camera plane defined by the angular coordinates used as parameters for FDL construction).
%		The angular coordinates can be accessed for reading with the properties obj.u0 and obj.v0 (both set to 0 for a default RenderModel object).
%
% - obj.setApShape(ApShapeName)
%		Set the aperture shape type by its name.
%		See below for a description of aperture types and corresponding names.
%		The aperture shape type can be accessed for reading by its index with the property obj.ApShapeId (set to 1 (polygon type) for a default RenderModel object).
%		The method obj.computeAperture must be called before rendering the image to take the new value into account.
%
% - obj.setApShapeId(ApShapeId)
%		Set the index of the aperture shape type.
%		See below for a description of aperture types and corresponding indices.
%		The aperture shape index can be accessed for reading with the property obj.ApShapeId (set to 1 (polygon type) for a default RenderModel object).
%		The method obj.computeAperture must be called before rendering the image to take the new value into account.
%
% - obj.setApThickness(apThickness)
%		Set the aperture thickness: value form 0 (only the border of the shape) to 1 (full shape).
%		The aperture thickness can be accessed for reading with the property obj.apThickness (set to 1 for a default RenderModel object).
%		The method obj.computeAperture must be called before rendering the image to take the new value into account.
%
% - obj.setApAngle(apAngle)
%		Set the rotation angle of the aperture shape (in radians).
%		The rotation angle can be accessed for reading with the property obj.apAngle (set to 0 for a default RenderModel object).
%		The method obj.computeAperture must be called before rendering the image to take the new value into account.
%
% - obj.setNumBlades(numBlades)
%		Set the number of blades parameter (number of sides of the polygon, in the case of a polygon aperture shape type). The minimum value is 3.
%		The rotation angle can be accessed for reading with the property obj.numBlades (set to 5 for a default RenderModel object).
%		The method obj.computeAperture must be called before rendering the image to take the new value into account.
%
%--------------------------------------------------------------------------
%
% Aperture types:
%	- 'polygon' (ApShapeId=1) : Regular polygon aperture.
%		The radius parameter coresponds to the radius of the circumscribed circle of the polygon.
%		Note 1: The 'polygon' type has similar effect to the 'disk' or 'ring' types when the number of blades is large.
%		Note 2: When numblades=4, the effect of the 'polygon' type differs from the 'rect' because of the definition of the radius. 
%
%	- 'disk' (ApShapeId=2): Perfect disk.
%		The aperture radius parameter coresponds to the radius of the disk.
%		Angle, thickness and number of blades parameters have no effect for this aperture type.
%
%	- 'ring' (ApShapeId=3): Ring aperture.
%		Similar to the disk aperture but the aperture thickness parameter can be used to control the thickness of the ring (small thickness values increase the radius of the inner circle).
%		The 'ring' type approximates the 'disk' type when the thickness parameter is 1.
%
%	- 'rect' (ApShapeId=4): Square aperture.
%		The radius parameter corresponds to the half side length of the square.
%		The thickness and number of blades parameters have no effect for this aperture type.
%
%--------------------------------------------------------------------------
%
% Usage example:
% rMod = RenderModel(FDL, fullSize, crop, Disps, 0);% Initialise the RenderModel object with a FDL (disable disparity estimation using a dummy value 0 for the disparity map argument).
% rMod.setApShape('disk');	%Configure disk aperture shape.
% rMod.setRadius(5);		%Configure aperture radius parameter.
% rMod.setFocus(2);			%Configure focus parameter.
% rMod.computeAperture();	%Computes the aperture shape (needed because of the call to setApShape).
% rMod.renderImage();  		%Render the image.
% figure,imshow(rMod.Image);%Display the final image (in uint8 format).
%--------------------------------------------------------------------------
%
% See also RenderAppMain

classdef RenderModel < Observable
%Note: inherits from the custom class Observable instead of the built-in
%handle class to overcome a bug in the event notification (in some cases
%the notification does not produce any effect).
    
    properties ( Access = private )
        FDL
        wx
        wy
        Disps
        
        usingGPU
        
        Image_
        crop
        
        Apfft, dWu, dWv, uC, vC
        
        xC, yC
        even_fft
        
        SpatialApRes = 50;       %spatial resolution of the aperture radius.
        FreqApResIncrease = 400; %zero-padding of the aperture (for a higher Fourier domain resolution).
        
        radCorrection
    end
    
    properties ( Access = public, Constant )
        ApShapes={'polygon','disk','ring','rect'}
        numApShapes = length(RenderModel.ApShapes);
    end
    
    properties ( SetAccess = private)
        ApShapeId = 1
        numBlades = 5
        apThickness = 1
        apAngle = 0;

        s = single(0)
        radius = single(0)
        trueRadius = 0
        u0 = 0
        v0 = 0
        
        skipInvTr = false
        
        Ap
        Image=0
        DispMap=[]
    end
    
    properties
        isLinear = false
        gammaOffset = 0
    end
%{
    events
      ChangeRadius
      ChangeFocus
      ChangePosition
      ToggleSkipInvTr
      ChangeApShape
      ChangeNumBlades
      ChangeApThickness
      ChangeApAngle
   end
%}
    
    methods
        function obj = RenderModel(FDL, fullSize, crop, Disps, DispMap, isLinear, gammaOffset, useGPU)
            
            if(nargin==0)
                return;
            end
            
            %Check if GPU should be used
            hasGPU = parallel.gpu.GPUDevice.isAvailable;
            if(~exist('useGPU','var') ||  isempty(useGPU) )
                useGPU = hasGPU;
            elseif(useGPU && ~hasGPU)
                warning('Impossible to use GPU. Only CPU will be used.');
            end
            obj.usingGPU = useGPU & hasGPU;
            
            
            if(exist('isLinear','var') && ~isempty(isLinear)), obj.isLinear=isLinear;end
            if(exist('gammaOffset','var') &&  ~isempty(gammaOffset)), obj.gammaOffset=gammaOffset;end
            
            tic;fprintf('Initializating data...');
            
            %Dimension variables
            [wx,wy] = meshgrid(1:fullSize(2),1:fullSize(1));
            obj.xC = ceil((fullSize(2)+1)/2);
            obj.yC = ceil((fullSize(1)+1)/2);
            wx=single((wx(:,1:obj.xC)-obj.xC)/(fullSize(2)));
            wy=single((obj.yC-wy(:,1:obj.xC))/(fullSize(1)));

            nChan = size(FDL,3);
            numDisp = length(Disps);
            obj.even_fft = 1-mod(fullSize(1:2),2);
            if(isa(crop,'BorderParams'))
                obj.crop = [crop.L,crop.R,crop.T,crop.B];
            else
                obj.crop = crop;
            end
            
            %Prepare data for single precision format and GPU if needed.
            if(obj.usingGPU)
                obj.FDL = gpuArray(single(FDL(:,1:obj.xC,:,:)));
                obj.wx = gpuArray(single(wx));
                obj.wy = gpuArray(single(wy));
                obj.Disps = zeros(1,1,1,numDisp,'single','gpuArray');
                obj.Image_ = zeros([fullSize(1:2), nChan],'single','gpuArray');
            else
                obj.FDL = single(FDL(:,1:obj.xC,:,:));
                obj.wx = single(wx);
                obj.wy = single(wy);
                obj.Disps = zeros(1,1,1,numDisp,'single');
                obj.Image_ = zeros([fullSize(1:2), nChan],'single');
            end
            obj.Disps(1,1,1,:)=Disps;
            
            
            
            t=toc;fprintf(['\b\b\b (done in ' num2str(t) 's)\n']);
            
            
            %Disparity map
            if(~exist('DispMap','var') || isempty(DispMap))
                tic;fprintf('Estimating disparity map...');
                %{
                %Estimate disparity map (quite inaccurate, but it is already useful for refocus on click)
                DL=zeros(size(FDL));for ch=1:nChan,for k=1:numDisp,DL(:,:,ch,k) = (ifft2(ifftshift(FDL(:,:,ch,k))));end,end;DL=real(DL);
                %[~,obj.DispMap]=max(mean(DL,3),[],4);
                %obj.DispMap=Disps(obj.DispMap(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2)));
                DL = squeeze(mean(DL(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2),:,:),3));
                obj.DispMap = sum(bsxfun(@times,DL.^2,permute(Disps(:),[3 2 1])),3)./sum(DL.^2,3);
                %}
                
                
                %%{
                %Disparity estimation based on a defocus cue (select for each pixel the refocus disparity minimizing the differences of the images obtained with varying aperture radius).
                ApShapeId0 = obj.ApShapeId;
                radius0 = obj.radius;
                s0 = obj.s;
                obj.ApShapeId=4;
                
                obj.computeAperture();
                if(obj.usingGPU)
                    MinErrMap = gpuArray(inf(size(obj.Image_,1),size(obj.Image_,2),2,'single'));
                else
                    MinErrMap = inf(size(obj.Image_,1),size(obj.Image_,2),2,'single');
                end
                obj.DispMap = zeros(size(obj.Image_,1),size(obj.Image_,2));
                for k=1:numDisp
                    obj.s = Disps(k);
                    obj.radius=0;
                    obj.renderImage();
                    SAI = mean(real(obj.Image_),3);
                    rList = .5:.5:3;
                    MinErrMap(:,:,2)=0;
                    for r = rList
                        obj.radius=r;
                        obj.renderImage();
                        MinErrMap(:,:,2) = MinErrMap(:,:,2) + (SAI-mean(real(obj.Image_),3)).^2;
                    end
                    [MinErrMap(:,:,1),Id] = min(MinErrMap,[],3);
                    obj.DispMap(Id==2) = Disps(k);
                end
                obj.DispMap = obj.DispMap(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2));
                obj.radius = radius0;
                obj.s = s0;
                obj.ApShapeId = ApShapeId0;
                %}
                t=toc;fprintf(['\b\b\b (done in ' num2str(t) 's)\n']);
                clear DL
            else
                obj.DispMap = DispMap;
            end
            
            %Initialize Aperture data
            obj.computeAperture();
            
            %Intialize Image
            obj.renderImage();
            
        end
        
        
        %Image render function
        function renderImage(obj)
            %Render first half of the spectrum.
            if(obj.usingGPU)
                obj.Image_(:,1:obj.xC,:) = RenderHalfFT(obj.FDL, obj.wx, obj.wy, obj.Disps, obj.Apfft, obj.dWu, obj.dWv, obj.uC, obj.vC, obj.s, obj.radius, obj.u0, obj.v0);
            else
                obj.Image_(:,1:obj.xC,:) = RenderHalfFT_cpu(obj.FDL, obj.wx, obj.wy, obj.Disps, obj.Apfft, obj.dWu, obj.dWv, obj.uC, obj.vC, obj.s, obj.radius, obj.u0, obj.v0);
            end
            
            if(obj.skipInvTr)
                %Reconstruct second half of the spectrum using symmetries.
                obj.Image_(obj.yC+1:end,obj.xC+1:end,:) = conj(obj.Image_(obj.yC-1:-1:1+obj.even_fft(1),obj.xC-1:-1:1+obj.even_fft(2),:,:));
                obj.Image_(1+obj.even_fft(1):obj.yC, obj.xC+1:end,:) = conj(obj.Image_(end:-1:obj.yC, obj.xC-1:-1:1+obj.even_fft(2),:,:));
                %render the fourier magnitude spectrum (only for display, the true fourier transform remains in the internal image representation (i.e. obj.Image_).
                obj.Image = gather( real(uint8(abs(obj.Image_(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2),:)*255*10000/size(obj.Image_,1)/size(obj.Image_,2)))));
            else
                %Apply inverse Fourier Transform directly from the half spectrum.
                obj.Image_ = ifft(conj(fliplr(ifft(ifftshift(obj.Image_(:,1:obj.xC,:),1),[],1))),size(obj.Image_,2),2,'symmetric');
                
                if(obj.isLinear)
                %Apply gamma correction if needed.
                    obj.Image_(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2),:) = RenderModel.BT709_gamma(obj.Image_(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2),:))-obj.gammaOffset;
                end
                %Format conversion.
                obj.Image = gather( uint8(255*obj.Image_(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2),:)));
            end
        end
        
        %Aperture update function
        function computeAperture(obj)
            if(obj.usingGPU)
                [obj.Ap, Apfft, obj.dWu, obj.dWv, obj.uC, obj.vC, obj.radCorrection] = buildAperture(RenderModel.ApShapes{obj.ApShapeId},obj.SpatialApRes,obj.FreqApResIncrease,[obj.apThickness obj.numBlades obj.apAngle]);
                obj.trueRadius = obj.radius * obj.radCorrection / obj.SpatialApRes;
                Apfft = Apfft/Apfft(obj.vC,obj.uC);
                obj.Apfft = gpuArray(complex(single(Apfft)));
            else
                [obj.Ap, obj.Apfft, obj.dWu, obj.dWv, obj.uC, obj.vC, obj.radCorrection] = buildAperture(RenderModel.ApShapes{obj.ApShapeId},obj.SpatialApRes,obj.FreqApResIncrease,[obj.apThickness obj.numBlades obj.apAngle]);
                obj.trueRadius = obj.radius * obj.radCorrection / obj.SpatialApRes;
                obj.Apfft = obj.Apfft/obj.Apfft(obj.vC,obj.uC);
            end
        end
        
        %Get internal rendered image representation either with or without
        %cropping the padded borders (no crop by default / if the inverse Fourier transform is skipped, the crop is never applied).
        function I = getInternalImage(obj,applyCrop)
            if(~exist('applyCrop','var')), applyCrop=false;end
            if(applyCrop && ~obj.skipInvTr)
                I = gather(obj.Image_(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2),:));
            else
                I = gather(obj.Image_);
            end
        end
        
        %Set methods
        function setRadius(obj,radius)
            if(radius ~= obj.radius)
                obj.radius = max(0,radius);
                obj.trueRadius = obj.radius * obj.radCorrection / obj.SpatialApRes;
                notify(obj,'ChangeRadius');
            end
        end
        
        function setFocus(obj,s)
            if(s ~= obj.s)
                obj.s = s;
                notify(obj,'ChangeFocus');
            end
        end
        
        function setPosition(obj,u0,v0)
            if(u0 ~= obj.u0 || v0 ~= obj.v0)
                obj.u0 = u0;
                obj.v0 = v0;
                notify(obj,'ChangePosition');
            end
        end
        
        function setSkipInvTr(obj, skipInvTr)
            if(obj.skipInvTr ~= skipInvTr && (obj.skipInvTr==false || skipInvTr==false)) %avoid recomputing if skipInvTr remains true or remains false.
                obj.skipInvTr = skipInvTr;
                notify(obj,'ToggleSkipInvTr');
            end
        end
        
        function setApShape(obj, ApShapeName)
            id = find(strcmp(ApShapeName,RenderModel.ApShapes));
            if(~isempty(id))
                if(obj.ApShapeId ~= id)
                    obj.ApShapeId = id;
                    notify(obj,'ChangeApShape');
                end
            else
                warning(['unknown aperture shape (', ApShapeName, ')'] );
            end
            
        end
        
        function setApShapeId(obj, ApShapeId)
            if(ApShapeId>0 && ApShapeId <= RenderModel.numApShapes && round(ApShapeId)==ApShapeId && obj.ApShapeId ~= ApShapeId)
                obj.ApShapeId = ApShapeId;
                notify(obj,'ChangeApShape');
            end
            
        end
        
        function setNumBlades(obj, numBlades)
            if(numBlades ~= obj.numBlades)
                obj.numBlades = round(max(numBlades,3));
                notify(obj,'ChangeNumBlades');
            end
        end
        
        function setApThickness(obj, apThickness)
            %if(apThickness >= 1 && apThickness <= obj.SpatialApRes && apThickness ~= obj.apThickness)
            if(apThickness ~= obj.apThickness)
                %obj.apThickness = round(min(max(apThickness,1),1+obj.SpatialApRes));
                obj.apThickness = min(max(apThickness,0),1);
                notify(obj,'ChangeApThickness');
            end
        end
        
        function setApAngle(obj, apAngle)
            if(apAngle ~= obj.apAngle)
                obj.apAngle = apAngle;
                notify(obj,'ChangeApAngle');
            end
        end
        

        function saveFDL(obj,filename, U,V)
            FDL = gather(obj.FDL);
            Disps = squeeze(gather(obj.Disps))';
            crop = obj.crop;
            fullSize = [size(obj.Image_,1), size(obj.Image_,2)];
            VarList = {'FDL','Disps','crop','fullSize'};
            if(obj.isLinear)
                isLinear = obj.isLinear;
                gammaOffset = obj.gammaOffset;
                VarList(end+1:end+2)={'isLinear','gammaOffset'};
            end
            if(exist('U','var') && exist('V','var'))
                VarList(end+1:end+2)={'U','V'};
            end
            save(filename, VarList{:});
        end
        
    end %methods
    
    
    methods ( Static )
        %Apply BT709 standard gamma correction to convert from linear RGB data to BT709.
        %(negative input values are clipped to 0).
        function Iout = BT709_gamma(Iin)
            Mask = Iin < 0.018;
            Iout = max(0,Iin * 4.5) .* Mask + (1.099*max(Iin,0).^0.45 - 0.099) .* (~Mask);
        end
        
    end %methods ( Static )
    
    
end
