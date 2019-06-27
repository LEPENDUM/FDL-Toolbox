%Main class for light field rendering from the FDL model.
%Usage example:
% rMod = RenderModel(FDL, fullSize, crop, Disps, 0);
% rMod.setRadius(r);
% rMod.setFocus(f);
% rMod.renderImage();
% figure,imshow(rMod.Image);

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
        
        %Constructor with parameters:
        % -FDL   : Fourier Disparity layers (only left half of the spectrum (i.e. negative or null horizontal frequencies) may be given).
        % -fullSize : full image size in the format [verticalResolution, horizontalResolution], including the padded borders.
        % -crop  : number of pixels to crop (to remove the padded borders) on each side of the image : [top, bottom, left, right].
        % -Disps : list of disparity values associated with each layer.
        % -DispMap: (Optional) Disparity map used for refocusing automatically when the user clicks on the image (if not given, a disparity map will be estimated).
        % -isLinear: (Optional) Set to true if gammaCorrection is needed.
        % -gammaOffset: (Optional) offset parameter to apply after gamma correction (only used if isLinear is true).
        % -useGPU : (Optional) true(default): use GPU acceleration if a CUDA Device is available and matlab's parallel processing toolbox is active. / false: CPU only.
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
            wx=single((wx(:,1:obj.xC)-obj.xC)/(fullSize(2)-1));
            wy=single((obj.yC-wy(:,1:obj.xC))/(fullSize(1)-1));

            nChan = size(FDL,3);
            numDisp = length(Disps);
            obj.even_fft = 1-mod(fullSize(1:2),2);
            obj.crop = crop;

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
            %Render first half of the spectrum
            if(obj.usingGPU)
                obj.Image_(:,1:obj.xC,:) = RenderHalfFT(obj.FDL, obj.wx, obj.wy, obj.Disps, obj.Apfft, obj.dWu, obj.dWv, obj.uC, obj.vC, obj.s, obj.radius, obj.u0, obj.v0);
            else
                obj.Image_(:,1:obj.xC,:) = RenderHalfFT_cpu(obj.FDL, obj.wx, obj.wy, obj.Disps, obj.Apfft, obj.dWu, obj.dWv, obj.uC, obj.vC, obj.s, obj.radius, obj.u0, obj.v0);
            end
            %Reconstruct second half of the spectrum using symmetries
            obj.Image_(obj.yC+1:end,obj.xC+1:end,:) = conj(obj.Image_(obj.yC-1:-1:1+obj.even_fft(1),obj.xC-1:-1:1+obj.even_fft(2),:,:));
            obj.Image_(1+obj.even_fft(1):obj.yC, obj.xC+1:end,:) = conj(obj.Image_(end:-1:obj.yC, obj.xC-1:-1:1+obj.even_fft(2),:,:));
             
            if(obj.skipInvTr)
                obj.Image_ = abs(obj.Image_)/50;
                %obj.Image_ = cat(3,real(obj.Image_(:,:,2)),imag(obj.Image_(:,:,2)),abs(obj.Image_(:,:,2)))/50;
            else
            %Apply inverse Fourier Transform
                obj.Image_ = ifft2(ifftshift(ifftshift(obj.Image_,1),2));
            end

%%%%%%%%%%%%%
            if(obj.isLinear)
                obj.Image_(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2),:) = RenderModel.BT709_gamma(obj.Image_(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2),:))-obj.gammaOffset;
            end
            obj.Image = gather( real(uint8(255*obj.Image_(1+obj.crop(3):end-obj.crop(4),1+obj.crop(1):end-obj.crop(2),:))));
            
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
