% Write YUV file with the specified chroma format ( 400 / 420 / 422 / 444 )
% For 420 and 422, by default, the downsampling is performed with the imresize function
% and the 'bicubic' option.
% Inputs :
% -YUV : Array of uint8 or uint16 containing yuv image or sequence (full chroma sampling).
%        The dimensions of the array should be (1.vertical dimension, 2.horizontal dimension, 3.colour channels, 4.frames).
%        Bitdepth for saving YUV file is based on input array type (8 for uint8, 16 for uint16).
% -filename : (string) file to write.
% -chromaFormat : (string) '444', '420', '422', or '400'.
% -append : (optional) true->append YUV image data to existing YUV file / false(default)->overwrite existing YUV file.
% -DnSampFunc : (optional) function handler for the Downsampling. By default, imresize is used in 'bicubic' mode.
%               The function must take as inputs:
%                1. I: Images or sequence of images with U and V components, i.e. array of size (resY,resX,2,number of frames).
%                2. sizeY : target vertical resolution.
%                3. sizeX : target horizontal resolution.
%               ex : DnSampFunc = @(I,sizeY,sizeX) imresize(I,[sizeY sizeX],'lanczos2');

function writeYUV( YUV, filename, chromaFormat, append, DnSampFunc )

[resY, resX, nChan, numFrames] = size(YUV);

datatype = class(YUV);
if(~isequal(datatype,'uint8') && ~isequal(datatype,'uint16'))
    error(['YUV image type ''' class(YUV) ''' is not supported. Use ''uint8'' or ''uint16'' data types.']);
end

if(nChan==1 && ~strcmp(chromaFormat,'400'))
    warning(['Chroma format ''' chromaFormat ''' cannot be used with a single channel image -> 400 used instead']);
    chromaFormat = '400';
elseif(nChan~=1 && nChan~=3)
    error('YUV input image must have 1 or 3 components')
end
    
if(~exist('append','var')), append=false;end

switch(chromaFormat)
    case '444'
            YUV = {YUV(:,:,1,:),YUV(:,:,2:3,:)};
    case '400'
            YUV = {YUV(:,:,1,:)};
    case '420'
        if(exist('DnSampFunc','var'))
            UVresiz = DnSampFunc(YUV(:,:,2:3,:), round(resY/2), round(resX/2));
        else
            UVresiz = imresize(YUV(:,:,2:3,:), .5, 'bicubic');
        end
        YUV = {YUV(:,:,1,:), UVresiz};
    case '422'
        if(exist('DnSampFunc','var'))
            UVresiz = DnSampFunc(YUV(:,:,2:3,:), resY , round(resX/2));
        else
            UVresiz = imresize(YUV(:,:,2:3,:), [resY round(resX/2)], 'bicubic');
        end
        YUV = {YUV(:,:,1,:), UVresiz};
    otherwise
        error('unknown chroma format');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Write to file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(append)
    f_permission='a';
else
    f_permission='w';
end

fid = fopen(filename,f_permission);
try
    for i=1:numFrames
        fwrite(fid,YUV{1}(:,:,i)',datatype);
        if(length(YUV)==2)
            fwrite(fid,permute(YUV{2}(:,:,:,i),[2,1,3]),datatype);
        end
    end
%close file even if an error is raised during execution.
    fclose(fid);
catch exception
    fclose(fid);
    rethrow(exception);
end
