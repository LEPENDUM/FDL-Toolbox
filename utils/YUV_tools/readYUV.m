% Reads a YUV sequence file into separate YUV files of each frame.
% Inputs :
%    - seq_filename  : filename (including path) of the sequence YUV file.
%    - resX, resY    : horizontal and vertical resolution.
%    - chromaFormat  : (string) '444', '420', '422', or '400'.
%    - bitdepth      : Bitdepth (supports any positive value up to 16):
%                      - If bitdepth <= 8, the return type is uint8.
%                      - If bitdepth > 8, the return type is uint16.
%    - numFrames     : (Optional) Number of frames to read (default=1).
%    - frameSkip     : (Optional) Number of frame to skip at the begining of the sequence file (default=0).
%    - UpSampFunc    : (Optional) function handler for the Upsampling. By default, imresize is used in 'bicubic' mode.
%                      The function must take as inputs:
%                       1. I: Images or sequence of images with U and V components, i.e. array of size (resY,resX,2,number of frames).
%                       2. sizeY : target vertical resolution.
%                       3. sizeX : target horizontal resolution.
%                      ex : UpSampFunc = @(I,sizeY,sizeX) imresize(I,[sizeY sizeX],'lanczos2');
%
% Output YUV :
% Single matrix containing all frames :
% For YUV400 format : YUV(resY,resX,numFrames).
% Otherwise         : YUV(resY,resX,3,numFrames).

function [YUV] = readYUV( seq_filename, resX, resY, chromaFormat, bitdepth, numFrames, frameSkip, UpSampFunc )

if(~exist('UpSampFunc','var'))
    UpSampFunc = @(I,sizeY,sizeX) imresize(I,[sizeY sizeX],'bicubic');
end
if(~exist('numFrames','var'))
    numFrames=1;
end
if(~exist('frameSkip','var'))
    frameSkip=0;
end

bitdepth=8*ceil(bitdepth/8);

if(bitdepth==8)
    dataType='uint8';
elseif (bitdepth==16)
    dataType='uint16';
else
    error('Data bitdepth should be between 1 and 16.');
end

switch(chromaFormat)
    case '444'
        resXC=resX; resYC=resY;
        img_filesize = (bitdepth+2*bitdepth)/8 * resX * resY;
    case '400'
        resXC=0; resYC=0;
        img_filesize = bitdepth/8 * resX * resY;
    case '420'
        resXC=round(resX/2); resYC=round(resY/2);
        img_filesize = bitdepth/8 * resX * resY + bitdepth/4 * resXC * resYC;
    case '422'
        resXC=round(resX/2); resYC=resY;
        img_filesize = bitdepth/8 * resX * resY + bitdepth/4 * resXC * resYC;
    otherwise
        error('Unknown chroma format');
end

if(strcmp(chromaFormat,'400'))
    YUV = zeros(resY,resX,numFrames,dataType);
else
    YUV = zeros(resY,resX,3,numFrames,dataType);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Read from file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen( seq_filename );
try
    fseek(fid, frameSkip*img_filesize, 'bof');

    for i = 0 : numFrames-1
        Y = fread(fid, [resX,resY], ['*' dataType])';
        UV = cat(3,fread(fid, [resXC,resYC], ['*' dataType])',fread(fid, [resXC,resYC], ['*' dataType])');

        switch(chromaFormat)
        case '444'
            YUV(:,:,:,i+1) = cat(3,Y,UV);
        case '400'
            YUV(:,:,i+1) = Y;
        case '420'
            UVresiz = UpSampFunc(UV, resY, resX);
            YUV(:,:,:,i+1) = cat(3,Y,UVresiz);
        case '422'
            UVresiz = UpSampFunc(UV, resY, resX);
            YUV(:,:,:,i+1) = cat(3,Y,UVresiz);
        end
    end
    
%close file even if an error is raised during execution.
    fclose(fid);
catch exception
    fclose(fid);
    rethrow(exception);
end
