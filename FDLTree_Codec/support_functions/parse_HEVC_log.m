%Retrieve information from HM decoder log files:
% nbBits_HEVC: Number of bits per frame.
% PSNR_HEVC: PSNR of each frame.
% PSNR_HEVC: Encoding time for each frame.
function [nbBits_HEVC, PSNR_HEVC, ET_HEVC] = parse_HEVC_log(filename,nbFrames)

fid = fopen(filename);

if(nargin<2)
    nbFrames=0;
end
tmpBits=zeros(1,nbFrames);  % Number of bits of the frame
tmpPSNRY=zeros(1,nbFrames); % Y-PSNR of the frame
tmpPSNRU=zeros(1,nbFrames); % U-PSNR of the frame
tmpPSNRV=zeros(1,nbFrames); % V-PSNR of the frame
tmpET=zeros(1,nbFrames);    % Encoding time of the frame
poc=zeros(1,nbFrames); % picture order count (not returned : only used to reorder frames in the display order)

numFrame = 0;
while (~feof(fid) )
    line=fgets(fid);
    C = textscan(line, '%s %d %s %s %s %s %s %s %s %s %s %d %s %s %f %s %s %f %s %s %f %s %s %d',1);
    if(~isempty(C{1,1}) && strcmp(C{1,1},'POC'))
        numFrame = numFrame+1;
        poc(numFrame)      = C{1,2};
        tmpBits(numFrame)  = C{1,12};
        tmpPSNRY(numFrame) = C{1,15};
        tmpPSNRU(numFrame) = C{1,18};
        tmpPSNRV(numFrame) = C{1,21};
        tmpET(numFrame) = C{1,24};
    end
end

fclose(fid);

%Reorder Frames using poc
nbFramesHEVC = numFrame;
nbBits_HEVC=zeros(1,nbFramesHEVC);
PSNR_HEVC=zeros(3,nbFramesHEVC);
ET_HEVC=zeros(1,nbFramesHEVC);

for i=1:nbFramesHEVC
    nbBits_HEVC(1+poc(i)) = tmpBits(i);
    PSNR_HEVC(1,1+poc(i)) = tmpPSNRY(i);
    PSNR_HEVC(2,1+poc(i)) = tmpPSNRU(i);
    PSNR_HEVC(3,1+poc(i)) = tmpPSNRV(i);
    ET_HEVC(1+poc(i))     = tmpET(i);
end