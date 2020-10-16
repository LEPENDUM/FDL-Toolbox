%Retrieve Decoding time from HM decoder log files.
function [DT_HEVC] = parse_HEVC_dec_log(filename,nbFrames)

fid = fopen(filename);

if(nargin<2)
    nbFrames=0;
end
tmpDT=zeros(1,nbFrames); % Decoding time of the frame
poc=zeros(1,nbFrames);   % picture order count (not returned : only used to reorder frames in the display order)

numFrame = 0;
while (~feof(fid) )
    line=fgets(fid);
    C = textscan(line, '%s %d %s %s %s %s QP %s %s %s %f',1);
    if(~isempty(C{1,1}) && strcmp(C{1,1},'POC'))
        numFrame = numFrame+1;
        poc(numFrame)      = C{1,2};
        tmpDT(numFrame) = C{1,10};
    end
end

fclose(fid);

%Reorder Frames using poc
nbFramesHEVC = numFrame;
DT_HEVC=zeros(1,nbFramesHEVC);

for i=1:nbFramesHEVC
    DT_HEVC(1+poc(i))     = tmpDT(i);
end