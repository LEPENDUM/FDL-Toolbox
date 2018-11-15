%Rescale the view position and disparity parameters assuming the views
%are located approximtely on an integer grid of size nU x nV.
%The function also centers the average view position to the position (0,0).
function [U,V,D,sigU,sigV,sigD] = ScaleParams(U,V,D,nU,nV,sigU,sigV,sigD)
U=U-mean(U);
V=V-mean(V);
maxU = max(abs(U));
maxV = max(abs(V));

scale = (floor(nU/2)/maxU+floor(nV/2)/maxV)/2;
U=U*scale;
V=V*scale;
D = D/scale;

if(exist('sigU','var')), sigU = sigU*scale; elseif(nargout>3), sigU=zeros(size(U));end
if(exist('sigV','var')), sigV = sigV*scale; elseif(nargout>4), sigV=zeros(size(V));end
if(exist('sigD','var')), sigD = sigD/scale; elseif(nargout>5), sigD=zeros(size(D));end