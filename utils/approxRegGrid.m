% Fits scaling and offset parameters between a set of 2D points and a theoretical grid.
%Inputs:
% - U: Horizontal coordinates of input set of 2D points.
% - V: Vertical coordinates of input set of 2D points.
% - GridU: Horizontal coordinates of set of 2D points to fit.
% - GridV: Vertical coordinates of set of 2D points to fit.
%          U, V, GridU and GridV must have the same number of elements.
%
% Output
% - delatU: Horizontal scale to apply to the (GridU,GridV) points to approximate the (U,V) points.
% - deltaV: Vertical scale to apply to the (GridU,GridV) points to approximate the (U,V) points.
% - u0: Horizontal offset to apply to the (GridU,GridV) points to approximate the (U,V) points.
% - v0: Vertical offset to apply to the (GridU,GridV) points to approximate the (U,V) points.
%
% The model (delatU*GridU + u0, deltaV*GridV + v0) is determined to approximate (U,V) in the least square sense.

function [delatU, deltaV, u0, v0] = approxRegGrid(U,V,GridU,GridV)

UV = [U(:); V(:)];
numPts = numel(U);
M = zeros(2*numPts,4);
M(1:numPts,1) = GridU(:);
M(1:numPts,3) = 1;
M(numPts+1:end,2) = GridV(:);
M(numPts+1:end,4) = 1;

x=M\UV;
delatU=x(1);
deltaV=x(2);
u0=x(3);
v0=x(4);
end