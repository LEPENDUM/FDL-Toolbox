% ElasticNet regression with ADMM method:
% X = argmin_{X} 0.5*||A*X-B||_{2}^{2} + .5*lambda*||X||_{2}^{2} + gamma*||X||_{1}
%
function [X,cost,stop,Z] = elasticNetADMM(B, A, gamma, lambda)

% Get dimensions
c = size(B,2);
r = size(A,2);

% Precompute matrices:
[Q,Sig2] = eig(A'*A);
ATB=A'*B;
Sig2 = diag(Sig2);

% Parameter initialization:
L = zeros(r,c); % Initialize Lagragian
rhoMult=1.01;   % Parameter controling convergence speed / accuracy tradeoff.
maxIter = 500;  % Maximum number of iterations
Z = zeros(r,c); % Initialize Z
if(lambda<1e-4) % Initialize rho
    rho=1e-4; 
else
    rho=0;
end

% Soft thresholding function
fast_sthresh = @(x,th) sign(x).*max(abs(x) - th,0);

% Norm functions
norm2 = @(x) x(:)'*x(:);
norm1 = @(x) sum(abs(x(:))); 

cost = [];
stop = [];

% First iteration out of the loop to initialize rho
X = Q*(diag(1./(Sig2+rho+lambda)))*Q'*(ATB + rho*Z - L);
rho = 1.2*gamma/max(abs(X(:)));
Z = fast_sthresh(X + L/rho, gamma/rho);
L = L + rho*(X - Z);
rho = rho*rhoMult;
cost(1) = 0.5 * ( norm2(B - A*X) + lambda*norm2(X) ) + gamma*norm1(X);
for n = 2:maxIter
%for n = 1:maxIter

    % Solve sub-problem to solve X
    X = Q*(diag(1./(Sig2+rho+lambda)))*Q'*(ATB + rho*Z - L);

    % Solve sub-problem to solve Z
    Z = fast_sthresh(X + L/rho, gamma/rho);
    
    % Update the Lagrangian
    L = L + rho*(X - Z);

    % rho update to accelerate convergence
    rho = rho*rhoMult;
    
    % get the current cost
    cost(n) = 0.5 * ( norm2(B - A*X) + lambda*norm2(X) ) + gamma*norm1(X);
    stop(n) = abs(cost(n)-cost(n-1))/cost(n-1);
    if(stop(n)<1e-3) break;end
end