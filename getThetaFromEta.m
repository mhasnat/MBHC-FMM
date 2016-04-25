function x = getThetaFromEta(D)
% Conversion from expectation parameter (eta) to natural parameter (theta)
% using the Newton–Raphson method.
% See Section 3.4 (Eq. 12) of [1] for details

% INPUT:
% D: value of eta (expectation parameter)

% OUTPUT
% x: value of theta (natural parameter)

% Reference:
% [1] Hasnat et al., Model-based hierarchical clustering with Bregman 
% divergences and Fishers mixture model: application to depth image analysis. 
% Statistics and Computing, 1-20, 2015.

% Author: Md Abul HASNAT

x = 1; % initial value

diff = 100;
x_v(1) = 500000;
ii=2;

while(diff>0.001)
    a = tanh(x).^(-1);
    b = x.^(-1);
    
    fx = a - b - D;
    f1x = 1 - a^2 + b^2;
    
    x = x - fx/f1x;
    
    x_v(ii) = x;
    
    diff = abs(x_v(ii-1) - x);
    
    ii = ii + 1;
end

end