function [cf,line1, line2]=wplr(x,y,wt)
%
% Weighted piecewise linear regression fit (WPLR − τ ) method for model
% selection (i.e., selecting number of clusters automatically). See Section
% 5.3 of Ref. [1]

% Reference:
% [1] Hasnat et al. Model-based hierarchical clustering with Bregman 
% divergences and Fishers mixture model: application to depth image analysis. 
% Statistics and Computing, 1-20, 2015.
% 

% Author: Md Abul HASNAT

deg = [1 1];

ss=Inf(1,size(x,2));

for c=2:(size(x,2)-1)
    x1=x(1:c);
    y1=y(1:c);
        
    linearCoef1 = polyfit(x1,y1,deg(1));
    linearFit1{c} = polyval(linearCoef1,x1);
    err1 = sum((linearFit1{c} - y1).^2);
    
    x2=x(c:end);
    y2=y(c:end);
    
    linearCoef2 = polyfit(x2,y2,deg(2));
    linearFit2{c} = polyval(linearCoef2,x2);
    
    err2 = sum((linearFit2{c} - y2).^2);
    
    if(nargin<3)
        wt(1) = 1;
        wt(2) = 1;
    end
    
    ss(c)=sum((wt(1)*err1) + (wt(2)*err2));
end

[~,cf]=min(ss);
line1 = linearFit1{cf};
line2 = linearFit2{cf};
