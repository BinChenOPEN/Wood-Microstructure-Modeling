function u = localDistort(xgrid,ygrid,x0,y0,k)
% The four parameter in k are very important to control the distortion

sigma = k(3)/k(2);
k0  = k(4)/k(2);

u = exp(k(1)*(xgrid-x0))./((1+exp(k(1)*(xgrid-x0)))) ...
    - exp(k(2)*(xgrid-x0))./((1+exp(k(2)*(xgrid-x0))));
weight = k0*exp(-((ygrid-y0).^2)./(2*sigma^2));
u   = u.*weight;
