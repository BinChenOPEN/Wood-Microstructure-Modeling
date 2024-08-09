[x,y] = ndgrid(1:1000,1:1000);
pointCoord = [500,500];

% k          = [0.02,0.015,5,3];
% k          = [0.1,0.08,2,15];
k = [0.04,0.03,2,6];

% k      = [0.04,0.035,7,9];
u_temp     = localDistort(x,y,pointCoord(1),pointCoord(2),k);

surf(u_temp),axis equal, axis tight, shading interp,axis on, colorbar,view([0,90])
