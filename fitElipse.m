function [P] = fitElipse(pointCoord)

A  = [pointCoord.^2,pointCoord];
b  = ones(length(A),1);

if size(A,1) == 4
    x  = A\b;
else
    x  = (A'*A)\A'*b;
end

x  = sign(x(1))*x;


h  = -x(3)/x(1)/2;
k  = -x(4)/x(2)/2;

r1 = sqrt(2*(x(1)*h^2+x(2)*k^2-1)/(x(1)+x(2)+abs(x(1)-x(2))));
r2 = sqrt(2*(x(1)*h^2+x(2)*k^2-1)/(x(1)+x(2)-abs(x(1)-x(2))));

P = [r1;r2;h;k];