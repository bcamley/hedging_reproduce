function h = contourf_better(X,Y,C,contours,factor)
if(nargin<5)
    factor =5;
end
C = interp2(C,factor,'linear');
X = interp2(X,factor,'linear');
Y = interp2(Y,factor,'linear');

h = contourf(X,Y,C,contours);
shading interp