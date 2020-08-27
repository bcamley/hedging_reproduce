function h = pcolor_better(X,Y,C,factor)
if(nargin<4)
    factor =5;
end
C = interp2(C,factor,'linear');
X = interp2(X,factor,'linear');
Y = interp2(Y,factor,'linear');

h = pcolor(X,Y,C);
shading interp