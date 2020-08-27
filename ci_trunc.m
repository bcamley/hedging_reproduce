function cis = ci_trunc(c)
% This function implements the chemotactic index function
% in terms of one way to write this generalized Laguerre polynomial
Lhalf = @(x) exp(x/2).*( (1-x).*besseli(0,-x/2)-x.*besseli(1,-x/2));
cis = sqrt(2/pi)*c./Lhalf(-0.5*c.^2);
cis(c>50) = 1-0.5./c(c>50).^2; % put in asymptotic value of function here
cis(c<-50) = -1+0.5./c(c<-50).^2; 

cis(isnan(cis)&(c>1)) = 1;
cis(isnan(cis)&(c<-1)) = -1;