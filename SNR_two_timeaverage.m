function SNR = SNR_two_timeaverage(c,fa,q)

caratio = c./(c+q.KA);
cbratio = c./(c+q.KB);

% use an asymptotic limit for c>>KA, c>>KB
thresh = 1e8;
caratio(c>thresh*q.KA) = 1-q.KA./c(c>thresh*q.KA);
cbratio(c>thresh*q.KB) = 1-q.KB./c(c>thresh*q.KB);

%SNR = ( fa*(c./(c+KA)) + (1-fa)*(c./(c+KB))*kminBAratio);
SNR = (q.kminusT)*((q.g^2)*q.nr/16)*( fa*caratio + (1-fa)*(cbratio)*q.kminBAratio);


