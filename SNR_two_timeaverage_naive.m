function SNR = SNR_two_timeaverage_naive(c,fa,q)

thr = 1e8; % threshold for asymptotic form
betaA = c*q.KA./(q.KA+c).^2;
% put in asymptotic form:
betaA(c>q.KA*thr) = q.KA./c(c>q.KA*thr);
betaB = c*q.KB./(q.KB+c).^2;
betaB(c>q.KB*thr) = q.KB./c(c>q.KB*thr);

tauA = 1./(1+c/q.KA);
tauB = (1/(q.kminBAratio))*(1./(1+c/q.KB));   % both scaled by k_-^A



%SNR = ( fa*(c./(c+KA)) + (1-fa)*(c./(c+KB))*kminBAratio);
SNR = (q.kminusT)*((q.g^2)*q.nr/32)*((fa*betaA+(1-fa)*betaB).^2)./(fa*betaA.*tauA + (1-fa)*betaB.*tauB);


