function SNR = SNR_two_snapshot(c,fa,q)

%SNR = NaN*ones(size(c));
%I = zeros(size(c));
thr = 1e8; % threshold for asymptotic form
betaA = c*q.KA./(q.KA+c).^2;
% put in asymptotic form:
betaA(c>q.KA*thr) = q.KA./c(c>q.KA*thr);
betaB = c*q.KB./(q.KB+c).^2;
betaB(c>q.KB*thr) = q.KB./c(c>q.KB*thr);

SNR = ((q.g^2)*q.nr/16)*(betaA*fa+(1-fa)*betaB);


end