function SNR = SNR_multi_snapshot(c,logKvec,q,fracs)

Ntypes = length(logKvec);

if(nargin<4)
    fracs = (1/Ntypes)*ones(size(logKvec));
else
    fracs = abs(fracs) / sum(abs(fracs));  % to ensure normalization
end

thr = 1e8; % threshold for asymptotic form

SNR = zeros(size(c));


for j = 1:length(logKvec)
    KD = exp(logKvec(j));
    beta = c*KD./(KD+c).^2;
    beta(c>KD*thr) = KD./c(c>KD*thr);
    SNR = SNR + fracs(j)*beta;   % assuming equal fraction of each type
end
SNR = ((q.g^2)*q.nr/16)*SNR;


end