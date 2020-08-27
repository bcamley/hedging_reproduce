function [fractions,maxCImean,deltaCI,percdeltaCI,SR,CS] = find_fraction_function(q)

opts = optimset('TolFun',1e-9,'TolX',1e-6);

switch(q.CIfunc)
    case 'ci'
        CIf = @(SNR) ci_trunc(sqrt(SNR));
    case 'hill'
        CIf = @(SNR) (SNR.^q.hilln)./(1+SNR.^q.hilln);
end

KA = q.KA;
KB = q.KB;

if(~isfield(q,'sigrels'))
    sigrels = logspace(-2,3,50); 
    cstars = logspace(-1,4,30);
else
    sigrels = q.sigrels;
    cstars = q.cstars;
end

Ksqrt = sqrt(KA*KB);


[SR,CS] = meshgrid(sigrels,cstars);

fractions = NaN*ones(length(sigrels),length(cstars)).';
maxCImean = NaN*ones(length(sigrels),length(cstars)).';
deltaCI = maxCImean;
percdeltaCI = maxCImean;

shiftfactor = q.shiftfactor;

for kk = 1:length(cstars)
    for jj = 1:length(sigrels)
        
        probc = @(logc) exp(-((logc-log(CS(kk,jj))).^2)/(2*SR(kk,jj)^2))./(sqrt(2*pi)*SR(kk,jj));

        CI = @(c,fa) CIf(feval(q.SNRfuncname,c,fa,q));

        % One issue with this integral is that its support can be a little
        % non-obvious... we force quadgk to include all of the relevant
        % points by including some waypoints, spanning a region around
        % cstar, and ensuring that we are also choosing points close to
        % Ksqrt
        waypoints1 = linspace(log(CS(kk,jj))-shiftfactor*SR(kk,jj),log(CS(kk,jj))+shiftfactor*SR(kk,jj),125);
        waypoints2 = linspace(log(Ksqrt)-shiftfactor,log(Ksqrt)+shiftfactor,25);
        waypoints3 = linspace(-shiftfactor,+shiftfactor,25);
        waypoints = sort(unique([waypoints1 waypoints2 waypoints3]));        
        
        negCImean = @(fa) quadgk( @(logc) -CI(exp(logc),fa).*probc(logc),-Inf,Inf,'Waypoints',waypoints);

        [f,negCIval] = fminbnd(negCImean,0,1,opts);
        
        fractions(kk,jj) = f;
        maxCImean(kk,jj) = -negCIval;
        ciallA = -negCImean(1);
        ciallB = -negCImean(0);
        deltaCI(kk,jj) = -negCIval - max(ciallA,ciallB);
        percdeltaCI(kk,jj) = (-negCIval - max(ciallA,ciallB))/(max(ciallA,ciallB));
    end
    kk
end
end