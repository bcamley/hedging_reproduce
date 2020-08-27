function [fractions,maxCImean,deltaCI,percdeltaCI,SR,CS,nr_float] = find_fraction_function_allow_N_to_float(q)

opts = optimset('TolFun',1e-6,'TolX',1e-5,'MaxFunEvals',1e3);

switch(q.CIfunc)
    case 'ci'
        CIf = @(SNR) ci_trunc(sqrt(SNR));
    case 'hill'
        CIf = @(SNR) (SNR.^q.hilln)./(1+SNR.^q.hilln);
end

KA = q.KA;
KB = q.KB;

if(~isfield(q,'sigrels'))
    %sigrels = logspace(-2,3,50); 
    %cstars = logspace(-1,4,30);
    sigrels = logspace(-2,3,25); 
    cstars = logspace(-1,4,15); % this is a bit slower, so by default we'll do less computation
else
    sigrels = q.sigrels;
    cstars = q.cstars;
end

Ksqrt = sqrt(KA*KB);


[SR,CS] = meshgrid(sigrels,cstars);

fractions = NaN*ones(length(sigrels),length(cstars)).';
nr_float = NaN*ones(length(sigrels),length(cstars)).';

maxCImean = NaN*ones(length(sigrels),length(cstars)).';
deltaCI = maxCImean;
percdeltaCI = maxCImean;

shiftfactor = q.shiftfactor;

for kk = 1:length(cstars)
    for jj = 1:length(sigrels)
        
        probc = @(logc) exp(-((logc-log(CS(kk,jj))).^2)/(2*SR(kk,jj)^2))./(sqrt(2*pi)*SR(kk,jj));

        %CI = @(c,fa) CIf(feval(q.SNRfuncname,c,fa,q));

        % One issue with this integral is that its support can be a little
        % non-obvious... we force quadgk to include all of the relevant
        % points by including some waypoints, spanning a region around
        % cstar, and ensuring that we are also choosing points close to
        % Ksqrt
        waypoints1 = linspace(log(CS(kk,jj))-shiftfactor*SR(kk,jj),log(CS(kk,jj))+shiftfactor*SR(kk,jj),125);
        waypoints2 = linspace(log(Ksqrt)-shiftfactor,log(Ksqrt)+shiftfactor,25);
        waypoints3 = linspace(-shiftfactor,+shiftfactor,25);
        waypoints = sort(unique([waypoints1 waypoints2 waypoints3]));        
        
        negCImean = @(fa,nr_new) quadgk( @(logc) -CI_with_varying_nr(exp(logc),fa,q,nr_new,CIf).*probc(logc),-Inf,Inf,'Waypoints',waypoints);
        negCIpenalty = @(xn) negCImean(0.5+0.5*tanh(xn(1)),xn(2)).*exp(-max(xn(2)-q.nr,0)/q.npenalty); % use this to constrain fa between 0 and 1 while x varies everywhere
        [xn,negCIval] = fminsearch(negCIpenalty,[0 q.nr],opts);
        x = xn(1);
        
        f = 0.5+0.5*tanh(x); 
        
        fractions(kk,jj) = f;
        nr_float(kk,jj) = xn(2);
        maxCImean(kk,jj) = -negCIval;
        ciallA = NaN;
        ciallB = NaN;
        deltaCI(kk,jj) = -negCIval - max(ciallA,ciallB);
        percdeltaCI(kk,jj) = (-negCIval - max(ciallA,ciallB))/(max(ciallA,ciallB));
    end
    kk
end
end

function CI = CI_with_varying_nr(c,fa,q,nr_new,CIf)
    qnew = q;
    qnew.nr = nr_new;
    CI = CIf(feval(qnew.SNRfuncname,c,fa,qnew));
end