function [logKvecs,maxCImean,SR,fracs] = get_optimal_KDs_and_fractions_symmetric(q)


switch(q.CIfunc)
    case 'ci'
        CIf = @(SNR) ci_trunc(sqrt(SNR));
    case 'hill'
        CIf = @(SNR) (SNR.^q.hilln)./(1+SNR.^q.hilln);
end

if(~isfield(q,'sigrels'))
    sigrels = logspace(-1,1,50);
    cstar = 1;
else
    sigrels = q.sigrels;
    cstar = 1;
end

Ntypes = 1:q.maxNtypes;

[SR,NT] = meshgrid(sigrels,Ntypes);

logKvecs = cell(length(sigrels),length(Ntypes)).';           % generate logKD for each uncertainty, number of types
maxCImean = NaN*ones(length(sigrels),length(Ntypes)).';
fracs = cell(length(sigrels),length(Ntypes)).';

warning off % during minimization, can get some Inf/NaN at extreme values
for kk = 1:length(Ntypes)
    for jj = 1:length(sigrels)
        
        probc = @(logc) exp(-((logc-log(cstar)).^2)/(2*SR(kk,jj)^2))./(sqrt(2*pi)*SR(kk,jj));
        waypoints1 = linspace(log(cstar)-q.shiftfactor*SR(kk,jj),log(cstar)+q.shiftfactor*SR(kk,jj),100);
        waypoints2 = linspace(-q.shiftfactor,q.shiftfactor,50);
        waypoints = sort(unique([waypoints1 waypoints2]));
        
        CI = @(c,logKvec,fracs) CIf(feval(q.SNRfuncname,c,logKvec,q,fracs));
        %if(mod(Ntypes(kk),2)==0)
        % The first half of z will be the logKDs, N/2 of them, with
        % reflection assumed, the second half the probabilities
        %   fh = (1:Ntypes(kk)/2); sh = ((Ntypes(kk)/2)+1):(Ntypes(kk));
        %   negCImean = @(z) quadgk( @(logc) -CI(exp(logc),[z(fh) -fliplr(z(fh)) z(sh) fliplr(z(sh))]).*probc(logc),-Inf,Inf,'Waypoints',waypoints);
        %else
        % with odd number, always makes sense to have logKD = 0
        %   fh = (1:Ntypes(kk)/2);
        %   sh_with_zero = ((Ntypes(kk)/2)+1):(Ntypes(kk)+1);
        %    sh_no_zero = ((Ntypes(kk)/2)+2):(Ntypes(kk)+1);
        %   negCImean = @(z) quadgk( @(logc) -CI(exp(logc),[z(fh) 0 -fliplr(z(fh)) z(sh_with_zero) fliplr(z(sh_no_zero))]).*probc(logc),-Inf,Inf,'Waypoints',waypoints);
        
        %end
        numv = floor(Ntypes(kk)/2); % number of variables in z's kD
        negCImean = @(z) quadgk( @(logc) -CI(exp(logc),symk(z(1:numv),Ntypes(kk)),symf(z(numv+1:end),Ntypes(kk))).*probc(logc),-Inf,Inf,'Waypoints',waypoints);
        %symk and symf functions symmetrize z
        
        opts = optimset('MaxFunEvals',4500*Ntypes(kk),'MaxIter',4500*Ntypes(kk),'TolX',1e-8,'TolFun',1e-8,'Display','off');
        
        
        logKstart_widespread = linspace(-15,15,Ntypes(kk));
        logKstart_defaultspread = linspace(-5,5,Ntypes(kk));
        logKstart_tinyspread = linspace(-0.1,0.1,Ntypes(kk));
        
        logKstartcell = {logKstart_widespread,logKstart_defaultspread,logKstart_tinyspread};
        %logKstartcell = {logKstart_defaultspread};
        
        xss = cell(size(logKstartcell));
        negCIvals = NaN*ones(size(logKstartcell));
        
        for ici = 1:length(logKstartcell)
            logKstart1 = logKstartcell{ici};
            %logKstart1 = linspace(-2,2,Ntypes(kk));
            
            if(mod(Ntypes(kk),2)==0)
                logKstart = sort(abs(logKstart1(1:Ntypes(kk)/2)));   % works when we have sorted, increasing logK
                fstart = ones(size(1:Ntypes(kk)/2)); fstart = fstart/sum(fstart);
            else
                %Ntypes(kk)
                logKstart = sort(abs(logKstart1(1:(Ntypes(kk)-1)/2)));
                fstart = ones(size(1:(Ntypes(kk)+1)/2)); fstart = fstart/sum(fstart);
            end
            [xs,negCIval] = fminsearch(negCImean,[logKstart fstart],opts);   % choose initial KDs linearly spaced
            xss{ici} = xs;
            negCIvals(ici) = negCIval;
        end
        [~,icibest] = min(negCIval);
        xs = xss{icibest};
        negCIval = negCIvals(icibest);
        
        logKunsort = symk(xs(1:numv),Ntypes(kk));
        fracunsort = symf(xs((numv+1):end),Ntypes(kk));
        [logKsort,sortind] = sort(logKunsort);
        fracsort = fracunsort(sortind);
        logKvecs{kk,jj} = logKsort;
        fracs{kk,jj} = fracsort;
        maxCImean(kk,jj) = -negCIval;
        
    end
    kk
end
warning on

end