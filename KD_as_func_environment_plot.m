function [h,logKopt,fracopt] = KD_as_func_environment_plot(SR,logKvecs,maxCImean,deltaCIthreshold,sigplots,fracs)
checkplot = false;
%clf
hold on

optNtrue = NaN*ones(1,size(maxCImean,2));
optNcost = NaN*ones(1,size(maxCImean,2));

Ntypes = 1:size(maxCImean,1);

for s = 1:size(maxCImean,2)
    [CImax,kk] = max(maxCImean(:,s));
    optNtrue(s) = Ntypes(kk);
    closeenough = abs(maxCImean(:,s)-CImax)<deltaCIthreshold;
    optNcost(s) = min(Ntypes(closeenough)); % the smallest n given nearly-optimal performance
end

logKopt = cell(1,size(maxCImean,2));
fracopt = cell(1,size(maxCImean,2));
for s = 1:size(maxCImean,2)
    correct = find(optNcost(s)==Ntypes);
    lkds = logKvecs{correct,s};
    try
        fracopts = fracs{correct,s};
        fracopts = abs(fracopts)/sum(abs(fracopts));
    catch err
        getReport(err)
        fracopts = ones(size(lkds));
    end
    %plot(SR(1,s)*ones(size(lkds)),lkds,'ko');
    logKopt{s} = lkds; % this assumes cstar = 1
    fracopt{s} = fracopts; % this assumes cstar = 1
end

%fracopt
cmap = bone(1024);

for s = 1:length(logKopt)
    lkds = logKopt{s};
    fo = fracopt{s};
    fo = abs(fo) / sum(abs(fo));
    for ss = 1:length(fo)
        %h(s) = plot(SR(1,s),exp(lkds(ss)),'ko','MarkerSize',40*fo(ss),'MarkerFaceColor',[1 0.4 0.7]);
        h(s) = plot(SR(1,s),exp(lkds(ss)),'o','MarkerSize',32*sqrt(fo(ss)),'MarkerFaceColor',interp1(linspace(0,1,1024),cmap,fo(ss)),'LineWidth',1,'Color',[0.5 0.5 0.5]);
        % scaling markersize as 40*sqrt(fo(ss)) = area of marker is
        % proportional to fraction of receptors that are this species
    end
end

set(gca,'xscale','log','yscale','log');
xlabel('Environment concentration uncertainty \sigma_\mu')
ylabel('Optimal receptor K_D / c_*')
set(gca,'FontSize',40)

yr = ylim;
cf = logspace(log10(yr(1)),log10(yr(2)),1e3);
logc = log(cf);
shrinkfactor = 0.2;
for j = 1:length(sigplots)
    plot(sigplots(j)*exp(-shrinkfactor*exp(-(logc.^2)/(2*sigplots(j)^2))),cf,'LineWidth',4);
    
end

box on
set(gca,'LineWidth',2);

if(checkplot)
    figure
    clf

for jj = 1:length(Ntypes)
    logKall = cell2mat(logKvecs(jj,:).');
    subplot(3,3,jj)
    hold on
    for j = 1:Ntypes(jj)
        plot(SR(1,:),logKall(:,j));
    end
    set(gca,'xscale','log')
    xlabel('Environment concentration uncertainty')
    ylabel('Receptor log K_D')
end

figure
plot(SR.',maxCImean.')
set(gca,'xscale','log')
xlabel('Environment concentration uncertainty')
ylabel('CI')
end

end