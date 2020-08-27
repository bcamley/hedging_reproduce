function h = plot_timeaverage_ratios_noncamel(fracs_as_func_ratio,ratios,q,CS,SR)

fs = 20;
lw = 8;

clf

KA = q.KA;
KB = q.KB;

h = cell(size(ratios));

for s = 1:4
subplot(2,4,s)

    title(sprintf('\\rho = %d',ratios(s)));
        if(round(ratios(s))-ratios(s)>1e-3)
            warning('Rounding label!!')
        end

q.kminBAratio = ratios(s);

%c = logspace(-3,6,100);
c = logspace(log10(min(CS(:))),6,100);

cA = c;
cB = c;


SNRboth = SNR_two_timeaverage(c,0.5,q);
SNRA = SNR_two_timeaverage(cA,1,q);
SNRB = SNR_two_timeaverage(cB,0,q);

%cla
hold on


plot(cA,SNRA,'-','LineWidth',lw/1.5,'MarkerSize',24,'color',[0.93 0.69 0.13]);
plot(cB,SNRB,'-','LineWidth',lw/1.2,'MarkerSize',24,'color',[0.85 0.33 0.1]);
plot(c,SNRboth,'-','LineWidth',lw,'color',[0 0.45 0.74]);
set(gca,'xscale','log','yscale','log');
%box on
set(gca,'LineWidth',2)

if(s==2)
xlabel('Concentration [nM]');
end
if(s==1)
ylabel('Signal-to-noise ratio');
end
set(gca,'FontSize',fs)
set(gca,'XTick',[1e-2 1 100 1e4])
axis tight
xlim([min(c(:)) 1e4])
switch(s)
    case 1
        ylim([1e-3 1e2]);
    case 2
        ylim([1e-2 1e3]);
    case 3
        ylim([1e-2 1e4]);
    case 4
        ylim([1e-2 1e5]);
end
%ylim([1e-3 2]*ratios(s))
box on
end

for ll = 1:length(ratios)
    subplot(2,4,ll+4)
    
    kminBAratio = ratios(ll);
    cbalance = (KB-KA*kminBAratio)/(kminBAratio-1);
    
    h{ll} = contourf_better(CS,SR,fracs_as_func_ratio{ll},linspace(0,1,15)); set(gca,'xscale','log'); set(gca,'yscale','log');
    %set(gca,'XTick',[0.1 1 10 100 1e3 1e4])
    set(gca,'XTick',[1 100 1e4])
    set(gca,'YTick',[0.01 1 100])
    shading interp; caxis([0 1]);

    if(ll==2)
        xlabel('Environment concentration c_* [nM]','FontSize',fs)
    end
    
    if(ll==1)
        ylabel('Environment variation \sigma_\mu','FontSize',fs)
    end
    hold on
    yr = ylim;
    
    
    plot(cbalance*ones(1,1e2),linspace(yr(1),yr(2),1e2),'w--','LineWidth',3);
    
    if( ll == 4)
        hh = colorbar; ylabel(hh,'Fraction of A receptors','FontSize',fs);
    end
    
    box on
    set(gca,'FontSize',fs,'LineWidth',2)

    drawnow
end
