function h = check_plot(fractions,maxCImean,deltaCI,percdeltaCI,CS,SR,q,timeaverage)

subplot(1,4,1)
h=pcolor(CS,SR,fractions)
set(gca,'xscale','log','yscale','log');
xlabel('Environment concentration [nM]')
ylabel('Environment variation \sigma_\mu')
%title('Fraction of A receptors')
if(timeaverage)
    %title(sprintf('Fraction of A receptors; \\rho = %d',q.kminBAratio))
    title(sprintf('Optimal A fraction; $\\rho = %d$',q.kminBAratio),'interpreter','latex')

else
    %sprintf('Fraction of A receptors; \\rho = %d',q.kminBAratio)
    title('Optimal A fraction','interpreter','latex')
end
caxis([0 1]);
colorbar

subplot(1,4,2)
pcolor(CS,SR,maxCImean)
set(gca,'xscale','log','yscale','log');
xlabel('Environment concentration [nM]')
ylabel('Environment variation \sigma_\mu')
%title('Maximum mean CI')
title('Maximum $\overline{CI}$','interpreter','latex')
caxis([0 1]);
colorbar

subplot(1,4,3)
pcolor(CS,SR,deltaCI)
set(gca,'xscale','log','yscale','log');
xlabel('Environment concentration [nM]')
ylabel('Environment variation \sigma_\mu')
colorbar
%title('Change in \bar{CI} from all-A/all-B')
title('Change in $\overline{CI}$ due to hedging','interpreter','latex')

subplot(1,4,4)
pcolor(CS,SR,percdeltaCI)
set(gca,'xscale','log','yscale','log');
xlabel('Environment concentration [nM]')
ylabel('Environment variation \sigma_\mu')
colorbar
title('Fractional change in $\overline{CI}$ due to hedging','interpreter','latex')

for s = 1:4
    subplot(1,4,s)
    set(gca,'FontSize',12);
    set(gca,'XTick',[1 100 1e4]);
    shading interp
end
set(gcf,'Position',[66 333 1855 318]);
