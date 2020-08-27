function h = plot_tradeoff_single(CS,SR,fractions,q,timeaverage)
KA = q.KA;
KB = q.KB;
Ksqrt = sqrt(KA*KB);

clf
h = contourf_better(CS,SR,fractions,linspace(0,1,16)); set(gca,'xscale','log'); set(gca,'yscale','log');
shading interp; caxis([0 1]);
set(gca,'FontSize',40,'LineWidth',4)

hh = colorbar; ylabel(hh,'Fraction of A receptors','FontSize',42);
xlabel('Environment concentration c_* [nM]','FontSize',42)
ylabel('Environment variation \sigma_\mu','FontSize',42)
hold on
yr = ylim;
if(~timeaverage)
    plot(Ksqrt*ones(1,1e2),linspace(yr(1),yr(2),1e2),'r','LineWidth',6);
    plot(KA*ones(1,1e2),linspace(yr(1),yr(2),1e2),'--','color',[0.3 0.3 0.3],'LineWidth',3);
    text(KA,19,'K_A','color',[0.3 0.3 0.3],'FontSize',60)
    
    text(5,0.1,'All-A','color','k','FontSize',60)
    
    plot(KB*ones(1,1e2),linspace(yr(1),yr(2),1e2),'--','color',[0.7 0.7 0.7],'LineWidth',3);
    text(KB,19,'K_B','color',[0.7 0.7 0.7],'FontSize',60)
    
    text(100,0.1,'All-B','color','w','FontSize',60)
    text(sqrt(KA*KB),19,'(K_A K_B)^{1/2}','color','r','FontSize',60)
else
    
    
end

box on

