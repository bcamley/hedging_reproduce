% script to generate figures showing tradeoffs between two receptor types.
% Run to reproduce Fig. 1, 4, 7-9
% Note there are sections for each of these figures which may be run
% individually
qd = struct;  % sets default parameters
qd.CIfunc = 'ci';
qd.KA = 1;
qd.KB = 1e3;
qd.g = 0.05;
qd.nr = 5e4;
qd.shiftfactor = 100; % numerical parameter for specifying the integration waypoints, shouldn't matter much
                      % if interested in this, dig into
                      % find_fraction_function.m
% qd.hilln = 1; % if CIfunc = 'hill', this controls the sharpness of the
% Hill function

plotchecks = true;

%% generate tradeoff plot for snapshots
qsnap = qd;
qsnap.SNRfuncname = 'SNR_two_snapshot';

[fractions_snap,maxCImean_snap,deltaCI_snap,percdeltaCI_snap,SR,CS] = find_fraction_function(qsnap);
timeaverage = false;
plot_tradeoff_single(CS,SR,fractions_snap,qsnap,timeaverage);
figure
if(plotchecks)
    check_plot(fractions_snap,maxCImean_snap,deltaCI_snap,percdeltaCI_snap,CS,SR,qsnap,timeaverage);
    figure
end

%% Receptor number floating 
% Do a tradeoff plot for A,B receptors but allow the number of receptors
% to change, but with a penalty expressing added costs to the cell for
% having more receptors than typical
qsnapfloat = qsnap;
qsnapfloat.npenalty = 50*qsnap.nr; % CI_m suppressed exponentially with this decay constant
                                   % i.e. CI_m is penalized to 1/e of its original value
                                   % at npenalty. We choose this to be
                                   % pretty large -- that way you can
                                   % increase beyond the nominal value
                                   % without much effect
[fractions_snap_float,maxCImean_snap_float,deltaCI_snap_float,percdeltaCI_snap_float,SR,CS,nr_float] = find_fraction_function_allow_N_to_float(qsnapfloat);

subplot(1,2,1);
pcolor_better(CS,SR,fractions_snap_float); set(gca,'xscale','log'); set(gca,'yscale','log');
shading interp; colorbar;
hold on
yr = ylim;
Ksqrt = sqrt(qsnapfloat.KA*qsnapfloat.KB);
plot(Ksqrt*ones(1,1e2),linspace(yr(1),yr(2),1e2),'r','LineWidth',6);
title('Optimal A receptor fraction');
set(gca,'FontSize',18,'LineWidth',4)
xlabel('Environment concentration c_* [nM]','FontSize',18)
ylabel('Environment variation \sigma_\mu','FontSize',18)
set(gca,'XTick',[1 100 1e4])
text(Ksqrt,19,'(K_A K_B)^{1/2}','color','r','FontSize',24)

subplot(1,2,2);
pcolor_better(CS,SR,nr_float); set(gca,'xscale','log'); set(gca,'yscale','log');
shading interp; colorbar;
title('Optimal receptor number');
set(gca,'FontSize',18,'LineWidth',4)
xlabel('Environment concentration c_* [nM]','FontSize',18)
ylabel('Environment variation \sigma_\mu','FontSize',18)
set(gca,'XTick',[1 100 1e4])


%% generate tradeoff plots for time average with different values of ratio rho
qta = qd;
%qta.SNRfuncname = 'SNR_two_timeaverage_nonasymptotic';
qta.SNRfuncname = 'SNR_two_timeaverage';
qta.kminusT = 2;               % k-A times the averaging time T
%qta.nr = 5e4;
qta.nr = 1e4;                  % benefit of time averaging is that SNR is appreciable at lower nr

ratios = [1 10 100 1000];      % ratios of off rates - rho
fracs_timeaverage = cell(size(ratios));

timeaverage = true;

for j = 1:length(ratios)

    qta.kminBAratio = ratios(j);
    [fractions,maxCImean,deltaCI,percdeltaCI,SR,CS] = find_fraction_function(qta);
    if(plotchecks)
        figure
        check_plot(fractions,maxCImean,deltaCI,percdeltaCI,CS,SR,qta,timeaverage);
        drawnow
    end
    fracs_timeaverage{j} = fractions;
    
end
figure
plot_timeaverage_ratios_noncamel(fracs_timeaverage,ratios,qta,CS,SR);
% 
% %% generate tradeoff plots for time average with rho < 1
% % this will give all-A - commented out because it's not very interesting!
% qta = qd;
% %qta.SNRfuncname = 'SNR_two_timeaverage_nonasymptotic';
% qta.SNRfuncname = 'SNR_two_timeaverage';
% qta.kminusT = 2;               % k-A times the averaging time T
% %qta.nr = 5e4;
% qta.nr = 1e4;                  % benefit of time averaging is that SNR is appreciable at lower nr
% 
% ratios2 = [0.001 0.01 0.1 1];      % ratios of off rates - rho
% fracs_timeaverage2 = cell(size(ratios2));
% 
% timeaverage = true;
% 
% for j = 1:length(ratios2)
% 
%     qta.kminBAratio = ratios2(j);
%     [fractions2,maxCImean2,deltaCI2,percdeltaCI2,SR,CS] = find_fraction_function(qta);
%     if(plotchecks)
%         figure
%         check_plot(fractions2,maxCImean2,deltaCI2,percdeltaCI2,CS,SR,qta,timeaverage);
%         drawnow
%     end
%     fracs_timeaverage2{j} = fractions2;
%     
% end
% figure
% plot_timeaverage_ratios_noncamel(fracs_timeaverage2,ratios2,qta,CS,SR);