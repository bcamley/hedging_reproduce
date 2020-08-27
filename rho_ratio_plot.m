% SNR_ratio_plot
% Plots the SNR as a function of rho -- Fig. 3 of the paper
lw = 8;
rhos = logspace(0,3,50);

q = struct;
q.KA = 1;
q.KB = 1e3;
q.g = 0.05;
q.nr = 5e4;
q.kminusT = 2;

c = sqrt(q.KA*q.KB);
fracA = 0.5;

SNR_MLE = NaN*ones(size(rhos));
SNR_naive = NaN*ones(size(rhos));
SNR_naive_fastonly = NaN*ones(size(rhos));

for j = 1:length(rhos)
    q.kminBAratio = rhos(j);
    SNR_MLE(j) = SNR_two_timeaverage(c,fracA,q);
    SNR_naive(j) = SNR_two_timeaverage_naive(c,fracA,q);
    qfastonly = q;
    qfastonly.nr = q.nr/2;
    SNR_naive_fastonly(j) = SNR_two_timeaverage_naive(c,0,qfastonly);
end


%[SNR_MLE,SNR_naive,SNR_naive_fastonly] = naive_MLE_as_function_rho(f,KD,rhos,c0);
clf
hold on
plot(rhos,SNR_MLE,'x-','LineWidth',lw*1.2);
plot(rhos,SNR_naive,'LineWidth',lw*0.7);
plot(rhos,SNR_naive_fastonly,'o-','LineWidth',lw*0.5);
plot(rhos,SNR_MLE/2,'--','LineWidth',lw);
xlabel('\rho')
ylabel('Signal-to-noise ratio')
set(gca,'FontSize',48,'LineWidth',2);
set(gca,'xscale','log','yscale','log')
%set(gca,'YTick',[0.01 0.1 1 10])
box on
axis tight
legend('ERT','Naive','Naive (fast only)','ERT / 2')