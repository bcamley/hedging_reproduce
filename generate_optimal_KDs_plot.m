% generate optimal KD branching plot
% this is Fig. 2 in the paper
q = struct;
q.g = 0.05;
q.nr = 5e4;
q.CIfunc = 'ci';
%q.CIfunc = 'hill';
%q.hilln = 1;
q.maxNtypes = 7;
q.shiftfactor = 100; % used in the numerical integrals
q.SNRfuncname = 'SNR_multi_snapshot';


[logKvecs,maxCImean,SR,fracs] = get_optimal_KDs_and_fractions_symmetric(q);

CI_difference_threshold = 0.01;   % choose the receptor configuration with the lowest number of types, within those configurations within 0.01 of the optimal

KD_as_func_environment_plot(SR,logKvecs,maxCImean,CI_difference_threshold,[0.2 1 5],fracs);