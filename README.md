# Code to reproduce "Chemotaxis in uncertain environments: hedging bets with multiple receptor types"
## Austin Hopkins and Brian Camley

This code is hacked-together scientific code; it is not pretty, and could certainly be more efficient, but it serves to illustrate some of the detailed methods that cannot be captured within a paper and allow for easier reproduction of the results. We cannot guarantee that it performs as predicted if you try to change things without understanding them! In this sense, it is released in the spirit of the CRAPL license.

Paper available at: https://arxiv.org/abs/2002.10441

1. To reproduce the figures showing tradeoffs between type-A and type-B receptors (Fig. 1, Fig. 4, and Fig. 7-9), run `generate_tradeoff_plots_new.m` - this takes about 15 minutes on my laptop. The script sections can be run individually using Run Section as well. 

2. To reproduce the figure showing the optimal evolved receptor configuration as a function of uncertainty (Fig. 2), run `generate_optimal_KDs_plot`. This takes about 20 minutes. 

3. To reproduce the plot of signal-to-noise ratios, run `rho_ratio_plot.m` - this should be fast.

This code has been tested in Matlab 2017b.
