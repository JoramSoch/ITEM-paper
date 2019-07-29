# ITEM-paper

<b>MATAB Code for the ITEM Paper, as submitted to NeuroImage 2019</b>

The code belongs to the paper on "Inverse Transformed Encoding Models" (ITEM) by Soch et al. (2019), <a href="https://www.biorxiv.org/content/10.1101/610626v1">publicly available</a> from <i>bioRxiv</i> and currently under review at <i>NeuroImage</i>. It consists of two sub-folders, “Simulation” (Section 3, Figures 4-5 in the paper) and “Application” (Section 4, Figures 6-8 in the paper).


<h3>Requirements</h3>

This code was developed and run under the following requirements:
- Windows 7 Enterprise N SP1 64-bit
- <a href="https://de.mathworks.com/help/matlab/release-notes-R2013b.html">MATLAB R2013</a> (Version 8.2)
- <a href="https://de.mathworks.com/products/statistics.html">MATLAB Statistics Toolbox</a> (Version 8.3)
- <a href="https://www.fil.ion.ucl.ac.uk/spm/software/spm12/">SPM12</a> (Revision 7219 as of 16/11/2017)
- <a href="https://github.com/JoramSoch/ITEM">ITEM Toolbox</a> (available from GitHub)


<h3>Simulation</h3>

This sub-folder contains the following code:
- `Simulation.m`: runs simulation, as described in Section 3 and Appendix E of the paper;
- `MD_unirnd.m`: samples random numbers from the continuous uniform distribution;
- `MD_matnrnd.m`: samples random numbers from the matrix normal distribution;
- `Figure_4.m` and `Figure_5.m`: generate Figures 4 and 5, as they appear in the paper.

Note that, other than in the paper, the simulation
- implements a fourth method, “least squares, whitened” (LS-W), which is based on whitening the trial-wise parameters estimates with `W = U^(-1/2)` (and then applying the two-sample t-test and logistic regression, respectively) instead of performing the whole inverse transformed encoding model approach (which implicitly also requires whitening the design matrix `T`). In the figures, LS-W results are given in cyan color. Because whitening only the trial-wise response amplitudes is not sufficient, LS-W is inferior to LS-T a.k.a. ITEM and these results were not included in the paper.
- produces more figures, some of which are just there to check that the simulation was done correct (Figures 1-2) or to explore several statistics of raw, separate and whitened trial-wise parameter estimates (Figures 3-6). Figures 7 and 8 are preliminary views of the results reported in the paper and can be reproduced by running `Figure_4.m` and `Figure_5.m`. Note that this requires to first run `Simulation.m` which produces a file `Simulation.mat` containing the simulation results to be plotted.

The simulation uses the routines `ttest2`, `mnrfit`, `normrnd` and `mvnrnd` from MATLAB’s Statistics Toolbox. Use of the normal distribution samplers `normrnd` and `mvnrnd` can be circumvented by using the supplied matrix normal sampler `MD_matnrnd.m` (which uses the MATLAB built-in `randn` and which is already implemented in place of `normrnd` and `mvnrnd`).
