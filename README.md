# ITEM-paper

<b>MATAB Code for the ITEM Paper, as submitted to NeuroImage 2019</b>

This code belongs to the paper on "Inverse Transformed Encoding Models" (ITEM) by Soch et al. (2019), publicly available from <i>bioRxiv</i> and currently under review at <i>NeuroImage</i>. It consists of two sub-folders, "Simulation" (Section 3, Figures 4-5 in the paper) and "Application" (Section 4, Figures 7-8 in the paper).

- Preprint: https://www.biorxiv.org/content/10.1101/610626v2
- Data: https://openneuro.org/datasets/ds002013 (Version 1.0.2)
- Code: https://github.com/JoramSoch/ITEM-paper
- Toolbox: https://github.com/JoramSoch/ITEM


<h3>Requirements</h3>

This code was developed and run using the following software:
- Windows 7 Enterprise N SP1 64-bit
- <a href="https://de.mathworks.com/help/matlab/release-notes-R2013b.html">MATLAB R2013</a> (Version 8.2)
- <a href="https://de.mathworks.com/products/statistics.html">MATLAB Statistics Toolbox</a> (Version 8.3)
- <a href="https://www.fil.ion.ucl.ac.uk/spm/software/spm12/">SPM12</a> (Revision 7219 as of 16/11/2017)
- <a href="https://github.com/JoramSoch/MACS">MACS Toolbox</a> (Version 1.3)
- <a href="https://github.com/JoramSoch/ITEM">ITEM Toolbox</a> (as on GitHub)


<h3>Simulation</h3>

This sub-folder contains the following code:
- `Simulation.m`: runs simulation, as described in Section 3 and Appendix E of the paper;
- `MD_unirnd.m`: samples random numbers from the continuous uniform distribution;
- `MD_matnrnd.m`: samples random numbers from the matrix normal distribution;
- `Figure_4.m` and `Figure_5.m`: generate Figures 4 and 5, as they appear in the paper.

Note that, other than in the paper, the simulation
- implements a fourth method, "least squares, whitened" (LS-W), which is based on whitening the trial-wise parameters estimates with `W = U^(-1/2)` (and then applying the two-sample t-test and logistic regression, respectively) instead of performing the whole inverse transformed encoding model approach (which implicitly also requires whitening the design matrix `T`). In the figures, LS-W results are given in cyan color. Because whitening only the trial-wise response amplitudes is not sufficient, LS-W is inferior to LS-T a.k.a. ITEM and these results were not included in the paper.
- produces more figures, some of which are just there to check that the simulation was done correct (Figures 1-2) or to explore several statistics of raw, separate and whitened trial-wise parameter estimates (Figures 3-6). Figures 7 and 8 are preliminary views of the results reported in the paper and can be reproduced by running `Figure_4.m` and `Figure_5.m`. Note that this requires to first run `Simulation.m` which produces a file `Simulation.mat` containing the simulation results to be plotted.

Findings denoted in the paper as “results not shown”, regarding the specificity of ITEM and other approaches, can be found in the sub-folder “Simulation/null_results/”. They can also be replicated by setting `μA = μB` in the simulation code (l. 35 in `Simulation.m`).

The simulation uses the routines `ttest2`, `mnrfit`, `normrnd` and `mvnrnd` from MATLAB’s Statistics Toolbox. Use of the normal distribution samplers `normrnd` and `mvnrnd` can be circumvented by using the supplied matrix normal sampler `MD_matnrnd.m` (which uses the MATLAB built-in `randn` and which is already implemented in place of `normrnd` and `mvnrnd`).


<h3>Application</h3>

For re-running analyses of the empirical data, you need to perform the following steps:
1. Create a folder on your computer that hereafter is referred to as the "study directory".
2. Download the empirical data set from https://openneuro.org/datasets/ds002013 and place it into a sub-folder of the study directory called "BIDS".
3. Download the application scripts from https://github.com/JoramSoch/ITEM-paper and place them into a sub-folder of the study directory called "tools".
4. Open MATLAB, set your current directory to this sub-folder "tools", edit the study directory in the `project_directories.m` and run this script.
5. Run the script `analyses_VisRec.m` located in the same folder. Ideally, run the code line by line to ensure that each step of the analysis succeeds.
6. When the entire analysis of the empirical data set has been performed, the scripts `Figure_7.m` and `Figure_8.m` from the sub-folder "VisRec_ITEM" of the tools directory should generate Figures 7 and 8, as they appear in the paper.

The sub-folder "VisRec_ITEM/LS-A_LS-S/" contains additional code that allows to combine LS-A and LS-S with support vector regression (SVR; requiring `libsvm`), to be called from the scripts `VisRec_ITEM_LS_A.m` and `VisRec_ITEM_LS_S.m` in the sub-folder "VisRec_ITEM". Because SVR makes the reconstruction performance of LS-A and LS-S even worse, compared with ITEM-style inversion, these results were not included in the paper.

Note that ROI images supplied with the BIDS-compatible data set on OpenNeuro are registered to the pre-processed data when using the pre-processing routines from <i>this</i> repository. If pre-processing routines are modified (e.g. to exclude spatial realignment or to include spatial normalization), the ROI images may <i>not</i> match the pre-processed data anymore in which case they must be resliced to the pre-processed data (i.e. the data to which statistical analysis will be applied).


<h3>Bonus: Graphical Abstract</h3>

<img src="https://github.com/JoramSoch/ITEM-paper/raw/master/Figure_GA.png" alt="Graphical Abstract" width=1000>

This graphical abstract illustrates the core idea of the paper: When multiplying the trial-wise design matrix with itself, weighted by the (inverse of the) scan-by-scan covariance matrix, this results in the (inverse of the) trial-by-trial covariance matrix which describes the distribution of the trial-wise parameter estimates. Click <a href="https://github.com/JoramSoch/ITEM-paper/raw/master/Figure_GA.pdf">here</a> for a PDF version.
