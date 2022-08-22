# Revised_Huang_Model
Theta-burst stimulation (TBS), a patterned form of repetitive transcranial magnetic stimulation (TMS) for rapid neuromodulation, is widely used in both experimental brain research and clinical applications. A decade ago, the first calcium-dependent model was developed by Y.Z. Huang to explain the bidirectional effects of TBS. We discover, however, that the published computer code is not consistent with the model formulations in the [initial paper](https://doi.org/10.1016/j.clinph.2010.08.016). Further analysis confirms that the computer model with the index confusion was used for fitting the experimental data. The Huang's model has distinct three stages: calcium dynamics (Stage I), the substances dynamics (Stage II), and the final after-effect curves (Stage III). We modified the computer code of Stage II of the initial model. Hence, we re-calibrated the parameters of Stage II of the revised model by creating a non-convex optimisation problem. The revised model successfully predicts the polarities and time courses of the TBS-induced after-effects and has an advantage over the old model in prediction accuracy.

## Folder Summary
1. ExperimentalMeasurements: this folder contains all experimental data used to re-calibrate the revised model. These experimental data are digitalized from plots in articles:
    * Huang YZ, Edwards MJ, Rounis E, Bhatia KP, Rothwell JC. Theta burst stimulation of the human motor cortex. Neuron. 2005 Jan 20;45(2):201-6.
    * Gentner R, Wankerl K, Reinsberger C, Zeller D, Classen J. Depression of human corticospinal excitability induced by magnetic theta-burst stimulation: evidence of rapid polarity-reversing metaplasticity. Cerebral cortex. 2008 Sep 1;18(9):2046-53.
    * Huang YZ, Rothwell JC, Chen RS, Lu CS, Chuang WL. The theoretical model of theta burst form of repetitive transcranial magnetic stimulation. Clinical Neurophysiology. 2011 May 1;122(5):1011-8.
2. Functions: this folder contains all MATLAB functions used by: `FourPara_Calibration_LM.m`, `PlotFigures.m`, `PredictionAccuracyAnalysis.m`.
    * `peakM.m`: it contains Stages I and II of the revised model and calculates the net after-effect M(t) for a given TBS protocol;
    * `AfterEffectFun.m`: it uses the calculated net after-effect and calculates the inhibitory or facilitatory after-effect curve for a given TBS protocol at Stage III of the revised model;
    * `HuangModel_Old.m`: it is the original model initially proposed by Huang et al,;
    * `HuangModel_V2_modified.m`: it is the revised model;
    * `ModifiedCostFun.m`: it calculates the sum of squared residuals between the predictions and the corresponding experimental data;
    * `PlotFittingCurve.m`: it plots the after-effect curves for a given set of parameters of Stage II.

## Code Summary
This code project implements a non-convex optimisation problem.
1. `FourPara_Calibration_LM.m`: this code file aims to optimise the parameters of Stage II by using the Levenberg-Marquardt algorithm as implemented in Matlab (v2021b,
The Mathworks, USA).
2. `FourPara_Calibration_Swarm.m`: this code file aims to optimise the parameters of Stage II by using the Particle swarm algorithm as implemented in Matlab (v2021b, The Mathworks, USA).
3. `PlotFigures_PC.m`, `PlotFigures_AC.m`, `PlotFigures_noPC.m`: these code files aim to plot the after-effect curves for both initial and revised models for given TBS protocols. `_PC` is for protocols with prior contraction plots; `_AC` is for protocols with post contraction plots; `_noPC` is for protocols without prior contraction plots.
4. `PredictionAccuracyAnalysis.m`: this code file aims to analyse the prediction accuracy of the revised model in terms of root-mean-square error (RMSE), pseudo R squared (R^2) for nonlinear models, and the distribution of prediction errors.
5. `HuangModel_Initial.m`: this code file is the initial model code programmed by Huang et al,.
6. `HuangModel_Revised.m`: this code file is the revised model with the optimised parameters of Stage II.
