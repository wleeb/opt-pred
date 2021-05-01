%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        This code contains the following user-callable functions for denoising
%        observations from the linearly transformed spiked model (LTSM). Read
%        the function descriptions for details.
%
%     lintr_fulshr - shrinkage on LTSM, general transformations
%     lintr_fulshr2 - shrinkage on LTSM, list of transformations
%     lintr_whit_matr - shrinkage on LTSM, general transformations
%     lintr_whit_matr2 - shrinkage on LTSM, list of transformations
%     lintr_mrows - shrinkage for transformed rows, estimated noise variance
%     lintr_bsn - backprojects, shrinks, normalizes
%     lintr_bns - backprojects, normalizes, shrinks
%     
%        The methods for these codes are described in the paper
%        ``Optimal Prediction in the Linearly Transformed Spiked Model'',
%        by E. Dobriban, W. Leeb, and A. Singer.
%
%     See the function ``run_...m'' for an example of how to run the test code,
%     necessary dependencies, etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
