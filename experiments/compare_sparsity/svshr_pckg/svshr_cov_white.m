        function [cov_est,ells,ells_fr,cos_out,cos_inn,errs] = ...
            svshr_cov_white(ys,m,n,k,sig)
%
%        This code applies the optimal eigenvalue shrinker for Frobenius loss
%        to estimate the covariance of data from the spiked covariance model
%        with white noise of specified variance.
%
%
%                                input parameters:
%
%   ys - the m-by-n matrix of data. The data model is iid columns of dimension m. 
%      It is assumed that the columns are NOT divided by sqrt(n).
%   m,n - the dimensions of the data. m is intepreted as the dimension of each 
%      observation, while n is the number of observations.
%   k - the rank of the low-rank signals.
%   sig - the standard deviation of the noise, which is assumed to be white noise.
%
%
%                                output parameters:
%
%   cov_est - the m-by-m symmetric estimated covariance matrix. The rank
%      does not exceed k.
%   ells - the k-dimensional vector of estimated population eigenvalues.
%   ells_fr - the k-dimensional vector of optimal eigenvalues for Frobenius loss 
%      (and the eigenvalues of cov_est).
%   cos_out - the k-dimensional vector of estimated cosines between the empirical and 
%      population PCs (i.e. the the outer singular vectors u)
%   cos_inn - the k-dimensional vector of estimated cosines between the empirical and 
%      population inner singular vectors v
%   errs - the k-dimensional vector of estimated errors from each component
%
%
        gam = m/n;

        [uy,sy,vy] = svshr_svdsmartc(ys,m,n,k);
        sy = sy / sig;


        for i=1:k
%
        rlam = sy(i)^2 / n;
        [ells(i),cos_out(i),cos_inn(i)] = svshr_emp2pop_white(rlam,gam);
        ells(i) = ells(i)*sig^2';
        ells_fr(i) = ells(i)*cos_out(i)^2;
%
        errs(i) = (1 - cos_out(i)^4) * ells(i)^2;
    end

        cov_est = uy * diag(ells_fr) * uy';


        end
%
%
%
%
%
