        function [ss_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_white2(ys,m,n,k,sig)
%
%
%        This code performs singular value shrinkage on spiked data with 
%        white noise of specified size. It uses the optimal Frobenius shrinker 
%        of Donoho and Gavish; unlike the code svshr_white, this code automatically 
%        rescales for a specified noise variance.
%
%
%                             input parameters:
%
%   ys - the m-by-n matrix of data. The data model is iid columns of dimension m. It
%      is assumed that the columns are NOT divided by sqrt(n).
%   m,n - the dimensions of the data. m is intepreted as the dimension of each 
%      observation, while n is the number of observations.
%   k - the rank of the low-rank signals.
%   sig - the standard deviation of the noise, which is assumed to be white noise.
%
%
%                             output parameters:
%
%   ss_est - the m-by-n estimated data matrix. Rank is no more than k.
%   s_op - the k-dimensional vector of estimated population singular values 
%      (also the optimal singular values, in operator norm)
%   s_fr - the k-dimensional vector of optimal singular values for Frobenius norm loss
%   cos_out - the k-dimensional vector of estimated cosines between the empirical and 
%      population PCs (i.e. the the outer singular vectors u)
%   cos_inn - the k-dimensional vector of estimated cosines between the empirical and 
%      population inner singular vectors v
%   uy - the m-by-k matrix of empirical PCs (i.e. outer singular vectors).
%   sy - the k-dimensional vector of singular values of the data matrix, after 
%      scaling by sig and sqrt(n)
%   vy - n-by-k matrix of empirical inner singular vectors
%   errs - the k-dimensional vector of estimated errors from each component
%
%
        gam = m/n;

        s_fr = zeros(1,k);
        s_fr88 = zeros(1,k);
        s_op = zeros(1,k);
        cos_out = zeros(1,k);
        cos_inn = zeros(1,k);
%
        [uy,sy,vy] = svshr_svdsmartc(ys,m,n,k);
        sy = sy / sqrt(n) / sig;

        for i=1:k
%
        rlam=sy(i)^2;
        [ell,cos_out(i),cos_inn(i)] = svshr_emp2pop_white(rlam,gam);
        s_op(i) = sqrt(ell);
%
        s_fr(i) = s_op(i) * cos_out(i) * cos_inn(i);
    end

        ells = s_op.^2;
        errs = svshr_errs(ells,cos_inn,cos_out,k);
        s_op = s_op * sig;
        s_fr = s_fr * sig;

        ells = ells*sig^2;
        errs = errs*sig^2;

        s_fr17 = s_fr * sqrt(n);
        ss_est = uy(:,1:k)*diag(s_fr17)*vy(:,1:k)';

        end
%
%
%
%
%
