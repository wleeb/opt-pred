        function [ss_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_stiel(ys,m,n,k)
%
%        This code performs singular value shrinkage using the sample Stieltjes
%        transform to estimate the statistics of the noise in the limit. This 
%        is exactly the OptShrink procedure of Raj Rao (without missing data), and 
%        it gives identical results to his optshrink.m code (but is faster)
%
%
%                             input parameters:
%
%   ys - the m-by-n matrix of data. The data model is iid columns of dimension m. 
%      It is assumed that the columns are NOT divided by sqrt(n).
%   m,n - the dimensions of the data. m is intepreted as the dimension of each 
%      observation, while n is the number of observations.
%   k - the rank of the low-rank signals.
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
%
        s_fr = zeros(1,k);
        s_fr88 = zeros(1,k);
        s_op = zeros(1,k);
        cos_out = zeros(1,k);
        cos_inn = zeros(1,k);
%
        [uy,sy,vy] = svshr_svdsmartc(ys,m,n,min(m,n));
        sy = sy / sqrt(n);


        bedge=0;
        for i=1:k
%
        [s_op(i),cos_out(i),cos_inn(i)] = svshr_emp2pop_stiel(sy,m,...
           n,k,i,bedge);

        s_fr(i) = s_op(i) * cos_out(i) * cos_inn(i);
%%%        s_fr88(i) = svshr_emp2frob(sy,m,n,k,i);
    end

        ells = s_op.^2;
        errs = svshr_errs(ells,cos_inn,cos_out,k);
        err_hat = sum(errs);

%%%        chk0 = norm(s_fr - s_fr88)

        s_fr17 = s_fr * sqrt(n);
        ss_est = uy(:,1:k)*diag(s_fr17)*vy(:,1:k)';

        end
%
%
%
%
%
