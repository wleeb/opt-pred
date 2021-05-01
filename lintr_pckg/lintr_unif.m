        function [xs_est,errs] = lintr_unif(ys,as,m,n,k)
%
%
%                            description:
%
%   This code performs singular value shrinkage on the linearly transformed
%   input matrix ys with missing values, under the assumption that the values
%   are missing uniformly at random. When used with missing data, this is the 
%   algorithm described in the OptShrink paper of Raj Rao Nadakuditi.
%
%
%                           input parameters:
%
%   ys - the m-by-n data matrix of reduced observations.
%   as - the m-by-n matrix of reduction matrices (each column contains the 
%      diagonal of the reduction matrix). For example, in a missing data
%      problem, ds would contain 1's in observed entries, 0's elsewhere.
%   m,n - the dimensionality and number of samples, respectively.
%   k - the rank of the signal.
%   var_ep - m-by-1 vector of variances of the noise
%
%
%                        output parameters:
%
%   xs_est - the estimated matrix. Size m-by-n, and rank (at most) k.
%      Optimally shrunk for minimizing asymptotic Frobenius norm loss.
%   errs - k-dimensional vector of (squared) error due to each spike; sum(errs)
%      is an estimate of the total error || xs_est - xs ||^2 / n (Frobenius norm)
%
%
        p_hat = mean(as(:))
        ys2 = ys / p_hat;
        [xs_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_stiel(ys2,m,n,k);


        end
%
