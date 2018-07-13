        function [xs_est,errs] = lintr_bns(ys,as,m,n,k,var_ep)
%
%
%                            description:
%
%   This code performs singular value shrinkage on the linearly transformed
%   input matrix ys using the D-transform by backprojecting the observations,
%   normalizing, and then shrinking.
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
        [ys3,var_ep3] = lintr_gen2std_color(ys,as,m,n,k,var_ep);
        [xs_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_stiel(ys3,m,n,k);


        end
%
%
%
%
%
