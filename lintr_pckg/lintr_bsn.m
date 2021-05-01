        function [xs_est,whts,errs] = lintr_bsn(ys,as,m,n,k,var_ep)
%
%
%                            description:
%
%   This code performs singular value shrinkage on the linearly transformed
%   input matrix ys using the D-transform by backprojecting the observations,
%   shrinking, and normalizing at the end.
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
%   whts - m-by-1 vector, containing the mean projection-backprojection matrix.
%   errs - k-dimensional vector of (squared) error due to each spike; sum(errs)
%      is an estimate of the total error || W*(xs_est - xs) ||^2 / n (Frobenius norm)
%      where W is diag(whts).
%

        [ys2,var_ep2,as2_mean] = lintr_gen2backp(ys,as,m,n,k,var_ep);
        [xs_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_stiel(ys2,m,n,k);
%
        xs_est = diag(1./as2_mean) * xs_est;
        whts = as2_mean;
        
        end
%
