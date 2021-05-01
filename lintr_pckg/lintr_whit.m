        function [xs_est,whts,errs] = lintr_whit(ys,as,m,n,k,var_ep)
%
%
%                            description:
%
%   This code performs singular value shrinkage on the linearly transformed
%   input matrix ys using the closed formula of Donoho and Gavish by back-
%   projecting, whitening, shrinking, and unnormalizing. The 
%   transformations are assumed to be diagonal. Complex values
%   are allowed for all inputs.
%
%
%                           input parameters:
%
%   ys - the m-by-n data matrix of reduced observations.
%   as - the m-by-n matrix of reduction matrices (each column contains the 
%      diagonal of the reduction matrix).
%   m,n - the dimensionality and number of observations, respectively
%   k - the rank of the signal.
%   var_ep - m-by-1 vector of variances of the noise
%
%
%                        output parameters:
%
%   xs_est - the estimated matrix. Size m-by-n, and rank (at most) k.
%   whts - m-by-1 vector defining the diagonal matrix with weights defining
%      the loss function
%   errs - k-dimensional vector of (squared) error due to each spike; sum(errs)
%      is an estimate of the total error || W*(xs_est - xs) ||^2 / n (Frobenius
%      norm) where W is diag(whts)
%
%
%
%        backproject the data
%
        [ys2,var_ep2,as2_mean] = lintr_gen2backp(ys,as,m,n,k,var_ep);
        wvals = sqrt(1./var_ep2);
        winv = 1./wvals;
%
%        now whiten the effective noise
%
        ys3 = repmat(wvals,1,n) .* ys2;
%
%        shrink this using Donoho-Gavish
%
        sig=1
        [xs_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_white2(ys3,m,n,k,sig);
%
%        finally, unnormalize
%
        whts_inv = (1./as2_mean) .* winv;
        xs_est = repmat(whts_inv,1,n) .* xs_est;
        whts = 1./whts_inv;

        end
%
