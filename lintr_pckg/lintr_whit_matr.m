        function [xs_est,whts,errs,spec] = lintr_whit_matr(ys,as,m,m2,n,k,cov_ep)
%
%
%                            description:
%
%   This code performs singular value shrinkage on the linearly transformed
%   input matrix ys using the closed formula of Donoho and Gavish by back-
%   projecting, whitening, shrinking, and unnormalizing. The 
%   transformations are assumed to be matrices of the same size. Complex values
%   are allowed for all inputs.
%
%
%                           input parameters:
%
%   ys - the m2-by-n data matrix of reduced observations.
%   as - the m2-by-m-by-n cube of reduction matrices (each slice contains the 
%      the reduction matrix). For example, in an image deblurring
%      problem, the point-spread functions (in image, not Fourier, domain).
%   m - the dimension of the original signal
%   m2 - the dimension of the linearly transformed signal; typically m2 <= m 
%   n - the number of observations
%   k - the rank of the signal.
%   cov_ep - the covariance matrix of the noise; of size m2-by-m2
%
%
%                        output parameters:
%
%   xs_est - the estimated matrix. Size m-by-n, and rank (at most) k.
%   whts - m-by-m symmetric matrix defining the loss function minimized by the
%      shrinkage procedure
%   errs - k-dimensional vector of (squared) error due to each spike; sum(errs)
%      is an estimate of the total error || W*(xs_est - xs) ||^2 / n (Frobenius
%      norm) where W is the matrix whts
%
%

%
%        backproject and whiten the data
%
        [ys3,wmat,winv,as2_mean,as2_inv,spec] = lintr_gen2whi_matr(ys,as,...
            m,n,cov_ep);

%
%        apply the Donoho Gavish shrinker
%
        sig=1;
        [xs_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_white2(ys3,m,n,k,sig);

%
%        unnormalize and build weights matrix
%
        whts_inv = as2_inv*winv;
%%%        whts_inv2 = winv*as2_inv;
        xs_est = whts_inv * xs_est;
        whts = wmat*as2_mean;


        end
%
