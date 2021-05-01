        function [xs_est,whts,errs,spec,xmean] = lintr_fulshr2(ys,as_uniq,ias,nas,...
            m,m2,n,k,cov_ep)
%
%
%                            description:
%
%   This code performs singular value shrinkage on the linearly transformed
%   input matrix ys using the closed formula of Donoho and Gavish by back-
%   projecting, whitening, shrinking, and unnormalizing. The 
%   transformations are assumed to be matrices of the same size. Complex values
%   are allowed for all inputs. The data is *not* assumed to be mean zero;
%   the least-squares estimator of the data is evaluated and subtracted off
%   before shrinkage, then added back. Unlike lintr_fulshr, this code assumes
%   that only a relatively small number of distinct transformations are used
%   relative to the number of data vectors n.
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
%   spec - the spectrum of the expected projection-backprojection
%   xmean - the estimated mean of the signal
%
%

%
%        estimate the mean
%
        [xmean,as2_mean,yback,fracs] = lintr_mean_matr2(ys,as_uniq,ias,nas,m,m2,n);

%
%        center the data and backproject the centered data
%
        ys2 = zeros(m,n);
        for i=1:n
%
        ys(:,i) = ys(:,i) - as_uniq(:,:,ias(i))*xmean;
        ys2(:,i) = as_uniq(:,:,ias(i))'*ys(:,i);
    end

%
%        construct normalization and whitening matrices
%

        cov_ep2 = zeros(m,m,1);
        for i=1:nas
%
        cov_ep2 = cov_ep2 + fracs(i)*as_uniq(:,:,i)'*cov_ep*as_uniq(:,:,i);
    end

        [uep2,sep2] = eig(cov_ep2);
        sep2=diag(sep2);
        wmat = uep2 * diag(1./sqrt(sep2)) * uep2';
        winv = uep2 * diag(sqrt(sep2)) * uep2';
%
        [ua2,sa2] = eig(as2_mean);
        spec=diag(sa2);
        [spec,isort] = sort(spec,'descend');
        ua2=ua2(:,isort);
        as2_inv = ua2 * diag(1./spec) * ua2';

%
%        whiten the data and shrink using Donoho-Gavish
%
        ys3 = wmat * ys2;
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

%
%        now add back the mean
%
        xs_est = xs_est + repmat(xmean,1,n);

        end
%
%
%
%
%
