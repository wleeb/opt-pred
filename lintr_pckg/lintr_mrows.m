        function [uhat,xhat,errhat,whts] = lintr_mrows(ys,as,m,n,k,sig,isig)
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
%   Importantly: unlike other codes in this package, it is assumed that the
%   transformations are applied to the *rows* of the matrix, not the columns.
%
%
%                           input parameters:
%
%   ys - the m-by-n data matrix of reduced observations.
%   as - the m-by-n matrix of reduction matrices (each column contains the 
%      diagonal of the reduction matrix).
%   m,n - the dimensionality and number of observations, respectively
%   k - the rank of the signal.
%   sig - the standard deviation of the noise
%   isig - set to 1 if noise level sig should be estimated automatically
%
%
%                        output parameters:
%
%   uhat - the estimated left singular vectors. Size m-by-k.
%   xhat - the estimated matrix. Size m-by-n, and rank (at most) k.
%   whts - m-by-1 vector defining the diagonal matrix with weights defining
%      the loss function
%   errhat - estimate of the weighted error || W*(xhat - xs) ||^2 (squared
%      Frobenius norm) where W is diag(whts)
%
%

%
%        transpose the data
%
        ys = ys';
        as = as';
        m2=m;
        m=n;
        n=m2;

%
%        center the data
%
        [xmean,as2_mean,yback] = lintr_mean_diag(ys,as,m,n);
        ys2 = ys - as.*repmat(xmean,1,n);

%
%        if desired, estimate noise level
%
        if (isig == 1)
%
        sig = sum(ys2(:).^2) / (sum(as(:).^2));
        sig = sqrt(sig)
    end

%
%        backproject and whiten
%
        [ys3,wvals2,as2_mean] = lintr_gen2whi3(ys2,as,m,n,k,sig);
        wvals = wvals2 ./ as2_mean;
        winv = 1./wvals;

%
%        shrink using Donoho-Gavish
%
        sig0=1;
        [xhat,s_op,s_fr,cos_out,cos_inn,vy,sy,uhat,errs] = ...
           svshr_white2(ys3,m,n,k,sig0);

        chk0 = norm(uhat*diag(s_fr)*vy' - xhat'/sqrt(n),'fro')

%
%        unnormalize the estimate
%
        whts_inv = (1./as2_mean) .* winv;
        xhat = repmat(whts_inv,1,n) .* xhat;
        whts = 1./whts_inv / sqrt(n);

        xhat = xhat + repmat(xmean,1,n);
        errhat = sum(errs);

%
%        finally, transpose back
%
        xhat = xhat';


        end
%
