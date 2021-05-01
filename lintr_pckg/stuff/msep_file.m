        function [ss_hat,etas,err_hat] = dumm(ys_out,m,n_out,k,gam,...
           dmu,ys,n)
%
%                          description:
%
%   This code takes in two diagonally reduced datasets from the same 
%   distribution, and uses the inner dataset to optimally denoise the outer
%   dataset. The diagonal entries of the reduction matrices are assumed to
%   be iid, and the rank of the model is assumed to be known. The noise is
%   NOT reduced in this model ('unre' means 'unreduced noise'.)
%
%                          input parameters:
%
%   ys_out - the m-by-n_out matrix of outer datapoints.
%
%   m - the dimension of the model
%
%   n_out - the number of outer samples
%
%   k - the rank of the model
%
%   dmu,demm - the first and second moments of the diagonal elements of the
%      reduction matrices.
%
%   ys - the m-by-n matrix of inner points (used to estimate the principal
%      components and singular values of the model)
%
%   n - the number of inner data points
%
%                          output parameters:
%
%   ss_hat - the m-by-n_out estimated matrix.
%
%   etas - the k-dimensional vector of shrinkage coefficients.
%
%   err_hat - the estimated error (in squared Frobenius norm).
%

        ells = zeros(1,k);
        cosl = zeros(1,k);
%
        [uy,sy,vy] = unre_svdsmart(ys,m,n,k);
%
        for i=1:k
%
        rlam_y = sy(i)^2 / n;
        [ells(i),cosl(i),cos_inn,ierr] = unre_yemp2pop(rlam_y,gam,dmu)
    end

        [ss_hat,etas,err_hat] = unre_eblp_cheat(ys_out,m,n_out,k,...
           gam,dmu,ells,cosl,uy);

        end
%
%
%
%
%



        function [ss_hat,etas,err_hat] = stoopid(ys_out,m,n_out,k,gam,...
           dmu,ys,n)
%
%                          description:
%
%   This code takes in two diagonally reduced datasets from the same 
%   distribution, and uses the inner dataset to optimally denoise the outer
%   dataset. The diagonal entries of the reduction matrices are assumed to
%   be iid, and the rank of the model is assumed to be known. The noise is
%   NOT reduced in this model ('unre' means 'unreduced noise'.)
%
%                          input parameters:
%
%   ys_out - the m-by-n_out matrix of outer datapoints.
%
%   m - the dimension of the model
%
%   n_out - the number of outer samples
%
%   k - the rank of the model
%
%   dmu,demm - the first and second moments of the diagonal elements of the
%      reduction matrices.
%
%   ys - the m-by-n matrix of inner points (used to estimate the principal
%      components and singular values of the model)
%
%   n - the number of inner data points
%
%                          output parameters:
%
%   ss_hat - the m-by-n_out estimated matrix.
%
%   etas - the k-dimensional vector of shrinkage coefficients.
%
%   err_hat - the estimated error (in squared Frobenius norm).
%

        ells = zeros(1,k);
        cosl = zeros(1,k);
%
        [uy,sy,vy] = unre_svdsmart(ys,m,n,k);
%
        for i=1:k
%
        rlam_y = sy(i)^2 / n;
        [ells(i),cosl(i),cos_inn,ierr] = unre_yemp2pop(rlam_y,gam,dmu)
    end

        [ss_hat,etas,err_hat] = unre_eblp_cheat(ys_out,m,n_out,k,...
           gam,dmu,ells,cosl,uy);

        end
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%             end of debugging code, beginning of real code below
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
        function [ss_hat,etas,err_hat] = unre_eblp(ys_out,m,n_out,k,gam,...
           dmu,ys,n)
%
%                          description:
%
%   This code takes in two diagonally reduced datasets from the same 
%   distribution, and uses the inner dataset to optimally denoise the outer
%   dataset. The diagonal entries of the reduction matrices are assumed to
%   be iid, and the rank of the model is assumed to be known. The noise is
%   NOT reduced in this model ('unre' means 'unreduced noise'.)
%
%                          input parameters:
%
%   ys_out - the m-by-n_out matrix of outer datapoints.
%
%   m - the dimension of the model
%
%   n_out - the number of outer samples
%
%   k - the rank of the model
%
%   dmu,demm - the first and second moments of the diagonal elements of the
%      reduction matrices.
%
%   ys - the m-by-n matrix of inner points (used to estimate the principal
%      components and singular values of the model)
%
%   n - the number of inner data points
%
%                          output parameters:
%
%   ss_hat - the m-by-n_out estimated matrix.
%
%   etas - the k-dimensional vector of shrinkage coefficients.
%
%   err_hat - the estimated error (in squared Frobenius norm).
%

        ells = zeros(1,k);
        cosl = zeros(1,k);
%
        [uy,sy,vy] = unre_svdsmart(ys,m,n,k);
%
        for i=1:k
%
        rlam_y = sy(i)^2 / n;
        [ells(i),cosl(i),cos_inn,ierr] = unre_yemp2pop(rlam_y,gam,dmu)
    end

        [ss_hat,etas,err_hat] = unre_eblp_cheat(ys_out,m,n_out,k,...
           gam,dmu,ells,cosl,uy);

        end
%
%
%
%
%
        function [ss_hat,etas,err_hat] = unre_eblp_cheat(ys_out,m,n_out,...
           k,gam,dmu,ells,cosl,uy);
%
%                          description:
%
%   This code takes in one diagonally reduced dataset, and the estimated
%   principal components and eigenvalues obtained from a difference
%   dataset drawn from the same distribution, and uses the inner dataset 
%   to optimally denoise the outer dataset. The diagonal entries of the 
%   reduction matrices are assumed to be iid, and the rank of the model is 
%   assumed to be known.
%
%                          input parameters:
%
%   ys_out - the m-by-n_out matrix of outer datapoints.
%
%   m - the dimension of the model
%
%   n_out - the number of outer samples
%
%   k - the rank of the model
%
%   dmu,demm - the first and second moments of the diagonal elements of the
%      reduction matrices.
%
%   ys - the m-by-n matrix of inner points (used to estimate the principal
%      components and singular values of the model)
%
%   n - the number of inner data points
%
%                          output parameters:
%
%   ss_hat - the m-by-n_out estimated matrix.
%
%   etas - the k-dimensional vector of shrinkage coefficients.
%
%   err_hat - the estimated error (in squared Frobenius norm).
%
%
        ss_hat = zeros(m,n_out);
        err_hats = zeros(1,k);
        etas = zeros(1,k);

%
%       find optimal shrinkage coefficients, and expected error
%
        for i=1:k
%
        [etas(i),err_hats(i)] = unre_eta(dmu,ells(i),cosl(i),gam);
    end

        err_hat = sum(err_hats);

%
%       denoise each new vector
%
        for i_out = 1:n_out
%
        ss_hat(1:m,i_out) = unre_eblp_apply(m,k,etas,uy,ys_out(:,i_out));

    end

        end
%
%
%
%
%
        function s_hat = unre_eblp_apply(m,k,etas,uy,y_out)
%
        s_hat = zeros(m,1);

        for i=1:k
%
        prod_i = sum(uy(1:m,i) .* y_out);
        coeff_i = etas(i)*prod_i;
%
        s_hat = s_hat + coeff_i*uy(1:m,i);
    end

        end
%
%
%
%
%
        function [eta,err_hat] = unre_eta(dmu,ell,cosl,gam)
%
        eta = dmu*ell*cosl^2 / (dmu^2*ell*cosl^2 + 1);
%
        aa = 1 + dmu^2*ell*cosl^2;
        bb = -2*dmu*ell*cosl^2;
        err_hat = ell + eta^2*aa + eta*bb;


        end
%
%
%
%
%
        function [s_est,err_hat,uy,sy,vy,sv_hats,smu] = unre_svshrink(...
           ys,ds,dmu,m,n,k,sig_hat,if_cen)
%
%
%                           description:
%
%   This code performs optimal singular value shrinkage on a matrix with
%   diagonally reduced entries, where the reduction entries are iid, and
%   the noise is not reduced ('unre' means unreduced noise).
%
%
%                         input parameters:
%
%   ys - the m-by-n matrix with zeros in the place of missing values
%
%   inds - the m-by-n matrix with a 1 for every observed value, 0 outside
%
%   m,n - the dimensions of the matrices ys and inds
%
%   k - the rank (estimated or true) of the clean matrix
%
%
%                        output parameters
%
%   s_est - the shrunken matrix
%
%   err_hat - the estimated error
%
%   uy - the m-by-k matrix of left singular vectors of ys
%
%   sy - the top k singular values of the matrix ys / sqrt(n)
%
%   vy - the n-by-k matrix of right singular vectors of ys
%
%   sv_hats - the shrunken singular values (top k s.v.'s of s_est)
%
%   if_cen - specifies if the mean is to be subtracted off before shrinkage
%      1 - estimate and subtract off mean, shrink, then add back mean
%      0 - do not estimate or adjust for mean
%
%
%


        ells = zeros(1,k);
%
        smu=0;
        if (if_cen == 1)
%
        [smu,ds_mean] = unre_mean_estim(ys,ds,m,n);
        for j=1:n
%
        ys(:,j) = ys(:,j) - ds(:,j) .* smu;
    end

    end


        if (sig_hat < 0)
           ifsvd = 1;
           sig_hat = sighat(ys,m,n,ifsvd);
        end

        ys = ys / sig_hat;
%
        gam = m/n;
        [uy,sy,vy] = unre_svdsmart(ys,m,n,k);
        sy = sy / sqrt(n);

        for i=1:k
%
        rlam_y = sy(i)^2;
        [ells(i),cos_out,cos_inn,ierr] = unre_yemp2pop(rlam_y,gam,dmu);

    end

        [s_est,err_hat,errs,sv_hats] = unre_svshrink_cheat(ys,m,n,k,dmu,...
           ells,uy,vy);
%
%       rescale data to have original noise level
%
        s_est = s_est * sig_hat;
        err_hat = err_hat * sig_hat^2;
        sv_hats = sv_hats * sig_hat;
        sy = sy * sig_hat;
%
%   add back mean, if necessary
%
        if (if_cen == 1)
%
        for j=1:n
%
        s_est(:,j) = s_est(:,j) + smu;
    end

    end
        end
%
%
%
%
%
        function [est_mean,ds_mean] = unre_mean_estim(ys,ds,m,n)
%
        ds_mean = mean(ds,2);
        est_mean = mean(ys,2) ./ ds_mean;


        return
%
%   the slow version:
%
        est_mean = zeros(m,1);
        ds_mean = zeros(m,1);
%
        for i = 1:m
%
        ds_sum=0;
        for j=1:n
%
        est_mean(i) = est_mean(i) + ys(i,j);
        ds_sum = ds_sum + ds(i,j);
    end

        est_mean(i) = est_mean(i)/ds_sum;
        ds_mean(i) = ds_sum / n;
    end


        end
%
%
%
%
%
        function [s_est,err_hat,errs,sv_hats] = unre_svshrink_cheat(ys,...
           m,n,k,dmu,ells,uy,vy)
%
%
%                         description:
%
%   This code performs asymptotically optimal singular value shrinkage on
%   the input matrix ys that has been diagonally reduced, with unreduced
%   noise. The eigenvalues ells of the population covariance and the
%   empirical singular vectors are provided by the user.
%
%                       input parameters:
%
%   ys - the m-by-n matrix of reduced observations
%
%   m,n - the dimensions of the input matrix (m = dimension of each vector,
%      n = the number of vectors)
%
%   k - the rank (or estimated rank) of the data
%
%   dmu - the mean of the diagonal entries of the reduction matrices
%
%   ells - the k-dimensional vector of eigenvalues of covariance 
%
%   uy - the m-by-k matrix of left empirical singular vectors of ys
%
%   vy - the n-by-k matrix of right empirical singular vectors of ys
%
%
%                        output parameters:
%
%   s_est - the estimated matrix
%
%   err_hat - the estimated error
%
%   errs - the individual errors for each singular value
%
%   sv_hats - the shrunken singular values
%
        sv_hats = zeros(k,1);
        errs = zeros(k,1);
        s_est = zeros(m,n);
        err_hat = 0;

        gam = m/n;

        for i=1:k
%
        ell = ells(i);
        [sv_hat,err_hat,cos_out,cos_inn] = unre_svshrink_one(ell,gam,dmu);
%
        sv_hat = sv_hat * sqrt(n);
        sv_hats(i) = sv_hat;
%
        s_est = s_est + sv_hat*uy(:,i) * vy(:,i)';
        errs(i) = err_hat;

    end

        err_hat = sum(errs);

        end
%
%
%
%
%
        function [sv_hat,err_hat,cos_out,cos_inn] = unre_svshrink_one(...
           ell,gam,dmu)
%
%       given eigenvalue ell of population covariance, returns the
%       shrunken singular value (for Frobenius loss), angles, and 
%       predicted error
%
        [rlam_s,rlam_y,cos_out,cos_inn,edge_t,edge_b,ierr] = ...
           unre_pop2emp(ell,gam,dmu);

        sv_hat = sqrt(ell) * cos_out * cos_inn;
        err_hat = (1 - cos_inn^2 * cos_out^2) * ell;

        end
%
%
%
%
%
        function [u_s,cov_op,cov_fr,err_op,err_fr,ells,cos_out,cos_inn,...
           ierrs] = unre_cov(ys,m,n,ds,k,if_cen)
        ells = zeros(1,k); cos_out = zeros(1,k); cos_inn = zeros(1,k);
        errs_fr = zeros(1,k); errs_op = zeros(1,k); ierrs = zeros(1,k);
%
%
%                         description:
%
%   This code returns an estimate of the covariance matrix of a random 
%   vector whose observations have been diagonally reduced by isotropic
%   reduction matrices. White noise is added AFTER reduction; that is, the
%   noise is UNREDUCED. An unbiased estimator of the covariance is 
%   constructed, and the eigenvalues of this matrix are shrunk in an 
%   optimal way. The rank k of the model us assumed known.
%
%                         input parameters:
%
%   ys - the m-by-n matrix of reduced observations.
%
%   ds - the m-by-n matrix whose columns are the diagonals of the reduction
%      matrices.
%
%   m,n - the dimension (number of rows) and number of samples (number of
%      columns), respectively.
%
%   k - the rank of the model (maximum rank of estimated covariance).
%
%   if_cen - indicates if the data should be centered or not
%      0 - do not center the data
%      1 - estimate the mean and center the data
%
%                       output parameters:
%
%   cov_op - the m-by-m optimally shrunken covariance matrix for operator 
%      norm loss
%
%   cov_fr - the m-by-m optimally shrunken covariance matrix for Frobenius 
%      norm loss
%
%   ells_op - the k-dimensional vector of eigenvalues of cov_op
%
%   ells_fr - the k-dimensional vector of eigenvalues of cov_fr
%
%   u_s - the m-by-k matrix of eigenvectors of ells_op and ells_fr (i.e.
%      the estimated principals components of the data)
%
%   err_op - the estimated error in operator norm between cov_op and the
%      population covariance matrix
%
%   err_fr - the estimated error in squared Frobenius norm between cov_fr 
%      and the population covariance matrix
%
%   sig_hat - the estimated standard deviation of the noise
%
%   ierrs - the k-dimensional vector that indicates which estimated
%      eigenvalues did not exceed the bulk edge (and hence were set to 0).
%      ierrs(j) = 0 if the j^th eigenvalue exceeds the bulk edge, and is 1
%      otherwise.
%
%
        gam = m/n;

        dmu = mean(ds(:));
        demm = mean(ds(:).^2);
        dvar = demm - dmu^2;
%

        if (if_cen == 1)
%
        [est_mean,ds_mean] = unre_mean_estim(ys,ds,m,n);
        for j=1:n
%
        ys(:,j) = ys(:,j) - ds(:,j) .* est_mean;
    end

    end


        ifsvd = 1;
        sig_hat = sighat(ys,m,n,ifsvd)




        sig_hat = 1



        ys = ys / sig_hat;

        
        [cov_s,cov_y] = unre_unbiased(ys,m,n,dmu,demm,dvar);
        [u_s,rlam_s] = unre_eigsmart(cov_s,m,k);

%
%       find the population eigenvalues and cosines
%
        for i=1:k
%
        [ells(i),cos_out(i),cos_inn(i),ierrs(i)] = unre_semp2pop(...
           rlam_s(i),gam,dmu)
%
%    Frobenius and operator norm optimal singular values
%
        [ells_op(i),ells_fr(i),errs_fr(i),errs_op(i)] = unre_evshrink(...
           ells(i),gam,dmu);

        ells_op(i) = ells_op(i)*sig_hat^2;
        ells_fr(i) = ells_fr(i)*sig_hat^2;

    end

        err_fr = sum(errs_fr(1:k)) * sig_hat^4;
        err_op = max(errs_op(1:k)) * sig_hat^2;

%
%   construct optimal covariance matrices
%
        cov_fr = u_s * diag(ells_fr) * u_s';
        cov_op = u_s * diag(ells_op) * u_s';



        end
%
%
%
%
%
        function [ell_op,ell_fr,err_fr,err_op] = unre_evshrink(ell,gam,dmu)
%
        [rlam_s,rlam_y,cos_out,cos_inn,edge_t,edge_b,ierr] = ...
           unre_pop2emp(ell,gam,dmu);

        ell_op = ell;
        ell_fr = ell*cos_out^2;
        err_fr = (1 - cos_out^4) * ell^2;

        sin_out = sqrt(1 - cos_out^2);
        err_op = ell_op * sin_out;

        end
%
%
%
%
%
        function [ell,cos_out,cos_inn,ierr] = unre_yemp2pop(rlam_y,gam,dmu)
%
%       takes in empirical eigenvalue of \hat{\Sigma}_s, and returns the 
%       corresponding population eigenvalue ir possible (if rlam exceeds
%       the bulk edge); otherwise, returns ierr = 1
%
        cos_out=0
        cos_inn=0
%
        bulk_edge = (1 + sqrt(gam))^2
        if (rlam_y <= bulk_edge + eps*100)
%
        ell = 0
        ierr = 1
        return
    end

        ierr = 0
        ell2 = unre_compute_ell(rlam_y,gam);
        ell = ell2/dmu^2

        cos_out = (1 - gam / ell2^2) / (1 + gam / ell2)
        cos_out = sqrt(cos_out)
%
        cos_inn = (1 - gam / ell2^2) / (1 + 1 / ell2)
        cos_inn = sqrt(cos_inn)


        end
%
%
%
%
%
        function [ell,cos_out,cos_inn,ierr,rlam_y] = unre_semp2pop(...
           rlam_s,gam,dmu)
%
%       takes in empirical eigenvalue of \hat{\Sigma}_s, and returns the 
%       corresponding population eigenvalue ir possible (if rlam exceeds
%       the bulk edge); otherwise, returns ierr = 1
%
        cos_out=0
        cos_inn=0

        rlam_y = dmu^2 * rlam_s + 1
%
        bulk_edge = (1 + sqrt(gam))^2
        if (rlam_y <= bulk_edge + eps*100)
           ell = 0
           ierr = 1
           return
        end

        ierr = 0
        ell2 = unre_compute_ell(rlam_y,gam);
        ell = ell2/dmu^2

        cos_out = (1 - gam / ell2^2) / (1 + gam / ell2)
        cos_out = sqrt(cos_out)
%
        cos_inn = (1 - gam / ell2^2) / (1 + 1 / ell2)
        cos_inn = sqrt(cos_inn)


        end
%
%
%
%
%
        function [rlam_s,rlam_y,cos_out,cos_inn,edge_t,edge_b,ierr] = ...
           unre_pop2emp(ell,gam,dmu)
%
%       takes the population eigenvalue ell in unreduced noise model, and 
%       returns the asymptotic empirical eigenvalue of the unbiased 
%       covariance estimator \hat{\Sigma}_s
%
        cos_out=0
        cos_inn=0
        ell2 = dmu^2 * ell
        aa = 1/dmu^2
%
        edge_t = (1 + sqrt(gam))^2  / dmu^2 - aa
        edge_b = (1 - sqrt(gam))^2  / dmu^2 - aa
        rlam_y = (ell2 + 1) * (1 + gam/ell2)
%
        thresh = sqrt(gam);
        ierr = 0;
%
        if (ell2 <= thresh)
           rlam_y = (1 + sqrt(gam))^2;
           rlam_s = edge_t;
           ierr = 1;
           return
        end
%
        rlam_s = rlam_y / dmu^2 - aa
%

        cos_out = (1 - gam / ell2^2) / (1 + gam / ell2)
        cos_out = sqrt(cos_out)
%
        cos_inn = (1 - gam / ell2^2) / (1 + 1 / ell2)
        cos_inn = sqrt(cos_inn)

        end
%
%
%
%
%
        function [scov,ycov] = unre_unbiased(ys,m,n,dmu,demm,dvar)
%
%       returns unbiased estimator of \Sigma_s from the raw data ys
%
        ycov = ys * ys' / n;
        di = diag(diag(ycov));

        scov = ycov / dmu^2 - dvar*di/(demm*dmu^2) - eye(m)/demm;



        end
%
%
%
%
%
        function [u,s,v] = unre_svdsmart(a,m,n,k)
%
        if (m/k > 20)
%
        [u,s,v] = svds(a,k);
        s=diag(s);
    end

        if (m/k <= 20)
%
        [u,s,v] = svd(a,'econ');
        s=diag(s(1:k,1:k));
%%%        s=diag(s);
        u = u(1:m,1:k);
        v = v(1:n,1:k);
    end

        end
%
%
%
%
%
        function ell=unre_compute_ell(rlam,gam)
%
        yy = rlam - 1 - gam;
        ell = .5*(yy + sqrt(yy^2 - 4*gam));


        end
%
%
%
%
%
        function [u,s]=unre_eigsmart(a,m,k)
%
%       returns the top k eigenvectors/values of the m-by-m symmetric 
%       matrix a. Note that s is the vector of eigenvalues, not a diagonal 
%       matrix.
%
        if (m/k > 10)
%
        [u,s] = eigs(a,k);
        s = sort(diag(s),'descend');
    end

        if (m/k <= 10)
%
        [u,s] = eig(a);
        [s,ii] = sort(diag(s),'descend');
        s=s(1:k);
        u = u(1:m,ii(1:k));
    end


        end


