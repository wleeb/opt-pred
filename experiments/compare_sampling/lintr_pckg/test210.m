        function main()
        startnow();
%
        randn(25,15);
        rand(25,15);
        randn(25,15);
        rand(25,15);
%
        m = 120
%%%        m2 = floor(m/2)

%%%        m2=m/2
        m2=floor(m/1.5)

        gam = .8
        n = floor( m / gam )

        k=2
%
%        set parameters
%
        sig = 2
        ells = sig^2*sqrt(gam) + [1:k];
        ells = sort(ells','descend')

%
%        population covariance
%
        [rcov,u_true] = rand_covc(m,k,ells);

        ells - eigs(rcov,k)
        u_true' * u_true

        imag1=sqrt(-1)

%
%        generate signal and noise matrices
%
        var_ep = (1 + rand(m2,1))/2;
        var_ep = linspace(.5,1,m2)';
        var_ep = sig^2 * var_ep;

        amu = rand(m,1) + imag1*rand(m,1);
        as = randn(m,n) + imag1*randn(m,n) + repmat(amu,1,n);

        [xs,wnoise,zs] = draw_gaussian3(u_true,ells,m,n,k);


        sig_ep = sqrt(var_ep);
        cov_sqrt = randn(m2,m2) + imag1*randn(m2,m2);
        cov_sqrt = cov_sqrt * diag(sig_ep) / sqrt(n);
        cov_ep = cov_sqrt*cov_sqrt';

        alpr=rand()
        betr=sqrt(1-alpr^2)
        alpr^2+betr^2


        epr = randn(m2,n);
        epc = randn(m2,n);
        ep = alpr*epr + betr*imag1*epc;
        ep=cov_sqrt*ep;

        ff=randn(m,m);
        [ff,rr]=qr(ff);
        chk0 = norm(ff*ff' - eye(m))


        nas = 20;

        as_cube = zeros(m2,m,n);
        as_uniq = zeros(m2,m,nas);

        hs = zeros(m2,m,nas);
        for i=1:nas
%
        hh=randn(m,m);
        [hh,rr]=qr(hh);
        hs(:,:,i) = hh(1:m2,:);
        as_uniq(:,:,i) = hs(:,:,i)*diag(as(:,i))*ff';
    end

        size(as_uniq)

        ias = randi(nas,[1,n])


        for i=1:n
%
        as_cube(:,:,i) = as_uniq(:,:,ias(i));
        ys(:,i) = as_cube(:,:,i)*xs(:,i) + ep(:,i);
    end


        [xs_est,whts,errs,spec] = lintr_whit_matr2(ys,as_uniq,ias,nas,m,m2,n,k,cov_ep);
        [xs_est2,whts2,errs2,spec2] = lintr_whit_matr(ys,as_cube,m,m2,n,k,cov_ep);


        sum(errs)
        norm(whts*(xs_est - xs),'fro')^2 / n


        chk0 = norm(xs_est - xs_est2,'fro')
        chk0 = norm(errs-errs2)
        chk0 = norm(whts-whts2)
        chk0 = norm(spec-spec2)


        max(spec)/min(spec)

        stopnow
        end
%
%
%
%
%
        function [ys,xs,ep,as,zs] = gendata_clintr(m,n,k,ells,u_true,var_ep)
%
%        generates complex data from spiked model
%
        imag1 = sqrt(-1)

        [xs,wnoise,zs] = draw_gaussian3(u_true,ells,m,n,k);
%
%        make complex gaussian noise
%
        alpr=rand()
        betr=sqrt(1-alpr^2)
        alpr^2+betr^2

        epr = randn(m,n);
        epc = randn(m,n);
        ep = alpr*epr + betr*imag1*epc;

        sig_ep = sqrt(var_ep);
        ep=repmat(sig_ep,1,n).*ep;

%
%        diagonal reduction matrices
%
        amu = rand(m,1) + imag1*rand(m,1);
        as = randn(m,n) + imag1*randn(m,n) + repmat(amu,1,n);
%
        ys = as.*xs + ep;

        end
%
%
%
%
%
        function [ys,xs,ep,as,zs] = gendata_rlintr(m,n,k,ells,u_true,var_ep)
%
%        generates real-valued data from spiked model
%
        [xs,wnoise,zs] = draw_gaussian3(u_true,ells,m,n,k);
%
%        make complex gaussian noise
%
        ep = randn(m,n);

        sig_ep = sqrt(var_ep);
        ep=repmat(sig_ep,1,n).*ep;

%
%        diagonal reduction matrices
%
        amu = rand(m,1)
        as = randn(m,n) + repmat(amu,1,n);
%
        ys = as.*xs + ep;

        end
%
%
%
%
%
        function [ys,xs,ep,zs] = gendata_cspike(m,n,k,ells,u_true,var_ep)
%
%        generates complex data from spiked model
%
        imag1 = sqrt(-1)

        [xs,wnoise,zs] = draw_gaussian3(u_true,ells,m,n,k);
%
%        make complex gaussian noise
%
        alpr=rand()
        betr=sqrt(1-alpr^2)
        alpr^2+betr^2

        epr = randn(m,n);
        epc = randn(m,n);
        ep = alpr*epr + betr*imag1*epc;

        sig_ep = sqrt(var_ep);
        ep=repmat(sig_ep,1,n).*ep;

        ys = xs + ep;

        end
%
%
%
%
%
        function [xs_opt,u_opt,s_opt,v_opt] = call_optspace(ys,...
           inds,m,n,k,niter,tol)

        if (nargin == 5)
           niter=50
           tol=1d-6
        end

        ys_sparse = sparse(ys);
        [u_opt,s_opt,v_opt,dist] = OptSpace(ys_sparse,k,niter,tol);

        xs_opt = u_opt * s_opt * v_opt';


        end
%
%
%
%
%
        function [xs,wnoise,zs] = draw_gaussian3(u,spec,m,n,k)
%
%       draws a rank k gaussian blob with covariance rcov = u * spec * u',
%       and then adds to it white gaussian noise of variance sig^2
%
        zs=randn(n,k);

        xs = u * diag(sqrt(spec)) * zs';
        wnoise = randn(m,n);
%

        end
%
%
%
%
%
        function [rcov,u] = rand_covc(m,k,ells)
%
%        Generates an m-by-m covariance matrix of rank k with 
%        user-specified spectrum.
%
        u = orthoset_smallc(m,k);
        rcov = u * diag(ells) * u';
        rcov = .5 * (rcov + rcov');

        end
%
%
%
%
%
        function u = orthoset_bigc(m)
%
        imag1 = sqrt(-1);
        g = randn(m,m) + imag1 * randn(m,m);
        [u,r] = qr(g);



        chk0 = norm(u*r - g)
        chk0 = norm(u'*u - eye(m))

        end
%
%
%
%
%
        function u = orthoset_smallc(m,k)
%
%       Produces an m-by-k matrix u with orthonormal columns.
%
        imag1 = sqrt(-1);

        u = zeros(m,k);
        g = randn(m,k) + imag1*randn(m,k);
%
        u(1:m,1) = g(1:m,1) / norm(g(1:m,1));

        for ijk = 2:k
%
        vec = g(1:m,ijk);

        for ii = 1:ijk-1
%
        vec2 = u(1:m,ii);

        prod = sum(vec .* conj(vec2));
        vec = vec - prod*vec2;

        prod = sum(vec .* conj(vec2));
        vec = vec - prod*vec2;

%%%        chk0 = sum(vec .* conj(vec2))

        u(1:m,ijk) = vec/norm(vec);

    end


    end


        end
%
%
%
%
%
        function [rcov,u] = rand_cov(m,k,ells)
%
%        Generates an m-by-m covariance matrix of rank k with 
%        user-specified spectrum.
%
        u = orthoset_small(m,k);
        rcov = u * diag(ells) * u';
        rcov = .5 * (rcov + rcov');

        end
%
%
%
%
%
        function u = orthoset_small(m,k)
%
%       Produces an m-by-k matrix u with orthonormal columns.
%
        u = zeros(m,k);
        g = randn(m,k);
%
        u(1:m,1) = g(1:m,1) / norm(g(1:m,1));

        for ijk = 2:k
%
        vec = g(1:m,ijk);

        for ii = 1:ijk-1
%
        vec2 = u(1:m,ii);

        prod = sum(vec2 .* vec);
        vec = vec - prod*vec2;

        prod = sum(vec2 .* vec);
        vec = vec - prod*vec2;

%%%        chk0 = sum(vec .* vec2)

        u(1:m,ijk) = vec/norm(vec);

    end

    end


        end
%
%
%
%
%
        function inds = rand_inds(m,n,ps)
%
        ff = rand(m,n);
        inds = real(ff <= ps);

        end
%
%
%
%
%
        function startnow()
%
        delete out13
        diary('out13')
        diary on
%
        format short E
%%%        format long E

        randn('state',2)
        randn('state')
        rand('state',1)
        rand('state')

        end
%
%
%
%
%
        function stopnow
%
        diary off
        error('stop')

        end
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This is the end of the test code and the beginning of the shrinkage
%   codes proper.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
        function [xs_in,xs_out,errs_hat,wvals] = lintr_inout(ys_in,ys_out,...
           as_in,as_out,m,nin,nout,k,var_ep)
%
%
%        whiten in-sample data 
%
        as2_in = conj(as_in).*as_in;
        as2_mean = mean(as2_in,2);
        [ys4_in,wvals] = lintr_gen2whi2(ys_in,as_in,m,nin,k,var_ep,as2_mean);

        winv = 1./wvals;
%
%        whiten out-of-sample data
%
        [ys4_out,wvals2] = lintr_gen2whi2(ys_out,as_out,m,nout,...
           k,var_ep,as2_mean);

%
%        apply predictor to both in-sample and out-of-sample, and unwhiten
%
        done=1;
        [xs_in,xs_out,errs_hat] = lintr_inout_spike(ys4_in,ys4_out,m,nin,nout,...
            k,done);

        xs_in = diag(winv) * xs_in;
        xs_out = diag(winv) * xs_out;

        end
%
%
%
%
%
        function [xs_in,xs_out,errs_hat] = lintr_inout_spike(ys_in,ys_out,...
           m,nin,nout,k,sig)
%
%        computes the optimal in-sample (shrinkage) denoiser and 
%        out-of-sample denoiser for data from spiked model with
%        white noise of specified variance
%
        [xs_in,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_white2(ys_in,m,nin,k,sig);

        ells = s_op.^2
        gam_in = m / nin;
        [xs_out,errs_hat] = lintr_out_spike(ys_out,uy,ells,m,nout,k,gam_in,sig);

        chk0 = norm(errs_hat - errs')

        end
%
%
%
%
%
        function [xs_hat,errs_hat] = lintr_out_spike(ys_out,uy_in,ells,...
           m,n_out,k,gam,sig)
%
%        uses PCs from in-sample data to denoise out-of-sample data in
%        standard spiked model with istropic noise
%
        etas = zeros(k,1);
        errs_hat = zeros(k,1);

        gam_out = m / n_out;

        for i=1:k
%
        [etas(i),errs_hat(i)] = lintr_out_eta(ells(i),gam,sig);
    end

        prods = ys_out'*uy_in;
        xs_hat = uy_in*diag(etas)*prods';


        return

%
%        explicit formula, with loops around out-of-sample data and
%        in-sample singular vectors
%
        xs_hat2 = zeros(m,n_out);
        for i=1:n_out

        for j=1:k
%
        p_ij = sum(ys_out(:,i) .* conj(uy_in(:,j)));
        xs_hat2(:,i) = xs_hat2(:,i) + p_ij * etas(j) * uy_in(:,j);
    end
    end

%%%        xs_hat(:,1:2)

        chk0 = norm(xs_hat2 - xs_hat)

%%%        stopnow


        end
%
%
%
%
%
        function [eta,err_hat] = lintr_out_eta(ell,gam,sig)
%
%        out-of-sample denoising coefficient for spike strength ell,
%        for white noise with variance sig^2. Also returns estimated
%        error.
%
        ell = ell / sig^2;
        done=1;
        [rlam,cos_out,cos_inn] = svshr_pop2emp_white2(ell,gam,done)
        eta = ell*cos_out^2 / (ell*cos_out^2 + 1);

        aa = 1 + ell*cos_out^2;
        bb = -2*ell*cos_out^2;
        err_hat = ell + eta^2*aa + eta*bb;

        ell = ell*sig^2;
        err_hat = err_hat * sig^2;

        end
%
%
%
%
%
        function [xs_est,whts,errs,spec] = lintr_whit_matr2(ys,as_uniq,ias,nas,...
            m,m2,n,k,cov_ep)
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
%   as_uniq - the m2-by-m-by-nas cube of reduction matrices (each slice contains the 
%      the reduction matrix), without repetitions. For example, in an image 
%      deblurring problem, the point-spread functions (in image, not Fourier, domain).
%   ias - the n-dimensional vector telling which observation came from which
%      transformation (in the order of matrices in as_uniq)
%   nas - the number of distinct reduction matrices
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

        icounts = zeros(1,nas);
        for i=1:nas
%
        icounts(i) = length(find(ias == i));
    end

        fracs = icounts / n;
        chk0 = sum(fracs) - 1

%
%        backproject the data
%
        ys2 = zeros(m,n);
        for i=1:n
%
        ys2(:,i) = as_uniq(:,:,ias(i))'*ys(:,i);
    end
%
%        construct normalization and whitening matrices
%
        cov_ep2 = zeros(m,m);
        as2_mean = zeros(m,m);

        for i=1:nas
%
        as2_mean = as2_mean + fracs(i)*as_uniq(:,:,i)'*as_uniq(:,:,i);
        cov_ep2 = cov_ep2 + fracs(i)*as_uniq(:,:,i)'*cov_ep*as_uniq(:,:,i);
    end

        [uep2,sep2] = eig(cov_ep2);
        sep2=diag(sep2);
        wmat = uep2 * diag(1./sqrt(sep2)) * uep2';
        winv = uep2 * diag(sqrt(sep2)) * uep2';

        [ua2,sa2] = eig(as2_mean);
        spec=diag(sa2);
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


        end
%
%
%
%
%
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
%        backproject the data
%
        ys2 = zeros(m,n);
        for i=1:n
%
        ys2(:,i) = as(:,:,i)'*ys(:,i);
    end

%
%        construct normalization and whitening matrices
%
        as2 = zeros(m,m,n);
        covs_ep2 = zeros(m,m,n);
        for i=1:n
%
        as2(:,:,i) = as(:,:,i)'*as(:,:,i);
        covs_ep2(:,:,i) = as(:,:,i)'*cov_ep*as(:,:,i);
    end
        as2_mean = mean(as2,3);
        cov_ep2 = mean(covs_ep2,3);

        [uep2,sep2] = eig(cov_ep2);
        sep2=diag(sep2);
        wmat = uep2 * diag(1./sqrt(sep2)) * uep2';
        winv = uep2 * diag(sqrt(sep2)) * uep2';

        [ua2,sa2] = eig(as2_mean);
        spec=diag(sa2);
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


        end
%
%
%
%
%
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
%
%
%
%
        function [xs_est,errs] = lintr_unif(ys,as,m,n,k)
%
%
%                            description:
%
%   This code performs singular value shrinkage on the linearly transformed
%   input matrix ys with missing values, under the assumption that the values
%   are missing uniformly at random. This is the algorithm for missing data
%   described in the OptShrink paper of Raj Rao Nadakuditi.
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
%
%
%
%
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
%
%
%
%
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
        function [ys4,wvals] = lintr_gen2whi(ys,as,m,n,k,var_ep)
%
%        converts general observations A*X + ep to standard spiked model
%        with whitened noise by backprojecting, normalizing, and whitening
%
%        returns the whitened data, and the diagonal whitening transformation
%
%
%         . . . make normalization matrix
%
        as2 = conj(as) .* as;
        as2_mean = mean(as2,2);

%
%        make back-projected data ys2 and normalized backprojected ys3
%
        ys2 = conj(as) .* ys;
        ys3 = ys2  ./ repmat(as2_mean,1,n);

        var_ep2 = as2_mean .* var_ep;
        var_ep3 = var_ep2 ./ as2_mean.^2;
        wvals = 1./sqrt(var_ep3);

%
%        apply transformation to whiten the noise
%
        ys4 = repmat(wvals,1,n) .* ys3;


        end
%
%
%
%
%
        function [ys4,wvals] = lintr_gen2whi2(ys,as,m,n,k,var_ep,as2_mean)
%
%        converts general observations A*X + ep to standard spiked model
%        with whitened noise by backprojecting, normalizing, and whitening
%
%        returns the whitened data, and the diagonal whitening transformation
%
%
%        make back-projected data ys2, normalized backprojected ys3,
%        and effective noise variances for each one
%
        ys2 = conj(as) .* ys;
        ys3 = ys2 ./ repmat(as2_mean,1,n);
%
        var_ep2 = as2_mean .* var_ep;
        var_ep3 = var_ep2 ./ as2_mean.^2;

%
%        now whiten the effective noise
%
        wvals = 1./sqrt(var_ep3);
        ys4 = repmat(wvals,1,n) .* ys3;


        end
%
%
%
%
%
        function [ys2,var_ep2,as2_mean] = lintr_gen2backp(ys,as,m,n,k,var_ep)
%
%        converts general observations A*X + ep to standard spiked model
%        MX + ep2 with colored noise by backprojecting; no normalization
%        or whitening is transformed
%
%        returns the transformed data, and the normalization matrix and 
%        the variances of the new noise term ep2
%
        as2 = conj(as) .* as;
        as2_mean = mean(as2,2);
%
%        make back-projected data ys2
%
        ys2 = conj(as) .* ys;
        var_ep2 = as2_mean .* var_ep;



        end
%
%
%
%
%
        function [ys3,var_ep3] = lintr_gen2std_color(ys,as,m,n,k,var_ep)
%
%        converts general observations A*X + ep to standard spiked model
%        with colored noise by backprojecting and normalizing
%
%        returns the transformed data, and the new noise variances
%
        ys2 = zeros(m,n);
        ys3 = zeros(m,n);

        as2 = conj(as) .* as;
        as2_mean = mean(as2,2);

%
%        make back-projected data ys2, normalized backprojected ys3
%
        ys2 = conj(as).*ys;
        ys3 = ys2 ./ repmat(as2_mean,1,n);
%
        var_ep2 = as2_mean .* var_ep;
        var_ep3 = var_ep2 ./ as2_mean.^2;


        end
%
%
%
%
%
        function [est_mean,ds_mean] = lintr_mean_est(ys,ds,m,n)
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

        return

        est_mean2 = sum(ys,2) ./ sum(inds,2);
        chk0 = norm(est_mean - est_mean2)

        end
%
%
%
%
%


