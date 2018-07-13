        function main()
        startnow();
%
        randn(25,15);
        randn(25,15);
        randn(25,15);
%
        m = 500
        gam = .2
        n = floor( m / gam )

        k=5
%
%        set parameters
%
        sig = 3
        ells = 10+sig^2*sqrt(gam) + [1:k].^2;
        ells = sort(ells','descend')

        ells(k) = sig^2 * sqrt(gam) 

%
%        population covariance
%
        [rcov,u_true] = rand_covc(m,k,ells);

        ells - eigs(rcov,k)
        u_true' * u_true

        imag1 = sqrt(-1)

%
%        generate signal and noise matrices
%
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
        ep=sig*ep;

        sum(abs(ep).^2,2)/n
        sum(abs(ep).^2,1)/m
        sum(abs(ep(:)).^2)/(m*n)

%
%        the observed signal plus noise matrix
%
        ys = xs + ep;

        bedge=sig^2 * (1 + sqrt(gam))^2 


        sy = svd(ys).^2 / n

        sy(k)
        bedge


        [x_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_stiel2(ys,m,n,k,bedge);


        [x_est2,s_op2,s_fr2,cos_out2,cos_inn2,uy2,sy2,vy2,errs2] = ...
           svshr_stiel(ys,m,n,k);


        err2 = norm(xs - x_est2,'fro')^2 / n
        err2_pred = sum(errs2)

        err = norm(xs - x_est,'fro')^2 / n
        err_pred = sum(errs)


        chk0 = norm(x_est-x_est2,'fro')

        svds(x_est,k)

        stopnow
        end
%
%
%
%
%
        function prod = scapro55(x,y,n)
%
        prod = sum(x.*conj(y));
        end
%
%
%
%
%
        function r = guessRank(M_E);
%
%       function from OptSpace_matlab/OptSpace.m to estimate the rank for
%       missing data problem; M_E is sparse
%  
        [n m] = size(M_E);
        epsilon = nnz(M_E)/sqrt(m*n);
        S0 = svds(M_E,100) ;
        S1=S0(1:end-1)-S0(2:end);
        S1_ = S1./mean(S1(end-10:end));
        r1=0;
        lam=0.05;

        while(r1<=0)
%
        for idx=1:length(S1_)
            cost(idx) = lam*max(S1_(idx:end)) + idx;
        end
        [v2 i2] = min(cost);
        r1 = max(i2-1);
        lam=lam+0.05;
    end



%%%        r = r1
%%%        return


        clear cost;
        for idx=1:length(S0)-1
%
        cost(idx) = (S0(idx+1)+sqrt(idx*epsilon)*S0(1)/epsilon  )/S0(idx);
    end
        [v2 i2] = min(cost);
        r2 = max(i2);

        r = max([r1 r2]);
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
        function [ss_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_white2(ys,m,n,k,sig)
%
%
%        This code performs singular value shrinkage on spiked data with 
%        white noise of specified size. It uses the optimal Frobenius shrinker 
%        of Donoho and Gavish; unlike the code svshr_white, this code automatically 
%        rescales for a specified noise variance.
%
%
%                             input parameters:
%
%   ys - the m-by-n matrix of data. The data model is iid columns of dimension m. It
%      is assumed that the columns are NOT divided by sqrt(n).
%   m,n - the dimensions of the data. m is intepreted as the dimension of each 
%      observation, while n is the number of observations.
%   k - the rank of the low-rank signals.
%   sig - the standard deviation of the noise, which is assumed to be white noise.
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
        gam = m/n;

        s_fr = zeros(1,k);
        s_fr88 = zeros(1,k);
        s_op = zeros(1,k);
        cos_out = zeros(1,k);
        cos_inn = zeros(1,k);
%
        [uy,sy,vy] = svshr_svdsmartc(ys,m,n,k);
        sy = sy / sqrt(n) / sig;

        for i=1:k
%
        rlam=sy(i)^2;
        [ell,cos_out(i),cos_inn(i)] = svshr_emp2pop_white(rlam,gam);
        s_op(i) = sqrt(ell);
%
        s_fr(i) = s_op(i) * cos_out(i) * cos_inn(i);
    end

        ells = s_op.^2;
        errs = svshr_errs(ells,cos_inn,cos_out,k);
        s_op = s_op * sig;
        s_fr = s_fr * sig;

        ells = ells*sig^2;
        errs = errs*sig^2;

        s_fr17 = s_fr * sqrt(n);
        ss_est = uy(:,1:k)*diag(s_fr17)*vy(:,1:k)';

        end
%
%
%
%
%
        function [ss_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_white(ys,m,n,k)
%
%
%        This code performs singular value shrinkage on spiked data with 
%        white noise of specified size. It uses the optimal Frobenius shrinker 
%        of Donoho and Gavish; unlike the code svshr_white2, this code assumes
%        the noise variance is 1.
%
%
%                             input parameters:
%
%   ys - the m-by-n matrix of data. The data model is iid columns of dimension m. 
%      It is assumed that the columns are NOT divided by sqrt(n).
%   m,n - the dimensions of the data. m is intepreted as the dimension of each 
%      observation, while n is the number of observations.
%   k - the rank of the low-rank signals.
%   sig - the standard deviation of the noise, which is assumed to be white noise.
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
        gam = m/n;

        s_fr = zeros(1,k);
        s_fr88 = zeros(1,k);
        s_op = zeros(1,k);
        cos_out = zeros(1,k);
        cos_inn = zeros(1,k);
%
        [uy,sy,vy] = svshr_svdsmartc(ys,m,n,k);
        sy = sy / sqrt(n);

        for i=1:k
%
        rlam=sy(i)^2;
        [ell,cos_out(i),cos_inn(i)] = svshr_emp2pop_white(rlam,gam);
        s_op(i) = sqrt(ell);
%
        s_fr(i) = s_op(i) * cos_out(i) * cos_inn(i);
    end
        ells = s_op.^2;
        errs = svshr_errs(ells,cos_inn,cos_out,k);

        s_fr17 = s_fr * sqrt(n);
        ss_est = uy(:,1:k)*diag(s_fr17)*vy(:,1:k)';

        end
%
%
%
%
%
        function [ss_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_stiel2(ys,m,n,k,bedge)
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
%   bedge - the bulk edge of the empirical noise eigenvalues
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
        function [cov_est,ells,ells_fr,cos_out,cos_inn,errs] = ...
            svshr_cov_white(ys,m,n,k,sig)
%
%        This code applies the optimal eigenvalue shrinker for Frobenius loss
%        to estimate the covariance of data from the spiked covariance model
%        with white noise of specified variance.
%
%
%                                input parameters:
%
%   ys - the m-by-n matrix of data. The data model is iid columns of dimension m. 
%      It is assumed that the columns are NOT divided by sqrt(n).
%   m,n - the dimensions of the data. m is intepreted as the dimension of each 
%      observation, while n is the number of observations.
%   k - the rank of the low-rank signals.
%   sig - the standard deviation of the noise, which is assumed to be white noise.
%
%
%                                output parameters:
%
%   cov_est - the m-by-m symmetric estimated covariance matrix. The rank
%      does not exceed k.
%   ells - the k-dimensional vector of estimated population eigenvalues.
%   ells_fr - the k-dimensional vector of optimal eigenvalues for Frobenius loss 
%      (and the eigenvalues of cov_est).
%   cos_out - the k-dimensional vector of estimated cosines between the empirical and 
%      population PCs (i.e. the the outer singular vectors u)
%   cos_inn - the k-dimensional vector of estimated cosines between the empirical and 
%      population inner singular vectors v
%   errs - the k-dimensional vector of estimated errors from each component
%
%
        gam = m/n;

        [uy,sy,vy] = svshr_svdsmartc(ys,m,n,k);
        sy = sy / sig;


        for i=1:k
%
        rlam = sy(i)^2 / n;
        [ells(i),cos_out(i),cos_inn(i)] = svshr_emp2pop_white(rlam,gam);
        ells(i) = ells(i)*sig^2';
        ells_fr(i) = ells(i)*cos_out(i)^2;
%
        errs(i) = (1 - cos_out(i)^4) * ells(i)^2;
    end

        cov_est = uy * diag(ells_fr) * uy';


        end
%
%
%
%
%
        function [rlam,cos_out,cos_inn] = svshr_pop2emp_white2(ell,gam,sig)
%
        ell = ell / sig^2;
        [rlam,cos_out,cos_inn] = svshr_pop2emp_white(ell,gam);
        rlam = sig^2 * rlam;

        end
%
%
%
%
%
        function [rlam,cos_out,cos_inn] = svshr_pop2emp_white(ell,gam)
%
%       takes the population eigenvalue ell in unreduced noise model, and 
%       returns the asymptotic empirical eigenvalue of the unbiased
%       covariance estimator \hat{\Sigma}_s
%

        if (ell <= sqrt(gam))
%
        rlam = (1 + sqrt(gam))^2;
        cos_out = 0;
        cos_inn = 0;
        return
    end


        rlam = (ell + 1) * (1 + gam/ell);

        cos_out = (1 - gam / ell^2) / (1 + gam / ell);
        cos_out = sqrt(cos_out);
%
        cos_inn = (1 - gam / ell^2) / (1 + 1 / ell);
        cos_inn = sqrt(cos_inn);

        end
%
%
%
%
%
        function [ell,cos_out,cos_inn] = svshr_emp2pop_white_flip(rlam,gam)
%
%       takes in empirical eigenvalue of \hat{\Sigma}_{Y/sqrt{delta}}, and 
%       returns the corresponding population eigenvalue if possible (if 
%       rlam exceeds the bulk edge)
%
%       Critically, this code assumes gam > 1
%
        [ell,cos_inn,cos_out] = svshr_emp2pop_white(rlam/gam,1/gam);
        ell = ell*gam;

        end
%
%
%
%
%
        function [ell,cos_out,cos_inn] = svshr_emp2pop_white2(rlam,gam,sig)
%
        rlam = rlam / sig^2;
        [ell,cos_out,cos_inn] = svshr_emp2pop_white(rlam,gam);
        ell = sig^2 * ell;



        end
%
%
%
%
%
        function [ell,cos_out,cos_inn] = svshr_emp2pop_white(rlam,gam)
%
%        computes the population spike and angles, given empirical spike
%        where noise is assumed white. Any values within the bulk are
%        set to zero.
%
        btop = (1+sqrt(gam))^2;
        if (rlam <= btop);
%
        ell=0;
        cos_out=0;
        cos_inn=0;
        return;
    end

        [d_hat,d_der,stra,stra_der,sbar,sbar_der] = ...
           svshr_integrs_white(rlam,gam);
%
        ell = 1/d_hat;
        cos_out = sqrt(stra / (d_der * ell));
        cos_inn = sqrt(sbar / (d_der * ell));


        return
%
%       alternative formulas, explicit in terms of ell and gam
%
        yy = rlam - 1 - gam;
        ell2 = .5*(yy + sqrt(yy^2 - 4*gam));
        chk0 = ell-ell2

        cos_out2 = (1 - gam / ell^2) / (1 + gam / ell);
        cos_out2 = sqrt(cos_out2);
        chk0 = cos_out-cos_out2

        cos_inn2 = (1 - gam / ell^2) / (1 + 1 / ell);
        cos_inn2 = sqrt(cos_inn2);
        chk0 = cos_inn-cos_inn2

        end
%
%
%
%
%
        function dinv = svshr_dinv_white(x,gam)
%
%       computes inverse D-transform at x, for aspect ratio gam,
%       for white noise (unity variance)
%
        dinv = (1 + (1/x)) * (1 + gam*x);
        end
%
%
%
%
%
        function [d_hat,d_der,stra,stra_der,sbar,sbar_der] = ...
           svshr_integrs_white(z,gam)
%
%       evaluates the Stieltjes transform of white noise with unit
%       covariance with aspect ratio gam at z -- only tested for real z
%
        top = (1 - gam- z) + sqrt( (z - 1 - gam)^2 - 4*gam);
        bot = 2*gam*z;
        stra = top/bot;
%
        top_der = -1 + (z-1-gam) / sqrt( (z - 1 - gam)^2 - 4*gam);
        bot_der = 2*gam;
        stra_der = (bot*top_der - top*bot_der) / bot^2;
%
        sbar = gam*stra - (1-gam)/z;
        sbar_der = gam*stra_der + (1-gam)/z^2;
%
        d_hat = stra*sbar*z;
        d_der = stra_der*sbar*z + stra*sbar_der*z + stra*sbar;


        end
%
%
%
%
%
        function [si_hat,cos_out,cos_inn] = svshr_emp2pop_stiel(s,m,n,k,i,bedge)
%
        if (m <= n)
%
        [si_hat,cos_out,cos_inn] = svshr_emp2pop_fat(s,m,n,k,i,bedge);
    end
%
        if (m > n)
%
        [si_hat,cos_inn,cos_out] = svshr_emp2pop_fat(s,n,m,k,i,bedge);
    end

        end
%
%
%
%
%
        function [si_hat,cos_out,cos_inn] = svshr_emp2pop_fat(s,m,n,k,i,bedge)
%
%        estimate i^th singular value using the Stieltjes transform
%
        si = s(i);
        rlam = si^2
%
        if (rlam < bedge)
%
        si_hat = 0;
        cos_out=0;
        cos_inn=0;
        return
    end


        [d_hat,d_der,stra,stra_der,sbar,sbar_der] = svshr_integrs_stiel(s,...
           k,m,n,rlam)

%
%        estimate i^th singular value
%
        si_hat = 1/sqrt(d_hat);
%
%
%        estimate i^th inner and outer angles between pop and emp singular
%        vectors
%
        cos_inn = sqrt(sbar / (d_der * si_hat^2))
        cos_out = sqrt(stra / (d_der * si_hat^2))


        end
%
%
%
%
%
        function [d_hat,d_der,stra,stra_der,sbar,sbar_der] = ...
           svshr_integrs_stiel(s,k,m,n,rlam)
%
%        computes empirical estimates of the stieljes and D-transforms 
%        at the value rlam, with empirical singular values s (already
%        normalized by dimension). Critically, code assumes that m < n
%
        stra = 0
        stra_der = 0

        for i=k+1:m
%
        stra = stra + 1 / (s(i)^2 - rlam);
        stra_der = stra_der + 1 / (rlam - s(i)^2)^2;        
    end

        stra = stra/(m-k);
        stra_der = stra_der/(m-k);

        sbar = (m-k)*stra/(n-k) - (n-m)/rlam/(n-k);
        sbar_der = (m-k)*stra_der/(n-k) + (n-m)/rlam^2/(n-k);

        d_hat = stra*sbar*rlam;
        d_der = stra_der*sbar*rlam + stra*sbar_der*rlam + stra*sbar;


        end
%
%
%
%
%
        function svshr_compare_mp(evals,gam,sig,ifig)
%
%        plots histogram of evals and overlays plot of continuous
%        Marchenko-Pastur density, white noise with variance 1,
%        with aspect ratio gam. ifig is the figure number of plot
%

        x0 = sig^2*(1-sqrt(gam))^2;
        x1 = sig^2*(1+sqrt(gam))^2;

        nvals=100;
        ts=linspace(x0,x1,nvals);

        for i=1:nvals
%
        vals_mp(i) = svshr_mp_eval(ts(i),gam,sig);
    end

        figure(ifig);
        subplot(1,2,1);
        hh = histogram(evals,100,'Normalization','pdf');
        hold on; plot(ts,vals_mp,'linewidth',3)


        m=length(evals);

        subplot(1,2,2);
        plot(evals,'*','Color','b')
        hold on; plot(0:m+1,x1*ones(m+2,1),'LineWidth',2,'Color','r')
        hold on; plot(0:m+1,x0*ones(m+2,1),'LineWidth',2,'Color','r')

        xlim([0,m+1])


        set(figure(ifig),'Position',[500,500,1300,500])

        end
%
%
%
%
%
        function val = svshr_mp_eval(t,gam,sig)
%
%        evaluates standard Marchenko-Pastur density with 
%        aspect ratio gam at t
%
        x0 = sig^2*(1-sqrt(gam))^2;
        x1 = sig^2*(1+sqrt(gam))^2;
        val = (x1 - t) * (t - x0);
        val = sqrt(val) / (2*pi*t*gam) / sig^2;

        end
%
%
%
%
%
        function w = svshr_joint_basis(u,v,m,k)
%
%
%       make the joint basis w
%
        g = zeros(m,m);
        g(1:m, 1:k) = u(1:m, 1:k);
        g(1:m, k+1:2*k) = v(1:m, 1:k);
        g(1:m, 2*k+1:m) = v(1:m, k+1:m-k);

        w = svshr_gs_dumb(g,m);

        chk0 = norm(w(1:m,1:k) - u(1:m, 1:k))

        iperm = 1:m;
        for j=1:k
%
        iperm(2*j-1) = j;
        iperm(2*j) = k+j;
    end


%%%        iperm
        w = w(1:m,iperm);


        end
%
%
%
%
%
        function u = svshr_gs_dumb(g,m)
%
        u = zeros(m,m);
%
        u(1:m,1) = g(1:m,1) / norm(g(1:m,1));

        for ijk = 2:m
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
%
        end
%
%
%
%
%
        function [est_mean,ds_mean] = svshr_mean_est(ys,ds,m,n)
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
        function s_hat = svshr_emp2frob(s,m,n,k,i)
%
%        estimates i^th singular value using RRN's formula
%

        if (m > n)
%
        m2=m;
        m=n;
        n=m2;
    end

        s_tail = s(k+1:m);
        si = s(i);
%
        n1=n-k;
        m1=m-k;

        t1 = sum(1./(si^2 - s_tail.^2))*si/n1 + (n1-m1)/si/n1;
        t2 = sum(1./(si^2 - s_tail.^2))*si/m1;
%
        t3 = sum(1./(si^2 - s_tail.^2) - 2*si^2./(si^2 - s_tail.^2).^2)/m1;
        t4 = sum(1./(si^2 - s_tail.^2) - 2*si^2./(si^2 - s_tail.^2).^2)/n1;
        t4 = t4 - 2*(n1-m1)/si^2/n1 + (n1-m1)/si^2/n1;
%
        dbot = t1*t3 + t2*t4;
        dtop = t1*t2;
%
        s_hat = -2*dtop/dbot;



        end
%
%
%
%
%
        function [err_fr,cos_out,cos_inn] = svshr_fro_err(ells,...
           u,u_hat,v,v_hat,k)
%
        cos_out = zeros(1,k);
        cos_inn = zeros(1,k);

        for i=1:k
%
        cos_out(i) = sum(u(:,i) .* u_hat(:,i));
        cos_inn(i) = sum(v(:,i) .* v_hat(:,i));
    end

        err_fr = svshr_err_fmla(ells,cos_inn,cos_out,k);
        end
%
%
%
%
%
        function err_fr = svshr_err_fmla(ells,cos_inn,...
           cos_out,k)
%
        err_fr = 0;
        for i=1:k
%
        err_fr = err_fr + ells(i) * (1 - cos_inn(i)^2 * cos_out(i)^2);
    end
        end
%
%
%
%
%
        function errs = svshr_errs(ells,cos_inn,...
           cos_out,k)
%
        errs = zeros(1,k);
        for i=1:k
%
        errs(i) = ells(i) * (1 - cos_inn(i)^2 * cos_out(i)^2);
    end
        end
%
%
%
%
%
        function [u,s,v] = svshr_svdsmart(a,m,n,k)
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
        function [u,s,v] = svshr_svdsmartc(a,m,n,k)
%
        if (k==1)
%
        [u,s,v] = svds(a,2);
        s=s(1);
        u=u(:,1);
        v=v(:,1);
        return;
    end

        if (m/k > 20)
%
        [u,s,v] = svds(a,k);
        s=diag(s);
        return;
    end

        if (m/k <= 20)
%
        [u,s,v] = svd(a,'econ');
        s=diag(s(1:k,1:k));
%%%        s=diag(s);
        u = u(1:m,1:k);
        v = v(1:n,1:k);
        return;
    end

        end
%
%
%
%
%
