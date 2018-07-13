        function main()
        startnow();
%
        randn(25,15);
        randn(25,15);
        randn(25,15);
%
        m = 400
        gam = 1.2
        n = floor( m / gam )

        k=3
%
%        set parameters
%
        sig = 3
        ells = 10+sig^2*sqrt(gam) + [1:k].^2;
        ells = sort(ells','descend')
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
        [x_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_white2(ys,m,n,k,sig);

%%%        [x_est,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
%%%           svshr_white(ys,m,n,k);


        s_op
        sqrt(ells)


        for i=1:k
%
        prodv(i)=scapro55(vy(:,i),zs(:,i))/sqrt(n);
        prodv(i)=abs(prodv(i));
        produ(i)=scapro55(uy(:,i),u_true(:,i));
        produ(i)=abs(produ(i));
    end



        [x_est2,s_op2,s_fr2,cos_out2,cos_inn2,uy2,sy2,vy2,errs2] = ...
           svshr_stiel(ys,m,n,k);

        for i=1:k
%
        prodv2(i)=scapro55(vy2(:,i),zs(:,i))/sqrt(n);
        prodv2(i)=abs(prodv2(i));
        produ2(i)=scapro55(uy2(:,i),u_true(:,i));
        produ2(i)=abs(produ2(i));
    end


        prodv2
        cos_inn2
        produ2
        cos_out2



        prodv
        cos_inn
        produ
        cos_out


        s_op2
        sqrt(ells)


        err = norm(xs - x_est,'fro')^2 / n
        err_pred = sum(errs)

        err2 = norm(xs - x_est2,'fro')^2 / n
        err2_pred = sum(errs2)


        [cov_est,ells_hat,ells_fr,cos_out,cos_inn,errs] = svshr_cov_white(...
            ys,m,n,k,sig);
        norm(rcov - cov_est,'fro')^2
        sum(errs)


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
