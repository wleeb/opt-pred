        function main()
        startnow();
%
        randn(25,12);
        rand(25,12);
%
        m = 500
        gam = .8
        n = floor( m / gam )

        k=4
%
%        set parameters
%
        sig = 1

        delta=.5
        ells = sig^2*sqrt(gam)/delta + [1:k] + 15;
        ells = sort(ells','descend')

%
%        population covariance
%
        [rcov,u_true] = rand_cov(m,k,ells);
%%%        chk0 = norm(u_true*diag(ells)*u_true' - rcov,'fro')

%
%        generates complex data from spiked model
%
        var_ep = sig^2 * ones(m,1);


        var_ep = 1 + rand(m,1).^2 + linspace(1,2,m)';
        ep = diag(sqrt(var_ep)) * randn(m,n);

        zs = randn(n,k);
        xs = u_true * diag(sqrt(ells)) * zs';

        pvec = linspace(.2,.8,m)';
        ps = repmat(pvec,1,n);
%%%        ps = delta*ones(m,n);
        inds = rand_inds(m,n,ps);
        ys = inds.*(xs + ep);

        sy = svd(ys);
%%%        plot(sy,'*')

        bedge = sig^2*delta*(1+sqrt(gam))^2 + .01
        norm(inds.*ep)^2 / n

        sy(k).^2 / n


        as = inds;
        nsims=10;
        bedge2 = lintr_edge_diag(pvec,var_ep,m,n,nsims);

        bedge
        bedge2
        norm(inds.*ep)^2/n

        [xs_est,errs] = lintr_unif2(ys,inds,m,n,k,bedge2);


        norm(xs_est - xs,'fro')^2/n
        sum(errs)

        svds(xs_est,k)


        stopnow
        end
%
%
%
%
%
        function xhat = wifi_diag(y,a,u_true,ells,m,k,sig)
%
        c = repmat(a,1,k).*(u_true*diag(sqrt(ells)));
        [uc,sc,vc] = svshr_svdsmartc(c,m,k,k);
        sc=sc.^2;
%
%        project y onto the column space of C
%
        yc = uc*(uc'*y);
        yc2 = y - yc;
%%%        chk0 = norm(y)^2 - (norm(yc)^2 + norm(yc2)^2)
%
%        apply the Wiener filter
%
        ucy = uc'*y;
        ducy = ucy ./ (sig^2 + sc);
        binv_y = uc*ducy + yc2/sig^2;
%
        xhat = u_true*diag(ells)*(u_true'*(conj(a).*binv_y));



        end
%
%
%
%
%
        function [xcov,ells,wvals,err_fr] = covest(ys,as,m,n,k,var_ep)
%
%        form backprojected and whitened data
%
        gam=m/n;
        [ys4,wvals,as2_mean] = lintr_gen2whi(ys,as,m,n,k,var_ep);
        [u4,s4,v4] = svshr_svdsmartc(ys4,m,n,k);
        evals = s4.^2 / n;
%%%        svshr_compare_mp(evals,m/n,1,1)
%%%        stopnow

%
%        estimate cosines and population spikes of whitened signal
%
        for i=1:k
%
        [ells(i),cos_out(i),cos_inn(i)] = svshr_emp2pop_white(evals(i),gam);
    end

        ells_fr = zeros(k,1);
        for i=1:k
%
        ells_fr(i) = ells(i) * cos_out(i)^2;
    end

%
%        transform to the original coordinates and form covariance
%
        b=diag(1./wvals) * u4(:,1:k);
        xcov =  b * diag(ells_fr) * b';

%
%        estimate error (in whitened Frobenius norm)
%
        err_fr = 0;
        for i=1:k
%
        err_fr = err_fr + ells(i)^2 * (1 - cos_out(i)^4);
    end

        err_fr = sqrt(err_fr);

        end
%
%
%
%
%
        function xhat = wifi_approx_diag2(ys,as,u_true,ells,m,m2,n,k,sig)
%
        yback = conj(as).*ys;
        ws = zeros(k,n);
        for i=1:k
%
        bbb = conj(as).*repmat(u_true(:,i),1,n);
        ws(i,:) = sum(abs(bbb).^2)';
    end

        ellsr = repmat(ells,1,n);
        etas = ellsr ./ (ws .* ellsr + sig^2);
%%%        etas = repmat(ells,1,n) ./ (ws .* repmat(ells,1,n) + sig^2);
        xhat = u_true*(etas.*(u_true'*yback));


        end
%
%
%
%
%
        function xhat = wifi_approx_diag(y,a,u_true,ells,m,m2,k,sig)
%
        yback = conj(a).*y;
        au = repmat(conj(a),1,k).*u_true;
        ws = sum(abs(au).^2)';
%
        etas = ells ./ (ws.*ells + sig^2);
        pros = u_true'*yback;
        xhat = u_true*diag(etas)*pros;

        end
%
%
%
%
%
        function xhat = wifi_approx(y,amat,u_true,ells,m,m2,k,sig)
%
        yback = amat'*y;
        au = amat*u_true;
        ws = sum(abs(au).^2)';
%
        etas = ells ./ (ws.*ells + sig^2);
        pros = u_true'*yback;
        xhat = u_true*diag(etas)*pros;


        return

%
%        another way: do it with a loop over the components
%
        yback = amat'*y;
        xhat = zeros(m,1);

        au = amat*u_true;
        ws = sum(abs(au).^2);

        for j=1:k
%
        eta_j = ells(j) / (ws(j)*ells(j) + sig^2);
        pro_j = u_true(:,j)'*yback;
        xhat = xhat + eta_j*pro_j*u_true(:,j);
    end

        end
%
%
%
%
%
        function xhat = wifi_smart(y,amat,u_true,ells,m,m2,k,sig)
%
%%%        [uc,sc] = eig(c*c');
%%%        [sc,isort]=sort(diag(sc),'descend');
%%%        uc=uc(:,isort);
%%%        uc=uc(:,1:k);
%%%        sc=sc(1:k);

        c = amat*(u_true*diag(sqrt(ells)));
        [uc,sc,vc] = svshr_svdsmartc(c,m2,k,k);
        sc=sc.^2;
%
%        project y onto the column space of C
%
        yc = uc*(uc'*y);
        yc2 = y - yc;
%%%        chk0 = norm(y)^2 - (norm(yc)^2 + norm(yc2)^2)
%
%        apply the Wiener filter
%
        ucy = uc'*y;
        ducy = ucy ./ (sig^2 + sc);
        binv_y = uc*ducy + yc2/sig^2;
%
        xhat = u_true*diag(ells)*(u_true'*amat'*binv_y);



        return

        xcov = u_true*diag(ells)*u_true';
        [xhat2,fmat] = wifi_dumb(y,amat,xcov,m,m2,k,sig);
        chk0 = norm(xhat-xhat2)

        end
%
%
%
%
%
        function [xhat,fmat] = wifi_dumb(y,amat,xcov,m,m2,k,sig)
%
        b = amat*xcov*amat' + sig^2*eye(m2);
        b = .5*(b+b');
        [ub,sb] = eig(b);
        sb=diag(sb);
%
        sinv=1./sb;
        binv = ub*diag(sinv)*ub';

%%%        chk0 = norm(b*binv - eye(m2),'fro')

        fmat = xcov*amat'*binv;
        xhat = fmat*y;

        end
%
%
%
%
%
        function [ys,xs,ep,as,zs] = gendata_clintr2(m,m2,n,k,ells,u_true,var_ep)
%
%        generates complex data from spiked model
%
        imag1 = sqrt(-1)

        [xs,wnoise,zs] = draw_gaussian3(u_true,ells,m,n,k);
%
%        make complex gaussian noise
%
        alpr=rand();
        betr=sqrt(1-alpr^2);
        epr = randn(m2,n);
        epc = randn(m2,n);
        ep = alpr*epr + betr*imag1*epc;
%
        sig_ep = sqrt(var_ep);
        ep=repmat(sig_ep,1,n).*ep;

%
%        diagonal reduction matrices
%
%%%        amu = (2+rand(m,1)).^2 + imag1*rand(m,1);
        amu = rand(m2,1) + imag1*rand(m2,1);
        as = randn(m2,n) + imag1*randn(m2,n) + repmat(amu,1,n);


        acube = zeros(m2,m,n);
        asz = zeros(m,n);
        inds = zeros(m2,n);
        for i=1:n
%
        iis = randperm(m);
        inds(:,i) = sort(iis(1:m2),'ascend')';
        asz(inds(:,i),i) = as(:,i);

        for j=1:m2
%
        acube(j,inds(j,i),i) = as(j,i);
    end
    end


        ysz = asz.*xs;
        chk0 = 0;

        ys = zeros(m2,n);
        for i=1:n
%
        ys(:,i) = acube(:,:,i)*xs(:,i);
        chk0 = chk0 + norm(ysz(inds(:,i),i) - ys(:,i));
    end

        chk0

%%%        inds

%%%        as
%%%        asz
%%%        inds


%
        ys = as.*xs + ep;

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
%%%        amu = (2+rand(m,1)).^2 + imag1*rand(m,1);
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
