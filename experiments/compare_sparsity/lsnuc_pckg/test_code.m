        function main
        startnow;

        randn(1,10)

        gam=.8
        m=200
        n=floor(m/gam)
        k=3
%
%       make a low-rank matrix a
%
        u = randn(m,k);
        v = randn(n,k);

        spec = zeros(1,k);
        for i=1:k
%
        for ijk=1:2
%
        for j=1:i-1
%
        puij = sum(u(:,i).*u(:,j));
        u(:,i) = u(:,i) - puij*u(:,j);

        pvij = sum(v(:,i).*v(:,j));
        v(:,i) = v(:,i) - pvij*v(:,j);
    end

        u(:,i) = u(:,i) / norm(u(:,i));
        v(:,i) = v(:,i) / norm(v(:,i));
    end

    end

        delta00=.4
        for i=1:k
%
        spec(i) = 100*sqrt(m/n)/delta00 + (k + 10  - i);
    end
        spec=sort(spec,'descend')
        a = u*diag(spec)*v';

%
%        generate indices to be observed and form subsampled matrix
%
        pmin=.5
        pmax=.5
        pvec=linspace(pmin,pmax,m)';
        ps = repmat(pvec,1,n);
%
        inds = rand_inds(m,n,ps);

%
%
%        make colored noise
%
%%%        var_ep = 1+2*rand(m,1);

        sig=2;
        var_ep = sig^2 * ones(m,1);

        ep=diag(sqrt(var_ep))*randn(m,n);

        y = (a+ep).*inds;

        y = ep.*inds;

        nsims=100;
        rlam = lsnuc_find_rlam(var_ep,ps,m,n,nsims)

        norm(ep.*inds)

        sig*(sqrt(n) + sqrt(m))*sqrt(pmin)


%%%        stopnow

        iis = find(inds == 1);
        len=length(iis);

%
%        weights matrix
%
        cs=ones(n,1);
        rs = pvec;
        rs = rs / min(rs);

        whts=zeros(m,n);
        for i=1:m
%
        for j=1:n
%
        whts(i,j)=1 / (rs(i)*cs(j));
    end
    end


%%%        [ww,nmax,whts2] = lsnuc_wrap_agm2(y,inds,m,n,rlam,w0,niter,thresh)


        niter=2000;
        thresh=1d-8
        w0=zeros(m,n);
        [ww,nmax] = lsnuc_wrap_agm(y,inds,m,n,w0,niter,thresh,sig);


%%%        check_kkt2(ww,y,m,n,iis,len,rlam,thresh*1000,whts2)
        check_kkt(ww,y,m,n,iis,len,rlam,thresh*1000)


        nmax
        niter

        stopnow
        end
%
%
%
%
%
        function chkbds_egm(w0,ww,vals,m,n,len,nmax,ell,ifig)
%
        err_sq=norm(w0 - ww,'fro')^2;

        bds=zeros(nmax-1,1);
        gaps=zeros(nmax-1,1);

        for i=2:nmax
%
        bds(i-1) = err_sq * ell / (2*i);
        gaps(i-1) = vals(i) - vals(nmax);
    end

        [chk_pos,imin] = min(bds-gaps);
        fprintf('Should be positive: \n     %d \n\n',chk_pos)

        if (ifig > 0)
%
        figure(ifig)
        plot(bds,'*')
        hold on;
        plot(gaps,'*')
    end
        end
%
%
%
%
%
        function chkbds_agm(w0,ww,vals,m,n,len,nmax,ell,ifig)
%
        err_sq=norm(w0 - ww,'fro')^2;

        bds=zeros(nmax-1,1);
        gaps=zeros(nmax-1,1);

        for i=2:nmax
%
        bds(i-1) = err_sq * 2 * ell / (i+1)^2;
        gaps(i-1) = vals(i) - vals(nmax);
    end

        [chk_pos,imin] = min(bds-gaps);
        fprintf('Should be positive: \n     %d \n\n',chk_pos)


        if (ifig > 0)
%
        figure(ifig)
        plot(bds,'*')
        hold on;
        plot(gaps,'*')
    end
        end
%
%
%
%
%
        function check_kkt2(ww,y,m,n,iis,len,rlam,thresh,whts)
%
%
%        check KKT conditions 
%
        [uw,sw,vw]=lsnuc_svdsmart(ww,m,n,min(m,n));
        chk0 = norm(uw*diag(sw)*vw' - ww)

        k=sum(sw>thresh);

        uw=uw(:,1:k);
        vw=vw(:,1:k);
        sw=sw(1:k);
        chk0 = norm(uw*diag(sw)*vw' - ww)

        grad = lsnuc_eval_grad2(ww,iis,y,len,m,n,whts);
        smat = grad/rlam + uw*vw';

        chk_lt1 = norm(smat);
        chk0_1 = norm(uw'*smat);
        chk0_2 = norm(smat*vw);

        fprintf('Should be zero: \n     %d \n\n',chk0_1)
        fprintf('Should be zero: \n     %d \n\n',chk0_2)
        fprintf('Should be less than 1: \n     %d \n\n',chk_lt1)

        gnorm = norm(grad);
        fprintf('Otherwise, should be zero: \n     %d \n\n',gnorm)

        end
%
%
%
%
%

        function check_kkt(ww,y,m,n,iis,len,rlam,thresh)
%
%
%        check KKT conditions 
%
        [uw,sw,vw]=lsnuc_svdsmart(ww,m,n,min(m,n));
        chk0 = norm(uw*diag(sw)*vw' - ww)

        k=sum(sw>thresh);

        uw=uw(:,1:k);
        vw=vw(:,1:k);
        sw=sw(1:k);
        chk0 = norm(uw*diag(sw)*vw' - ww)

        grad = lsnuc_eval_grad(ww,iis,y,len,m,n);
        smat = grad/rlam + uw*vw';

        chk_lt1 = norm(smat);
        chk0_1 = norm(uw'*smat);
        chk0_2 = norm(smat*vw);

        fprintf('Should be zero: \n     %d \n\n',chk0_1)
        fprintf('Should be zero: \n     %d \n\n',chk0_2)
        fprintf('Should be less than 1: \n     %d \n\n',chk_lt1)

        gnorm = norm(grad);
        fprintf('Otherwise, should be zero: \n     %d \n\n',gnorm)

        end
%
%
%
%
%
        function w = cvx_klopp(y,m,n,inds,iis,len,rlam,val_max)
%
        cvx_begin
        cvx_precision('low')
        variable  w(m,n)

        minimize( sum_square_abs( w(iis)-y(iis) ) + rlam*norm_nuc(w) );
        subject to 
            max(abs(w(:))) <= val_max;

        cvx_end


        end
%
%
%
%
%
        function w = cvx_solve2(y,m,n,inds,iis,len,rlam,whts)
%
        cvx_begin
        cvx_precision('high')
%%%        cvx_precision('low')
        variable  w(m,n)
        minimize( .5*sum_square_abs( w(iis).*whts(iis)-y(iis) ) + rlam*norm_nuc(w) )

        cvx_end


        end
%
%
%
%
%
        function w = cvx_solve(y,m,n,inds,iis,len,rlam)
%
        cvx_begin
%%%        cvx_precision('high')
        cvx_precision('low')
        variable  w(m,n)
        minimize( sum_square_abs( w(iis)-y(iis) ) + rlam*norm_nuc(w) )

        cvx_end


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
        function startnow
%
        delete out13
        diary('out13')
        diary on
%
        format short E
%%%        format long E


%
        rng('default');

        end
%
%
%
%
%
        function stopnow
%
        diary off
        stop

        end
%
%
%
%
%
