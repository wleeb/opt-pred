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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This is the end of the testing code and the beginning of the nuclear norm
%   minimization codes proper. This file contains two user-callable functions:
%   lsnuc_agm and lsnuc_egm. They both minimize the objective
%
%                     || P(A - Y) ||_F^2 + L * ||A||_*
%
%   where P projects an m-by-n matrix onto a given collection of coordinates, whose
%   values are observed. The two methods are extended gradient method (EGM) and
%   accelerated gradient method (AGM), both described in the paper ``An Accelerated 
%   Gradient Method for Trace Norm Minimization'' by Ji and Ye. The EGM has error on the
%   k^th iteration decaying like 1/k, while the AGM has error on the k^th iteration 
%   decaying like 1/k^2.
%
%   Additionally, there is a function lsnuc_agm2, which minimizes the objective
%
%                     || P(RAC - Y) ||_F^2 + L * ||A||_*
%
%   where R and C are specified diagonal matrices (specified by their inverse values).
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
        function [ww,nmax,rlam] = lsnuc_wrap_agm(y,inds,m,n,w0,niter,thresh,sig)
%
        iis = find(inds == 1);
        len=length(iis);

        pp = len / (m*n);

        rlam = sig*(sqrt(n) + sqrt(m))*sqrt(pp);
        [ws,nmax,vals] = lsnuc_agm(y,m,n,iis,len,rlam,w0,niter,thresh^2);
        ww=ws(:,:,nmax);


        end
%
%
%
%
%
        function [ww,nmax,whts] = lsnuc_wrap_agm2(y,inds,m,n,rlam,w0,niter,thresh)
%
        iis = find(inds == 1);
        len=length(iis);

        rs=mean(inds,2);
        rs = rs / min(rs);

%%%        cs=mean(inds,1);
%%%        cs = cs / min(cs);
        cs = ones(n,1);

        [ws,nmax,vals,whts] = lsnuc_agm2(y,m,n,iis,len,rlam,w0,niter,thresh^2,rs,cs);
        ww=ws(:,:,nmax);

        end
%
%
%
%
%
        function [ws,nmax,vals,whts] = lsnuc_agm2(y,m,n,iis,len,rlam,w0,niter,thresh,...
            rs,cs)
%
%
%                            description:
%
%   This code minimizes the convex objective
%
%                     || P(A - Y) ||_F^2 + L * ||A||_*
%
%   where P projects an m-by-n matrix onto the coordinates iis, defined columnwise.
%   It uses the accelerated gradient method defined in the paper ``An Accelerated 
%   Gradient Method for Trace Norm Minimization'' by Ji and Ye.  The error on the
%   k^th iteration decays like 1/k^2.
%
%
%                           input parameters:
%
%   y - the m-by-n data matrix of observations, with zeros filled in for
%      the unobserved entries
%
%   m,n - the dimensions
%
%   iis - a vector of length len, containing the indices of the observed values,
%      stored columnwise.
%
%   len - the number of observed entries
%
%   rlam - the weight placed on the nuclear norm term in the optimization
%
%   ell - a parameter defining the Lipschitz constant of the smooth part. Should
%      be greater than or equal to 2.
%
%   w0 - an m-by-n matrix. The initial guess of the optimization.
%
%   niter - the maximum number of iterations
%
%   thresh - the precision. When the norm between consecutive terms is below
%      thresh, the procedure terminates
%
%
%                        output parameters:
%
%   ws - the m-by-n-by-nmax cube of matrices in the iteration, starting with w0
%
%   nmax - the number of iterations until convergence
%
%   vals - the values of the objective function (squared error plus nuclear norm) at
%      the matrices in ws
%
%
        alps=zeros(niter,1);
        alps(1)=1;

        ws=zeros(m,n,niter);
        zs=zeros(m,n,niter);

        z=w0;
        w=w0;

%
%        weights matrix
%
        whts=zeros(m,n);
        for i=1:m
%
        for j=1:n
%
        whts(i,j)=1 / (rs(i)*cs(j));
    end
    end
        ell = max(whts(iis))^2

        [vals(1),fvals(1),rnucs(1)] = lsnuc_eval_obj2(w0,iis,len,y,m,n,rlam,whts);

        ws(:,:,1)=w0;
        zs(:,:,1)=w0;

        nmax=niter;
        vvs=zeros(1,niter);
        uus=zeros(1,niter);

        for ijk=2:niter
%
        alps(ijk) = 1 + sqrt(1 + 4*alps(ijk-1)^2);
        alps(ijk) = alps(ijk)/2;

        w = lsnuc_prox_min2(z,ell,rlam,y,iis,len,m,n,whts);
        ws(:,:,ijk) = w;

        coef = (alps(ijk-1) - 1) / alps(ijk);
        z = w + coef*(w - ws(:,:,ijk-1));
        zs(:,:,ijk) = z;
        coefs(ijk) = coef;

%
%        check convergence
%
        [vals(ijk),fvals(ijk),rnucs(ijk)] = lsnuc_eval_obj2(w,iis,len,y,m,n,rlam,whts);
%%%        dif = norm(ws(:,:,ijk) - ws(:,:,ijk-1),'fro') / norm(ws(:,:,ijk),'fro');
        dif2 = abs(vals(ijk) - vals(ijk-1)) / vals(ijk-1);

        if (dif2 < thresh)
           nmax=ijk;
           break;
        end
    end

        ws=ws(:,:,1:nmax);
        vals=vals(1:nmax);


        end
%
%
%
%
%
        function wmin = lsnuc_prox_min2(w1,tt,rlam,y,iis,len,m,n,whts)
%
        grad = lsnuc_eval_grad2(w1,iis,y,len,m,n,whts);
        bb = w1 - grad / tt;
        thr = rlam / tt;
        wmin = lsnuc_softhresh(bb,m,n,thr);


        end
%
%
%
%
%
        function fval = lsnuc_eval_fsmall2(a,iis,y,len,m,n,whts)
%
%        evaluates the smooth term in the objective function.
%        y is the observed data, a is the variable matrix
%
        dd = a(iis) .* whts(iis);
        fval = norm(dd - y(iis))^2 / 2;

        end
%
%
%
%
%
        function grad = lsnuc_eval_grad2(a,iis,y,len,m,n,whts)
%
%        evaluates gradient of little f (fsmall), the smooth term
%        in the objective function (the Frobenius norm)
%        y is the observed data, a is the variable matrix
%
        grad = zeros(m,n);

        dd=whts(iis).*a(iis);
        grad(iis) = (dd - y(iis)) .*  whts(iis);

%%%        grad(iis) = 2*(a(iis) - y(iis));

        end
%
%
%
%
%
        function [val,fval,rnuc] = lsnuc_eval_obj2(a,iis,len,y,m,n,rlam,whts)
%
%        evaluates the objective function of the optimization; squared frobenius
%        loss on observed entries plus rlam times nuclear norm.
%
        val=0;
        rnuc=0;
        fval = lsnuc_eval_fsmall2(a,iis,y,len,m,n,whts);

        k=min(m,n);
        [ua,sa,va] = lsnuc_svdsmart(a,m,n,k);

%%%        chk0 = norm(ua*diag(sa)*va' - a,'fro')

        rnuc = sum(sa);
        val = fval + rlam*rnuc;

        end
%
%
%
%
%
        function rlam = lsnuc_find_rlam(var_ep,ps,m,n,nsims)
%
%
%        approximates rlam as the operator norm of noise term,
%        depending on variances and sampling probabilities
%
        rlam = 0;
        rnorms=zeros(1,nsims);
        for i=1:nsims
%
        inds00 = lsnuc_rand_inds(m,n,ps);
        ep00 = diag(sqrt(var_ep))*randn(m,n);
        ep00 = ep00.*inds00;
        rnorms(i) = norm(ep00);
    end
        rlam=median(rnorms);

        end
%
%
%
%
%
        function inds = lsnuc_rand_inds(m,n,ps)
%
        ff = rand(m,n);
        inds = real(ff <= ps);

        end
%
%
%
%
%
        function [ws,nmax,vals] = lsnuc_agm(y,m,n,iis,len,rlam,w0,niter,thresh)
%
%
%                            description:
%
%   This code minimizes the convex objective
%
%                     || P(A - Y) ||_F^2 + L * ||A||_*
%
%   where P projects an m-by-n matrix onto the coordinates iis, defined columnwise.
%   It uses the accelerated gradient method defined in the paper ``An Accelerated 
%   Gradient Method for Trace Norm Minimization'' by Ji and Ye.  The error on the
%   k^th iteration decays like 1/k^2.
%
%
%                           input parameters:
%
%   y - the m-by-n data matrix of observations, with zeros filled in for
%      the unobserved entries
%
%   m,n - the dimensions
%
%   iis - a vector of length len, containing the indices of the observed values,
%      stored columnwise.
%
%   len - the number of observed entries
%
%   rlam - the weight placed on the nuclear norm term in the optimization
%
%   ell - a parameter defining the Lipschitz constant of the smooth part. Should
%      be greater than or equal to 2.
%
%   w0 - an m-by-n matrix. The initial guess of the optimization.
%
%   niter - the maximum number of iterations
%
%   thresh - the precision. When the norm between consecutive terms is below
%      thresh, the procedure terminates
%
%
%                        output parameters:
%
%   ws - the m-by-n-by-nmax cube of matrices in the iteration, starting with w0
%
%   nmax - the number of iterations until convergence
%
%   vals - the values of the objective function (squared error plus nuclear norm) at
%      the matrices in ws
%
%
        alps=zeros(niter,1);
        alps(1)=1;

        ws=zeros(m,n,niter);
        zs=zeros(m,n,niter);

        z=w0;
        w=w0;

        [vals(1),fvals(1),rnucs(1)] = lsnuc_eval_obj(w0,iis,len,y,m,n,rlam);

        ws(:,:,1)=w0;
        zs(:,:,1)=w0;

        ell=1
        nmax=niter;

        vvs=zeros(1,niter);
        uus=zeros(1,niter);

        for ijk=2:niter
%
        alps(ijk) = 1 + sqrt(1 + 4*alps(ijk-1)^2);
        alps(ijk) = alps(ijk)/2;

        w = lsnuc_prox_min(z,ell,rlam,y,iis,len,m,n);
        ws(:,:,ijk) = w;

        coef = (alps(ijk-1) - 1) / alps(ijk);
        z = w + coef*(w - ws(:,:,ijk-1));
        zs(:,:,ijk) = z;
        coefs(ijk) = coef;
        [vals(ijk),fvals(ijk),rnucs(ijk)] = lsnuc_eval_obj(w,iis,len,y,m,n,rlam);

%%%        dif = norm(ws(:,:,ijk) - ws(:,:,ijk-1),'fro') / norm(ws(:,:,ijk),'fro');
        dif2 = abs(vals(ijk) - vals(ijk-1)) / vals(ijk);

        if (dif2 < thresh)
           nmax=ijk;
           break;
        end

    end

        ws=ws(:,:,1:nmax);
        vals=vals(1:nmax);


        end
%
%
%
%
%
        function [ws,nmax,vals] = lsnuc_egm(y,m,n,iis,len,rlam,...
           w0,niter,thresh)
%
%
%                            description:
%
%   This code minimizes the convex objective
%
%                     || P(A - Y) ||_F^2 + L * ||A||_*
%
%   where P project an m-by-n matrix onto the coordinates iis, defined columnwise.
%   It uses the extended gradient method defined in the paper ``An Accelerated 
%   Gradient Method for Trace Norm Minimization'' by Ji and Ye. The error on the
%   k^th iteration decays like 1/k.
%
%
%                           input parameters:
%
%   y - the m-by-n data matrix of observations, with zeros filled in for
%      the unobserved entries
%
%   m,n - the dimensions
%
%   iis - a vector of length len, containing the indices of the observed values,
%      stored columnwise.
%
%   len - the number of observed entries
%
%   rlam - the weight place on the nuclear norm term in the optimization
%
%   ell - a parameter defining the Lipschitz constant of the smooth part. Should
%      be greater than or equal to 2.
%
%   w0 - an m-by-n matrix. The initial guess of the optimization.
%
%   niter - the maximum number of iterations
%
%   thresh - the precision. When the norm between consecutive terms is below
%      thresh, the procedure terminates
%
%
%                        output parameters:
%
%   ws - the m-by-n-by-nmax cube of matrices in the iteration, starting with w0
%
%   nmax - the number of iterations until convergence
%
%   vals - the values of the objective function (squared error plus nuclear norm) at
%      the matrices in ws
%
%

        ws=zeros(m,n,niter);
        vals=zeros(niter,1);
        fvals=zeros(niter,1);
        rnucs=zeros(niter,1);
%
        w=w0;
        ws(:,:,1) = w;
        [vals(1),fvals(1),rnucs(1)] = lsnuc_eval_obj(w,iis,len,y,m,n,rlam);
        
        ell=1;

        nmax=niter;
        for ijk=2:niter
%
        w = lsnuc_prox_min(w,ell,rlam,y,iis,len,m,n);
        ws(:,:,ijk) = w;
        [vals(ijk),fvals(ijk),rnucs(ijk)] = lsnuc_eval_obj(w,iis,len,y,m,n,rlam);

        dif = norm(ws(:,:,ijk) - ws(:,:,ijk-1),'fro');
        if (dif < thresh)
           nmax=ijk;
           break;
        end

    end

        ws=ws(:,:,1:nmax);
        vals = vals(1:nmax);

        end
%
%
%
%
%
        function fval = lsnuc_eval_fsmall(a,iis,y,len,m,n,rlam)
%
%        evaluates the smooth term in the objective function.
%        y is the observed data, a is the variable matrix
%
        fval = norm(a(iis) - y(iis))^2 / 2;

        end
%
%
%
%
%
        function grad = lsnuc_eval_grad(a,iis,y,len,m,n)
%
%        evaluates gradient of little f (fsmall), the smooth term
%        in the objective function (the Frobenius norm)
%        y is the observed data, a is the variable matrix
%
        grad = zeros(m,n);
        grad(iis) = a(iis) - y(iis);

        end
%
%
%
%
%
        function [val,fval,rnuc] = lsnuc_eval_obj(a,iis,len,y,m,n,rlam)
%
%        evaluates the objective function of the optimization; squared frobenius
%        loss on observed entries plus rlam times nuclear norm.
%
        val=0;
        rnuc=0;
        fval = lsnuc_eval_fsmall(a,iis,y,len,m,n,rlam);

        k=min(m,n);
        [ua,sa,va] = lsnuc_svdsmart(a,m,n,k);

%%%        chk0 = norm(ua*diag(sa)*va' - a,'fro')

        rnuc = sum(sa);
        val = fval + rlam*rnuc;

        end
%
%
%
%
%
        function [u,s,v] = lsnuc_svdsmart(a,m,n,k)
%
        if (min(m,n)/k > 20)
%
        [u,s,v] = svds(a,k);
        s=diag(s);
    end

        if (min(m,n)/k <= 20)
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
        function wmin = lsnuc_prox_min(w1,tt,rlam,y,iis,len,m,n)
%
        grad = lsnuc_eval_grad(w1,iis,y,len,m,n);
        bb = w1 - grad / tt;
        thr = rlam / tt;
        wmin = lsnuc_softhresh(bb,m,n,thr);


        end
%
%
%
%
%
        function a_thr = lsnuc_softhresh(a,m,n,thr)
%
        k_thr = min(m,n);
        k=min(m,n);
        [ua,sa,va] = lsnuc_svdsmart(a,m,n,k);
%%%        chk0 = norm(ua*diag(sa)*va' - a,'fro')

        s_thr = 0*sa;

        s_thr = max(sa-thr,0);
        k_thr = sum(s_thr>0);

        s_thr=s_thr(1:k_thr);
        a_thr = ua(:,1:k_thr) * diag(s_thr) * va(:,1:k_thr)';



        return


        for i=1:k
%
        s_thr(i) = sa(i) - thr;

        if (s_thr(i) < 0)
           s_thr(i) = 0;
           k_thr = i - 1;
           break;
        end

    end
%
        s_thr = s_thr(1:k_thr);
        a_thr = ua(:,1:k_thr) * diag(s_thr) * va(:,1:k_thr)';

        end
