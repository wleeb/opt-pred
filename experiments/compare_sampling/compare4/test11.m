        function main
        startnow;

        ifplot=0
        if(ifplot==1)
%
        make_plot()
        stopnow
    end

        gam=.8
        m=300
        n=floor(m/gam)

        k=10

        sig0 = 1;
        delta0 = .5;
        ells = floor(sig0*sqrt(gam) / delta0 + [k:-1:1]) - 1

        delta=delta0
%
%        sampling mask
%
        nruns = 40
        nsigmas = 10

        errs_shr = zeros(nruns,nsigmas);
        errs_opt = zeros(nruns,nsigmas);
        errs_nuc = zeros(nruns,nsigmas);

        sigmas = exp(linspace(log(1),log(7),nsigmas))

        tic

        nsims=20
        pvec = linspace(.4,.6,m)';
        var_ep = ones(m,1);
        bedge = lintr_edge_diag(pvec,var_ep,m,n,nsims);

        parfor irun = 1:nruns
%
        rng(irun+1);
        [x,ep0,imask] = make_data(m,n,ells,k,pvec);

        bedge = lintr_edge_diag(pvec,var_ep,m,n,nsims);
        bedge=1.1*bedge

        for ijk=1:nsigmas
%
        sig=sigmas(ijk);

        xnoisy = x+sig*ep0;
        y = xnoisy .* imask;

        [x_est,x_est2,x_opt,x_nuc,x_raj] = test_all(y,imask,m,n,k,sig,delta,pvec,bedge);

        errs_shr(irun,ijk) = norm(x_est - x,'fro') / norm(x,'fro');
        errs_shr2(irun,ijk) = norm(x_est2 - x,'fro') / norm(x,'fro');
        errs_raj(irun,ijk) = norm(x_raj - x,'fro') / norm(x,'fro');
        errs_opt(irun,ijk) = norm(x_opt - x,'fro') / norm(x,'fro');
        errs_nuc(irun,ijk) = norm(x_nuc - x,'fro') / norm(x,'fro');


%%%        stopnow

    end
    end

        timing=toc

        errs_opt
        errs_raj
        errs_shr
        errs_shr2
        errs_nuc


        save('compare4.mat','errs_shr','errs_shr2','errs_raj','errs_opt',...
            'errs_nuc','m','n','sigmas','nsigmas','nruns','delta','k','timing')


        make_plot()

        stopnow
        end
%
%
%
%
%

        function make_plot()
%
        load compare4.mat
        who
        timing

        ifig = figure()
        hold on;

        errs_shr_mean = mean(errs_shr,1)
        errs_shr2_mean = mean(errs_shr2,1)
        errs_raj_mean = mean(errs_raj,1)
        errs_opt_mean = mean(errs_opt,1)
        errs_nuc_mean = mean(errs_nuc,1)


        xvals = log(sigmas)


        plot(xvals,log(errs_opt_mean),'x-','LineWidth',2,'Color','r')
        plot(xvals,log(errs_nuc_mean),'o-','LineWidth',2,'Color','m')
        plot(xvals,log(errs_raj_mean),'o-','LineWidth',2,'Color','g')
        plot(xvals,log(errs_shr_mean),'<-','LineWidth',2,'Color','b')
%%%        plot(xvals,log(errs_shr2_mean),'<-','LineWidth',2,'Color','c')

        xlim([min(xvals)-.05,max(xvals)+.05])
        ylim([min(log(errs_opt_mean)),max(log(errs_opt_mean))])

        xlabel('$\log(\sigma)$','Interpreter','latex','FontSize',20)
        ylabel('$\log$-error','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','OptShrink','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')

        set(figure(ifig),'Position',[500,500,600,500])



        end
%
%
%
%
%
        function [x,ep0,imask] = make_data(m,n,ells,k,pvec)
%

        u = randn(m,k);


        u(:,1)=u(:,1) / norm(u(:,1));

        for i=2:k

        for ijk=1:2
        for j=1:i-1
%
        pij = sum(u(:,i).*u(:,j));
        u(:,i) = u(:,i) - pij*u(:,j);
    end
        u(:,i) = u(:,i) / norm(u(:,i));

    end
    end
        chk0 = norm(u'*u - eye(k))
        zs = randn(n,k);

        ps = repmat(pvec,1,n);

        ff = rand(m,n);
        imask = real(ff <= ps);

        x = u*diag(ells)*zs';

        xmu = 10*randn(m,1) / sqrt(m);
        x = x + repmat(xmu,1,n);
        ep0=randn(m,n);



        end
%
%
%
%
%
        function [x_est,x_est2,x_opt,x_nuc,x_raj] = test_all(y,imask,m,n,k,...
            sig,delta,pvec,bedge)
%
%
%        run shrinkage with known rank
%
        as = imask;
        [xmean,as2_mean,yback] = lintr_mean_diag(y,as,m,n);
        y2=y-repmat(xmean,1,n).*as;


%%%        y2=y;
%%%        xmean=zeros(m,1);

        var_ep = sig^2*ones(m,1);

        [x_est,whts,errs] = lintr_whit(y2,as,m,n,k,var_ep);
        x_est = x_est + repmat(xmean,1,n);


%
%        run shrinkage with unknown rank
%
        kest77 = min(m,n);
        [x_est2,whts,errs] = lintr_whit(y2,as,m,n,kest77,var_ep);
        x_est2 = x_est2 + repmat(xmean,1,n);


%
%        run OptShrink
%
        bedge=sig^2*bedge;
        [x_raj,errs] = lintr_unif2(y2,as,m,n,k,bedge);
        x_raj = x_raj + repmat(xmean,1,n);



        x_nuc=zeros(m,n);
        x_opt=zeros(m,n);

%%%        return

%
%        nuclear norm minimization
%
        thresh=1d-8;
        w0=x_est;
        niter=2000;
        [x_nuc,nmax,rlam,whts] = lsnuc_wrap_agm2(y,imask,m,n,w0,niter,thresh,pvec,sig);

%
%        OptSpace
%
        niter=[]
        tol=[]
        kest=k+1;
        [x_opt,u_opt,s_opt,v_opt] = call_optspace(y,...
           m,n,kest,niter,tol);




        end
%
%
%
%
%
        function [xs_opt,u_opt,s_opt,v_opt] = call_optspace(ys,...
           m,n,k,niter,tol)

        ys_sparse = sparse(ys);
        [u_opt,s_opt,v_opt,dist] = OptSpace(ys_sparse,k,niter,tol);

        xs_opt = u_opt * s_opt * v_opt';


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
