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
