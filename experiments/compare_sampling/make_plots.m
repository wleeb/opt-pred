        function make_plots
%
        ifig = figure()


%
%        second box: half sparsity
%
        subplot(1,3,1)
        hold on;
        box on;

%%%        clear
        load compare3/compare3.mat
        who

        errs_shr_mean = mean(errs_shr)
        errs_raj_mean = mean(errs_raj)
        errs_opt_mean = mean(errs_opt)
        errs_nuc_mean = mean(errs_nuc)
        
        xvals = log(sigmas)

        plot(xvals,log(errs_opt_mean),'x-','LineWidth',2,'Color','k')
        plot(xvals,log(errs_nuc_mean),'o-','LineWidth',2,'Color','m')
        plot(xvals,log(errs_raj_mean),'>-','LineWidth',2,'Color','r')
        plot(xvals,log(errs_shr_mean),'<-','LineWidth',2,'Color','b')

        xlim([min(xvals),max(xvals)])
        ylim([min(log(errs_opt_mean)),max(log(errs_opt_mean))])

        xlabel('$\log(\sigma)$','Interpreter','latex','FontSize',20)
        ylabel('$\log$-error','Interpreter','latex','FontSize',20)

        title('$\delta = 0.3$','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','OptShrink','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')

%
%        third box: quarter sparsity
%
        subplot(1,3,2)
        hold on;
        box on;

%%%        clear
        load compare2/compare2.mat
        who

        errs_shr_mean = mean(errs_shr)
        errs_raj_mean = mean(errs_raj)
        errs_opt_mean = mean(errs_opt)
        errs_nuc_mean = mean(errs_nuc)
        
        xvals = log(sigmas)

        plot(xvals,log(errs_opt_mean),'x-','LineWidth',2,'Color','k')
        plot(xvals,log(errs_nuc_mean),'o-','LineWidth',2,'Color','m')
        plot(xvals,log(errs_raj_mean),'>-','LineWidth',2,'Color','r')
        plot(xvals,log(errs_shr_mean),'<-','LineWidth',2,'Color','b')

        xlim([min(xvals),max(xvals)])
        ylim([min(log(errs_opt_mean)),max(log(errs_opt_mean))])

        xlabel('$\log(\sigma)$','Interpreter','latex','FontSize',20)
        ylabel('$\log$-error','Interpreter','latex','FontSize',20)

        title('$\delta = 0.2$','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','OptShrink','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')


%
%        third box: quarter sparsity
%
        subplot(1,3,3)
        hold on;
        box on;

%%%        clear
        load compare1/compare1.mat
        who

        errs_shr_mean = mean(errs_shr)
        errs_raj_mean = mean(errs_raj)
        errs_opt_mean = mean(errs_opt)
        errs_nuc_mean = mean(errs_nuc)
        
        xvals = log(sigmas)

        plot(xvals,log(errs_opt_mean),'x-','LineWidth',2,'Color','k')
        plot(xvals,log(errs_nuc_mean),'o-','LineWidth',2,'Color','m')
        plot(xvals,log(errs_raj_mean),'>-','LineWidth',2,'Color','r')
        plot(xvals,log(errs_shr_mean),'<-','LineWidth',2,'Color','b')

        xlim([min(xvals),max(xvals)])
        ylim([min(log(errs_opt_mean)),max(log(errs_opt_mean))])

        xlabel('$\log(\sigma)$','Interpreter','latex','FontSize',20)
        ylabel('$\log$-error','Interpreter','latex','FontSize',20)

        title('$\delta = 0.1$','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','OptShrink','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')


%%%        set(figure(ifig),'Position',[500,500,1100,1000])
%%%        set(figure(ifig),'Position',[500,500,600,500])


        set(figure(ifig),'Position',[500,500,1475,400])

        set(figure(ifig),'PaperPositionMode','auto')
        print('sampling','-dpng','-r0')
        print(figure(ifig),'sampling','-depsc','-r0')
        print(figure(ifig),'sampling','-djpeg','-r0')
        print(figure(ifig),'sampling','-dpdf','-r0')



        end
