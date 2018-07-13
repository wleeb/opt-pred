        function make_plots
%
        ifig = figure()

%
%        first box: no sparsity
%
        subplot(2,2,1)
        hold on;
        box on;

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

        title('$\kappa = 2$','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','OptShrink','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')


%
%        second box: half sparsity
%
        subplot(2,2,2)
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

        title('$\kappa = 4$','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','OptShrink','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')

%
%        third box: quarter sparsity
%
        subplot(2,2,3)
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

        title('$\kappa = 8$','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','OptShrink','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')


%
%        third box: quarter sparsity
%
        subplot(2,2,4)
        hold on;
        box on;

%%%        clear
        load compare4/compare4.mat
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

        title('$\kappa = 16$','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','OptShrink','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')


        set(figure(ifig),'Position',[500,500,1100,1000])


%%%        set(figure(ifig),'Position',[500,500,600,500])

        set(figure(ifig),'PaperPositionMode','auto')
        print('variances','-dpng','-r0')
        print(figure(ifig),'variances','-depsc','-r0')
        print(figure(ifig),'variances','-djpeg','-r0')
        print(figure(ifig),'variances','-dpdf','-r0')



        end
