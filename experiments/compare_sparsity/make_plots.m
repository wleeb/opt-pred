        function make_plots
%
        ifig = figure()

%
%        first box: no sparsity
%

        load compare1/compare1.mat
        who


        subplot(1,3,1)
        hold on;
        box on;

        errs_shr_mean = mean(errs_shr)
        errs_opt_mean = mean(errs_opt)
        errs_nuc_mean = mean(errs_nuc)
        
        xvals = log(sigmas)

        plot(xvals,log(errs_opt_mean),'x-','LineWidth',2,'Color','k')
        plot(xvals,log(errs_nuc_mean),'o-','LineWidth',2,'Color','m')
        plot(xvals,log(errs_shr_mean),'<-','LineWidth',2,'Color','b')


        xlim([min(xvals),max(xvals)])
        ylim([min(log(errs_opt_mean)),max(log(errs_opt_mean))])

        xlabel('$\log(\sigma)$','Interpreter','latex','FontSize',20)
        ylabel('$\log$-error','Interpreter','latex','FontSize',20)

        title('Full vector','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')

%
%        second box: quarter sparsity
%
        subplot(1,3,2)
        hold on;
        box on;

%%%        clear
        load compare5/compare5.mat
        who

        errs_shr_mean = mean(errs_shr)
        errs_opt_mean = mean(errs_opt)
        errs_nuc_mean = mean(errs_nuc)
        
        xvals = log(sigmas)

        plot(xvals,log(errs_opt_mean),'x-','LineWidth',2,'Color','k')
        plot(xvals,log(errs_nuc_mean),'o-','LineWidth',2,'Color','m')
        plot(xvals,log(errs_shr_mean),'<-','LineWidth',2,'Color','b')

        xlim([min(xvals),max(xvals)])
        ylim([min(log(errs_opt_mean)),max(log(errs_opt_mean))])

        xlabel('$\log(\sigma)$','Interpreter','latex','FontSize',20)
        ylabel('$\log$-error','Interpreter','latex','FontSize',20)

        title('$p/4$-sparse','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')


%
%        third box: 10 sparsity
%
        subplot(1,3,3)
        hold on;
        box on;

%%%        clear
        load compare4/compare4.mat
        who

        errs_shr_mean = mean(errs_shr)
        errs_opt_mean = mean(errs_opt)
        errs_nuc_mean = mean(errs_nuc)
        
        xvals = log(sigmas)

        plot(xvals,log(errs_opt_mean),'x-','LineWidth',2,'Color','k')
        plot(xvals,log(errs_nuc_mean),'o-','LineWidth',2,'Color','m')
        plot(xvals,log(errs_shr_mean),'<-','LineWidth',2,'Color','b')

        xlim([min(xvals),max(xvals)])
        ylim([min(log(errs_opt_mean)),max(log(errs_opt_mean))])

        xlabel('$\log(\sigma)$','Interpreter','latex','FontSize',20)
        ylabel('$\log$-error','Interpreter','latex','FontSize',20)

        title('$10$-sparse','Interpreter','latex','FontSize',20)

        legend({'OptSpace','Nuclear norm','EBLP'},'FontSize',15,...
           'Interpreter','latex','Location','NorthWest')


%%%        set(figure(ifig),'Position',[500,500,1100,1000])

        set(figure(ifig),'Position',[500,500,1475,400])

%%%        set(figure(ifig),'Position',[500,500,600,500])

        set(figure(ifig),'PaperPositionMode','auto')
        print('sparsity','-dpng','-r0')
        print(figure(ifig),'sparsity','-depsc','-r0')
        print(figure(ifig),'sparsity','-djpeg','-r0')
        print(figure(ifig),'sparsity','-dpdf','-r0')



        end
