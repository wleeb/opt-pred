        function make_plots
%
%%%        plot1
        plot2

        end
%
%
%
%
%
        function plot1
%
        ifig = figure()

%
%        first box: density=1
%
        subplot(2,2,1)
        hold on;
        box on;

        load oos1/oos1.mat
        who

        xmin=min(min(errs_out),min(errs_in))
        xmax=max(max(errs_out),max(errs_in))

        xx=linspace(xmin,xmax,300);
        plot(xx,xx,'LineWidth',2,'Color','g')

        scatter(errs_out,errs_in,100,'k.')

        xlim([xmin,xmax])
        ylim([xmin,xmax])


        xlabel('Out-of-sample RMSE','Interpreter','latex','FontSize',20)
        ylabel('In-sample RMSE','Interpreter','latex','FontSize',20)

        title('Density = 1','Interpreter','latex','FontSize',20)

%
%        second box: density=.7
%
        subplot(2,2,2)
        hold on;
        box on;

        load oos2/oos2.mat
        who

        xmin=min(min(errs_out),min(errs_in))
        xmax=max(max(errs_out),max(errs_in))


        xx=linspace(xmin,xmax,300);
        plot(xx,xx,'LineWidth',2,'Color','g')

        scatter(errs_out,errs_in,100,'k.')

        xlim([xmin,xmax])
        ylim([xmin,xmax])



        xlabel('Out-of-sample RMSE','Interpreter','latex','FontSize',20)
        ylabel('In-sample RMSE','Interpreter','latex','FontSize',20)

        title('Density = 0.7','Interpreter','latex','FontSize',20)


%
%        third box: density=.3
%
        subplot(2,2,3)
        hold on;
        box on;

        load oos3/oos3.mat
        who

        xmin=min(min(errs_out),min(errs_in))
        xmax=max(max(errs_out),max(errs_in))


        xx=linspace(xmin,xmax,300);
        plot(xx,xx,'LineWidth',2,'Color','g')

        scatter(errs_out,errs_in,100,'k.')

        xlim([xmin,xmax])
        ylim([xmin,xmax])



        xlabel('Out-of-sample RMSE','Interpreter','latex','FontSize',20)
        ylabel('In-sample RMSE','Interpreter','latex','FontSize',20)

        title('Density = 0.3','Interpreter','latex','FontSize',20)



%
%        fourth box: density=.1
%
        subplot(2,2,4)
        hold on;
        box on;

        load oos4/oos4.mat
        who

        xmin=min(min(errs_out),min(errs_in))
        xmax=max(max(errs_out),max(errs_in))


        xx=linspace(xmin,xmax,300);
        plot(xx,xx,'LineWidth',2,'Color','g')

        scatter(errs_out,errs_in,100,'k.')

        xlim([xmin,xmax])
        ylim([xmin,xmax])



        xlabel('Out-of-sample RMSE','Interpreter','latex','FontSize',20)
        ylabel('In-sample RMSE','Interpreter','latex','FontSize',20)

        title('Density = 0.1','Interpreter','latex','FontSize',20)



        set(figure(ifig),'Position',[500,500,1100,1000])
%%%        set(figure(ifig),'Position',[500,500,900,800])



%%%        set(figure(ifig),'Position',[500,500,600,500])

        set(figure(ifig),'PaperPositionMode','auto')
        print('oos_rank10','-dpng','-r0')
        print(figure(ifig),'oos_rank10','-depsc','-r0')
        print(figure(ifig),'oos_rank10','-djpeg','-r0')
        print(figure(ifig),'oos_rank10','-dpdf','-r0')



        end
%
%
%
%
%
        function plot2
%
        ifig = figure()

%
%        first box: no subsampling
%
        subplot(1,3,1)
        hold on;
        box on;

        load oos8/oos8.mat
        who

        xmin=min(min(errs_out),min(errs_in))
        xmax=max(max(errs_out),max(errs_in))

        xx=linspace(xmin,xmax,300);
        plot(xx,xx,'LineWidth',2,'Color','g')

        scatter(errs_out,errs_in,100,'k.')

        xlim([xmin,xmax])
        ylim([xmin,xmax])


        xlabel('Out-of-sample RMSE','Interpreter','latex','FontSize',20)
        ylabel('In-sample RMSE','Interpreter','latex','FontSize',20)

        title('Density = 1','Interpreter','latex','FontSize',20)



        if(0)
%
%        second box: density=.7
%
        subplot(2,2,2)
        hold on;
        box on;

        load oos7/oos7.mat
        who

        xmin=min(min(errs_out),min(errs_in))
        xmax=max(max(errs_out),max(errs_in))


        xx=linspace(xmin,xmax,300);
        plot(xx,xx,'LineWidth',2,'Color','g')

        scatter(errs_out,errs_in,100,'k.')

        xlim([xmin,xmax])
        ylim([xmin,xmax])



        xlabel('Out-of-sample RMSE','Interpreter','latex','FontSize',20)
        ylabel('In-sample RMSE','Interpreter','latex','FontSize',20)

        title('Density = 0.7','Interpreter','latex','FontSize',20)
    end

%
%        third box: density=.3
%
        subplot(1,3,2)
        hold on;
        box on;

        load oos6/oos6.mat
        who

        xmin=min(min(errs_out),min(errs_in))
        xmax=max(max(errs_out),max(errs_in))


        xx=linspace(xmin,xmax,300);
        plot(xx,xx,'LineWidth',2,'Color','g')

        scatter(errs_out,errs_in,100,'k.')

        xlim([xmin,xmax])
        ylim([xmin,xmax])



        xlabel('Out-of-sample RMSE','Interpreter','latex','FontSize',20)
        ylabel('In-sample RMSE','Interpreter','latex','FontSize',20)

        title('Density = 0.3','Interpreter','latex','FontSize',20)



%
%        fourth box: density=.1
%
        subplot(1,3,3)
        hold on;
        box on;

        load oos5/oos5.mat
        who

        xmin=min(min(errs_out),min(errs_in))
        xmax=max(max(errs_out),max(errs_in))


        xx=linspace(xmin,xmax,300);
        plot(xx,xx,'LineWidth',2,'Color','g')

        scatter(errs_out,errs_in,100,'k.')

        xlim([xmin,xmax])
        ylim([xmin,xmax])



        xlabel('Out-of-sample RMSE','Interpreter','latex','FontSize',20)
        ylabel('In-sample RMSE','Interpreter','latex','FontSize',20)

        title('Density = 0.1','Interpreter','latex','FontSize',20)



%%%        set(figure(ifig),'Position',[500,500,1100,1000])


        set(figure(ifig),'Position',[500,500,1475,400])


%%%        set(figure(ifig),'Position',[500,500,600,500])

        set(figure(ifig),'PaperPositionMode','auto')
        print('oos_rank1','-dpng','-r0')
        print(figure(ifig),'oos_rank1','-depsc','-r0')
        print(figure(ifig),'oos_rank1','-djpeg','-r0')
        print(figure(ifig),'oos_rank1','-dpdf','-r0')



        end
