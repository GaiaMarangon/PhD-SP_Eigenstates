function [rmax,fmax,m,q,r2,qbis,r2bis,c,cbis]=fitPeaks(r,absf,k,eps,ifPeak,ifLin,savefigPeak,savefigLin,rejend,outputPath)
    %fit with (ordinary) least squares method to function fmax(rmax)= c*rmax^a
    %fit with (ordinary) least squares method to function fmax(rmax)= c*rmax^-1 (a=-1)
    %----------------------------------------------------
    %setting
    rejbeg=1; %min=1 (r=0 sempre rigettato)
    mark='o';
    markSize=4;
    %----------------------------------------------------
    %find peaks (local maxima) of f(r)
    [rmax,fmax,~]=findLocMaxima(r,absf,k);
    %move to log quantities (fit is now linear)
    nmax=length(rmax);
    xfit=log(rmax(rejbeg+1:nmax-rejend));
    yfit=log(fmax(rejbeg+1:nmax-rejend));
    %----------------------------------------------------
    %least squares coefficients (param a,c)
    n=length(xfit);
    sumx = sum(xfit);
    sumy = sum(yfit);
    avy = sumy/n;
    m = (n*(sum(xfit.*yfit))-sumx*sumy)/(n*sum(xfit.^2)-sumx^2);
    q = (sumy-m*sumx)/n;
    %goodness of ls method: chi square test and R-squares (param a,c)
    % redchi = sum((yfit-(m.*xfit+q)).^2)/(length(xfit)-2); %reduced chi square, should ~1
    r2 = 1-sum((yfit-(m.*xfit+q)).^2)/sum((yfit-avy).^2);
    % %print reduced chi-square for goodness
    % fprintf('reduced chi-square (param a,c): %f\n',redchi);
    % fprintf('r-squared (param a,c): %f\n\n',r2);
    %---------------------------------------------------------
    %least squares coefficients (param c)
    qbis = avy +sumx/n;
    %goodness of ls method: chi square test and R-squares (param a,c)
    % redchibis = sum((yfit-(-xfit+qbis)).^2)/(length(xfit)-2); %reduced chi square, should ~1
    r2bis = 1-sum((yfit-(-xfit+qbis)).^2)/sum((yfit-avy).^2);
    % %print reduced chi-square for goodness
    % fprintf('reduced chi-square (param c): %f\n',redchibis);
    % fprintf('r-squared (param c): %f\n\n',r2bis);
    %---------------------------------------------------------
    %plot log fit
    if(ifLin=='y')
        %setting
        labmq = sprintf("fit: y = %.3fx %+.1f",m,q);
        labq  = sprintf("fit: y = -x %+.1f",q);
        %figure
        figure()
        hold on
        plot(xfit, m.*xfit+q,'-b','LineWidth',1.5);               %,'DisplayName',labmq);
        %plot(xfit, -xfit+qbis,'--b','LineWidth',1.5);            %,'DisplayName',labq);
        plot(xfit,yfit,mark,'color','k','MarkerSize',markSize); %,'DisplayName','data, n = 25');
        plot(log(rmax(nmax-rejend+1:nmax)),log(fmax(nmax-rejend+1:nmax)), ...
            mark,'color','k','MarkerSize',markSize);            %,'DisplayName','rejected data')
        xl=xlim;
        yl=ylim;
        xpos=xl(1)+(xl(2)-xl(1))*0.67;
        ypos1=yl(1)+(yl(2)-yl(1))*0.93;
        ypos2=yl(1)+(yl(2)-yl(1))*0.88;
        text(xpos,ypos1,"n="+k);
        text(xpos,ypos2,labmq);
        % legend('location','northeast');
        xlabel('x = ln(r)')
        ylabel('y = ln(f)')
        %xticks([5 7 9 11 13 15])
        %save linear figure
        if savefigLin=='y'
            namefile = outputPath+"/fitPeaks_log";
            if(k<10)
              saveas(gcf,namefile+"_k00"+k+".png");
            elseif k<100
              saveas(gcf,namefile+"_k0"+k+".png");
            else
              saveas(gcf,namefile+"_k"+k+".png");
            end
        end
    end
    %---------------------------------------------------------------
    %back to pre-log coefficients and curve
    c=exp(q);
    a=m;
    cbis=exp(qbis);
    fit=c.*(r+eps*(r==0)).^a;
    fitbis=cbis./(r+eps*(r==0));
    %plot
    if (ifPeak=='y')
        figure()
        plot(r,absf,'-','LineWidth',1.5)
        hold on 
        plot(r(2:length(r)),fit(2:length(r)),'-r','LineWidth',1.5)
        %plot(r(2:length(r)),fitbis(2:length(r)),'--r','LineWidth',1.5)
        plot(rmax,fmax,'.k')
        xlabel('r')
        ylabel('f(r)')
        xlim([0,rmax(nmax)*1.1])
        ylim([0,fmax(1)*1.2])
        labca = sprintf("fit: y = %.1er^{%+.3f}",c,a);
        xl=xlim;
        yl=ylim;
        xpos=xl(1)+(xl(2)-xl(1))*0.67;
        ypos1=yl(1)+(yl(2)-yl(1))*0.93;
        ypos2=yl(1)+(yl(2)-yl(1))*0.88;
        text(xpos,ypos1,"n="+k);
        text(xpos,ypos2,labca);
        %save figure
        if savefigPeak=='y'
            namefile = outputPath+"/fitPeaks";
            if(k<10)
              saveas(gcf,namefile+"_k0"+k+".png");
            else
              saveas(gcf,namefile+"_k"+k+".png");
            end
        end
    end
    %---------------------------------------------------
end