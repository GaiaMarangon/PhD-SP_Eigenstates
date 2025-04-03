function [vc,ia,ib,m,q,r2]=fitVc(r,f,epsPhi,k,iffig,savefig,rejbeg,rejend,outputPath)
    %-----------------------------------------------
    %get vc
    [~,vc]=massVelCurv(r,f,epsPhi,k,'n','n');
    %find peaks (local maxima) of vc(r)
    [~,~,imax]=findLocMaxima(r,vc,k);
    imax(1)=[];
    %set data for (linear) fit
    ia=imax(rejbeg+1);   
    ib=imax(length(imax)-rejend);
    xfit=r(ia:ib);
    yfit=vc(ia:ib);
    nfit=length(xfit);
    %compute param with (ordinary) least square
    sumx = sum(xfit);
    sumy = sum(yfit);
    m = (nfit*(sum(xfit.*yfit))-sumx*sumy)/(nfit*sum(xfit.^2)-sumx^2);
    q = (sumy-m*sumx)/nfit;
    %check goodness of fit (reduced chi square,coeff of determin)
    redchi  = sum((yfit-(m.*xfit+q)).^2)/(nfit-2);
    r2 = 1-sum((yfit-(m.*xfit+q)).^2)/sum((yfit-sumy/nfit).^2);
    %plot and print info
    if iffig=='y'
        %figure setting
        textFit=sprintf("y = (%.3e)x %+.1e",m,q);
        xlimVal="auto"; %[0 250000]
        legpos='southeast';
        %figure 
        figure()
        hold on
        L = legend('Location',legpos);
        L.AutoUpdate = 'off';
        plot(r,vc,'-k','LineWidth',1.5)
        L.AutoUpdate = 'on';
        plot(xfit,m.*xfit+q,'-b','DisplayName',textFit,'LineWidth',1.5) 
        xlabel('r')
        ylabel('v(r)')
        xlim(xlimVal);
        %save figure
        if savefig=='y'
            namefile = outputPath+"/fitVelFunc";
            if(k<10)
              saveas(gcf,namefile+"_k00"+k+".png");
            elseif k<100
              saveas(gcf,namefile+"_k0"+k+".png");
            else
              saveas(gcf,namefile+"_k"+k+".png");
            end
        end
    end
    % %print fit parameters
    % fprintf('fit parameters: m=%e, q=%e\n',m,q)
    % fprintf('goodness parameters: redchi=%e, r2=%f\n',redchi,r2)
    % fprintf('fit domain: ra=%f, rb=%f\n\n',r(ia),r(ib))
end