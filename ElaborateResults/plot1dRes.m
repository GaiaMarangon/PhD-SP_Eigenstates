function []=plot1dRes(r,f,phi,k,savef,savephi,savelog,outputPath)
    %---print 1d results----------------------------------
    % f function
    figure()
    plot(r,f,'-','LineWidth',1.5)
    xlabel('r')
    ylabel('f(r)')
    if (savef=='y')
        namefile = outputPath+"/fFunc";
        if(k<10)
          saveas(gcf,namefile+"_k00"+k+".png");
        elseif k<100
          saveas(gcf,namefile+"_k0"+k+".png");
        else
          saveas(gcf,namefile+"_k"+k+".png");
        end
    end
    % phi function
    figure()
    plot(r,phi,'-','LineWidth',1.5)
    xlabel('r')
    ylabel('\phi(r)')
    if (savephi=='y')
        namefile = outputPath+"/phiFunc";
        if(k<10)
          saveas(gcf,namefile+"_k00"+k+".png");
        elseif k<100
          saveas(gcf,namefile+"_k0"+k+".png");
        else
          saveas(gcf,namefile+"_k"+k+".png");
        end
    end
    % log f function
    figure()
    ngrid = length(r);
    axis equal
    loglog(r(2:ngrid),abs(f(2:ngrid)),'-','LineWidth',1.5)
    xlabel('r')
    ylabel('f(r)')
    if (savelog=='y')
        namefile = outputPath+"/logfFunc";
        if(k<10)
          saveas(gcf,namefile+"_k00"+k+".png");
        elseif k<100
          saveas(gcf,namefile+"_k0"+k+".png");
        else
          saveas(gcf,namefile+"_k"+k+".png");
        end
    end
    %----------------------------------------------------
end