function [mass,vc]=massVelCurv(r,f,eps,k,iffig,savefig)
    %compute mass and velocities
    h=r(2)-r(1);
    mass = zeros(length(r),1);
    for i=2:length(r)
        if i==length(r)
            rule2 = 'b';
        else
            rule2 = 'c';
        end
        mass(i) = CTintegrate(1,i,f.^2.*r.^2,h,'f',rule2);
    end
    vc = sqrt(mass./(r+eps*(r==0)));
    %--------------------------------------------------
    %plot vc figure
    if iffig=='y'
        figure()
        plot(r,vc,'-','LineWidth',1.5','DisplayName',"n="+k)
        hold on 
        xlabel('r')
        ylabel('v(r)')
        %xlim([0 400])
        legend('Location','southeast')
        %optional: save vc figure
        if (savefig=='y')
            namefile = "Figures/vFunc";
            if(k<10)
              saveas(gcf,namefile+"_k00"+k+".png");
            elseif k<100
              saveas(gcf,namefile+"_k0"+k+".png");
            else
              saveas(gcf,namefile+"_k"+k+".png");
            end
        end
    end
end