function [rmax,fmax,imax]=findLocMaxima(r,absf,k)
    nmax=k+1;
    rmax=zeros(nmax,1);
    fmax=zeros(nmax,1);
    imax=zeros(nmax,1); %index of maxima in r,f
    inc=0; %before, increasing? 0=no, 1=yes
    j=0;    %counting nr of maxima
    for i=1:length(r)
        if j==nmax
            break
        end
        if i==1
            fmax(1)=absf(1);
            rmax(1)=r(1);
            imax(1)=1;
            j=1;
            continue
        end
        if i==2
            continue    %qui = f(1) per costruzione mia, ignoro
        end
        if(absf(i)<absf(i-1)) %decreasing detected
            if inc==1 
                j=j+1;
                fmax(j)=absf(i-1);
                rmax(j)=r(i-1);
                imax(j)=i-1;
                inc=0;
            end
        else %increasing detected
            inc=1;
        end
    end
end


%--main to test findLocMaxima-----------
% r=linspace(0,4*pi)';
% absf=abs(cos(r));
% [rmax,fmax]=findLocMaxima(r,absf,3);
% plot(r,absf,'-')
% hold on
% plot(rmax,fmax,'or')
%---------------------------------------