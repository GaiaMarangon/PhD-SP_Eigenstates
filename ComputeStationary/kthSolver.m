function [e1,r,u1,phi1,f] = kthSolver(k,hin,epsPhi,tolInt,tolExt,maxiterInt,maxiterExt,fileID)
    %input data:
    %   k = nodes for desired eigenvector (Lieb: k=0)
    % -domain:
    %   h = grid stepsize
    %   dR = radius increment
    % -tolerances:
    %   epsPhi = for regularizing 1/r at r=0
    %   tolInt = for finding const eigval, in internal iter
    %   tolExt = for finding const eigval, in external iter
    %   maxiterInt maximum number of internal iterations
    %   maxiterExt maximum number of external iterations

    %output data:
    %   eigval = eigenvalue epsilon, = -1/2 * computed eigenvalue
    %   r = final domain (may vary in extension)
    %   u1 = eigenvector (solution for matter field after change of variable)
    %   phi1 = solution for potential
    %   f = solution for matter field
    
    %initial parameters:
    flagPrec=1;
    R=0;        %initial radius
    eout0 = 1;  %initial (fake) eigenvalue
    eout1 = 0;  %initial (fake) eigenvalue
    rold = linspace(0,5,50)';
    aConst = (5*pi)^(-3/4);
    uold = aConst .*rold.* exp(-rold.^2./10);
    %counter
    iterExt = 0;

    %outer loop: extending the domain
    while ( and(abs(eout0-eout1)>tolExt, iterExt<maxiterExt) )
        %select threshold for fine resolution, based on k
        if k==0
            thr = 4.0;
        elseif k==1
            thr = 2.5;
        elseif k==2
            thr = 2.3;
        elseif k==4
            thr = 2.0;
        elseif k<12
            thr = 1.8;
        else 
            thr = 1.2;
        end
        %select h, dR, based on fine or rough resolution
        if  R > (10*k^2+17*k+8)*thr
            % fine resolution
            %----------
            if (flagPrec) %solo prima volta
                tolInt = tolExt;
                fprintf('Applying fine resolution\n')
                fprintf(fileID,'Applying fine resolution\n');
            end
            flagPrec=0; %once in finer, stays in finer
            %----------
            dR = (10*k^2+17*k+8)/24;
            h = hin;
        else
            % rough resolution
            dR = (10*k^2+17*k+8)/8;
            h = (R+dR)/max(500,(k+1)*5); 
        end 
        iterExt = iterExt+1;
        %rename old variables
        eout0 = eout1;
        %extend domain;
        R=R+dR;
        %h= R/npoints; 
        r = (0:h:R)';
        ngrid = size(r,1);
        %start solving the problem with fixed radius
        %initial guess for u(r):
        u1 = interp1(rold,uold,r,'linear',0);
        %initial (fake) eigenvalues
        e0 = 1;
        e1 = 0;
        %fixed quantities, not to be updated
        beta = ones(ngrid-3,1)/(h^2);   
        alpha0 = -2*ones(ngrid-2,1)/(h^2);
        %counter
        iterInt=0;
        %solve for fixed domain [0,R]
        while ( and(abs(e1-e0)>tolInt, iterInt<maxiterInt) )
            iterInt = iterInt+1;
            %rename previous results
            u0 = u1;
            e0 = e1;
            %update phi: solve Poisson with previous iterations
            phi1 = solvePoisson(epsPhi,r,u0);
            %update u: solve eig problem at k-th node eigvect
            alpha = alpha0 -2*phi1(2:length(phi1)-1);
            A = diag(alpha) + diag(beta,1) + diag(beta,-1);
            [V,D] = eig(A);
            [e,indexV] = sort(diag(D),'descend');
            V= V(:,indexV);
            e1=e(k+1);
            u1=[0; V(:,k+1); 0];
            %adjust normalization
            u1 = u1 ./ sqrt( CTintegrate(1,length(u1),u1.^2,h,'f','b') );
            %fprintf('Res internal: %d\n',abs(e1-e0))
        end
        if (iterInt == maxiterInt)
            fprintf('Error: reaching maxiterInt\n')
            fprintf(fileID,'Error: reaching maxiterInt\n');
        else
        fprintf('iterInt: %d\t\teigval: %.10e\tR: %f\n',iterInt,abs(eout1-e1),R)
        fprintf(fileID,'iterInt: %d\t\teigval: %.10e\tR: %f\n',iterInt,abs(eout1-e1),R);
        end
        %update the eigenvalue to extend domain (outer loop)
        eout1=e1;
        rold=r;
        uold=u1;
    end
    if (iterExt == maxiterExt)
        fprintf('Error: reaching maxiterExt\n')
        fprintf(fileID,'Error: reaching maxiterExt\n');
    else
        fprintf('iterExt: %d\n\n',iterExt)
        fprintf(fileID,'iterExt: %d\n\n',iterExt);
    end
    
    %compute f
    f= u1./ (r+epsPhi*(r==0)) ;
    f(1) = f(2); %assuming zero derivative at r=0 (otherwise=0)
    %adjust sign
    if f(1)<0
        f=-f;
    end

end