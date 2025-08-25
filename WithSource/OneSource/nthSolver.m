function [e1,r,phi1,f] = nthSolver(n,hin,phiParams,extRhoNames,epsPhi,tolInt,tolExt,maxiterInt,maxiterExt,fileID)
    %--- DESCRIPTION ----------------------------------------------------
    %input data:
    %   n = nodes for desired eigenvector 
    % -domain:
    %   h = grid stepsize
    % -source:
    %   phiParams = array Nx2, with each row = [r0,a] is a source 
    % -tolerances:
    %   epsPhi = for regularizing 1/r at r=0
    %   tolInt = for finding const eigval, in internal iter
    %   tolExt = for finding const eigval, in external iter
    %   maxiterInt maximum number of internal iterations
    %   maxiterExt maximum number of external iterations

    %output data:
    %   e1 = eigenvalue 
    %   r = final domain (may vary in extension)
    %   phi1 = solution for potential
    %   f = solution for matter field

    %NOTE: 
    % input parameters (phiParams) need to be in 4pi-Norm normalization
    % the output is in 4pi-Norm normalization
    %---------------------------------------------------------------------
    
    %--- PARAMETERS ------------------------------------------------------
    flagPrec=1;
    R=0;        %initial radius
    eout0 = 1;  %initial (fane) eigenvalue
    eout1 = 0;  %initial (fane) eigenvalue
    rOld = linspace(0,5,50)';
    aConst = (5*pi)^(-3/4);
    gOld = aConst .*rOld.* exp(-rOld.^2./10);
    %counter
    iterExt = 0;

    %--- OUTER LOOP, EXTENDING THE DOMAIN --------------------------------
    while ( and(abs(eout0-eout1)>tolExt, iterExt<maxiterExt) )
        %select threshold for fine resolution, based on n
        if n==0
            thr = 4.0;
        elseif n==1
            thr = 2.5;
        elseif n==2
            thr = 2.3;
        elseif n==4
            thr = 2.0;
        elseif n<12
            thr = 1.8;
        else 
            thr = 1.2;
        end
        %select h, dR, based on fine or rough resolution
        if  R > (10*n^2+17*n+8)*thr
            % fine resolution
            %----------
            if (flagPrec) %solo prima volta
                tolInt = tolExt;
                fprintf('Applying fine resolution\n')
                fprintf(fileID,'Applying fine resolution\n');
            end
            flagPrec=0; %once in finer, stays in finer
            %----------
            dR = (10*n^2+17*n+8)/24;
            h = hin;
        else
            % rough resolution
            dR = (10*n^2+17*n+8)/8;
            h = (R+dR)/max(500,(n+1)*5); 
        end 
        iterExt = iterExt+1;
        %rename old variables
        eout0 = eout1;
        %extend domain;
        R=R+dR;
        %h= R/npoints; 
        r = (0:h:R)';
        ngrid = size(r,1);
        %--- SOLVING PROBLEM WITH FIXED RADIUS ---------------------------
        %initial guess for g(r):
        g1 = interp1(rOld,gOld,r,'linear',0);
        %initial (fake) eigenvalues
        e0 = 1;
        e1 = 0;
        %fixed quantities, not to be updated
        beta = ones(ngrid-3,1)/(h^2);   
        alpha0 = -2*ones(ngrid-2,1)/(h^2);
        phiExt = phiExtSource(r,phiParams,extRhoNames,epsPhi); 
        %counter
        iterInt=0;
        %solve for fixed domain [0,R]
        while ( and(abs(e1-e0)>tolInt, iterInt<maxiterInt) )
            iterInt = iterInt+1;
            %rename previous results
            g0 = g1;
            e0 = e1;
            %update phi: solve Poisson with previous iterations
            phi1 = phiGenericRho(r, (g0./(r+epsPhi*(r==0))).^2 ,epsPhi) + phiExt;
            %update g: solve eig problem at n-th node eigvect
            alpha = alpha0 -2*phi1(2:length(phi1)-1);
            A = diag(alpha) + diag(beta,1) + diag(beta,-1);
            [V,D] = eig(A);
            [e,indexV] = sort(diag(D),'descend');
            V= V(:,indexV);
            e1=e(n+1);
            g1=[0; V(:,n+1); 0];
            %adjust normalization
            g1 = g1 ./ sqrt( CTintegrate(1,length(g1),g1.^2,h,'f','b') );
            %fprintf('Res internal: %d\n',abs(e1-e0))
        end
        if (iterInt == maxiterInt)
            fprintf('Error: reaching maxiterInt\n')
            fprintf(fileID,'Error: reaching maxiterInt\n');
        else
        fprintf('iterInt: %d\t\teigval: %.10e\tR: %f\n',iterInt,abs(eout1-e1),R)
        fprintf(fileID,'iterInt: %d\t\teigval: %.10e\tR: %f\n',iterInt,abs(eout1-e1),R);
        end
        %--- END OF FIXED RADIUS SOLVER ----------------------------------

        %update the eigenvalue to extend domain (outer loop)
        eout1=e1;
        rOld=r;
        gOld=g1;
    end
    if (iterExt == maxiterExt)
        fprintf('Error: reaching maxiterExt\n')
        fprintf(fileID,'Error: reaching maxiterExt\n');
    else
        fprintf('iterExt: %d\n\n',iterExt)
        fprintf(fileID,'iterExt: %d\n\n',iterExt);
    end
    %--- END OF OUTER LOOP; EXTENDING THE DOMAIN ------------------------



    %--- GO BACK TO DESIRED OUTPUT NOTATION ------------------------------
    % compute f from g
    f= g1./( (r+epsPhi*(r==0)));
    %adjust central value (assuming zero derivative at r=0)
    f(1) = f(2); 
    %adjust sign
    if f(1)<0
        f=-f;
    end

end