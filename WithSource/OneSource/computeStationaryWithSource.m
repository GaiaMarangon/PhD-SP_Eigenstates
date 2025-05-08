function computeStationaryWithSource(r0,a,nameFile,...
    n,h,rhoName,epsPhi,tolInt,tolExt,maxiterInt,maxiterExt)
    % Given source parameters (r0,a) in 4pi-normalization
    % computes the n-nodes stationary state with that sources, and the relative info:
    % - eigenvalue, eigenfunction, potential, velocity;
    % - for each source (D,B,G): density, velocity;
    % - total velocity (stationary state + sources);
    % saves all info to file 1-normalization

    % compute stationary state
    phiParams = [r0,a];
    fileID = fopen(nameFile+"_nthSolver.txt",'w');
    [eigval,rDM_Num,phi_Num,f_Num] = nthSolver(n,h,phiParams,rhoName,epsPhi,tolInt,tolExt,maxiterInt,maxiterExt,fileID);
    % compute external density
    rhoExt = rhoExtSource(rDM_Num,phiParams,rhoName);
    %switch from 4pi-norm to 1-norm
    [eigval, f_Num, phi_Num, rDM_Num, rhoExt]=Norm4piToNormOne_extSource(eigval, f_Num, phi_Num, rDM_Num, rhoExt);   
    r0 = r0 * (4*pi);
    a  = a  / (4*pi)^4;
    % compute DM velocity curve
    [~,vDM_Num] = massVelCurv(rDM_Num,f_Num,epsPhi,n,'n','n');
    % compute source velocity
    [~,vSource_Num] = massVelCurv(rDM_Num,sqrt(rhoExt),epsPhi,n,'n','n');
    % compute combined velocity DM+source
    vTot_Num = sqrt(vDM_Num.^2 + vSource_Num.^2);
    % save 
    save(nameFile+".mat","n","r0","a",...
        "eigval","rDM_Num","phi_Num","f_Num","vDM_Num",...
        "rhoExt","vSource_Num","vTot_Num");
end