clearvars, clc, close all

%--- PARAMETERS ---------------------------------------------------------
% - nodes for desired eigenvector 
k = 11;  
% - tolerances:
epsPhi = 1e-015;        %for regularizing 1/r at r=0
tolInt = 1e-08;         %for finding const eigval, in internal iter
tolExt = 1e-10;         %for finding const eigval, in external iter
maxiterInt = 100;       %maximum number of internal iterations
maxiterExt = 100;       %maximum number of external iterations
% - grid stepsize, based on domain (rough, for solving time of about 60s)
expDom95 = 10*k^2+17*k+8; 
if k==0 
    h = expDom95*0.0012;
elseif k==1
    h = expDom95*0.0008;
elseif k==2
    h = expDom95*0.0007;
elseif k==3
    h = expDom95*0.0006;
elseif k==4
    h = expDom95*0.0006;
elseif k<11
    h = expDom95*0.00055;
elseif k<15
    h = expDom95*0.0005;
elseif k<40            
    h = expDom95*0.0004;
else    
    %for high nr of nodes, higher refinement with longer solving time is suggested
    h = expDom95*0.0003;
end
%settings for data saving
pathfile = "results_NormOne\";
nametest = "results" ;
if k<10
    namek = "_k00";
elseif k<100
    namek = "_k0";
else
    namek = "_k";
end

%--- SOLVING ------------------------------------------------------------
%create output folder, if not existing
createSubfolder(pathfile);
%opening file for saving logs
fileID = fopen(pathfile+nametest+"_log"+namek+k+".txt",'w');
%solving the coupled problem
tic
[eigval,r,u1,phi,f] = kthSolver(k,h,epsPhi,tolInt,tolExt,maxiterInt,maxiterExt,fileID);
elapsedTime = toc;
% print time to screen
fprintf("elapsed time: %.4f\n",elapsedTime)

%--- OPTIONAL --------------------------------------------------------
% print time to screen and plot f(r) for visual check
fprintf("elapsed time: %.4f\n",elapsedTime)
plot(r,f)
%compare with literature (Bernstein,Giladi,Jones-1998)
fB= f./ sqrt(4*pi) ;    %scale to match literature notation
eigvalB = -eigval/2;    %scale to match literature notation
fprintf("eigenvalue: %.12e\n", eigvalB)
[~,~,imax]=findLocMaxima(r,abs(fB),k);
roeMy = r(imax(length(imax)));
fprintf('roe: %e\n',roeMy);
%
%scale for 1 normalization (comment to stay with 4pi normalization)
[eigval, f, phi, r, u1]=Norm4piToNormOne(eigval, f, phi, r, u1);

%--- SAVING DATA --------------------------------------------------------
%save all results to mat file, for future loading
save(pathfile+nametest+namek+k+".mat");



