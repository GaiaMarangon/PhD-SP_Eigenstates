clear all, close all, clc

set(groot,'defaultAxesFontSize',14)     % figures font size
set(groot,'DefaultTextFontSize',14)     % figures font size

%--- PARAMETERS ------------------------------------------------------
%node number
k=8;
%settings for data loading
inputPath = "../TestSensibility/results_NormOne/";
nametest = "results";
outputPath = "Figures_NormOne";

%--- LOADING DATA -------------------------------------------------------
%create folder for figure (if needed)
createSubfolder(outputPath);
%name setting
if k<10
    namek = "_k00";
elseif k<100
    namek = "_k0";
else
    namek = "_k";
end
%data loading
load(inputPath+nametest+namek+k+".mat",'r','f','epsPhi','phi','eigval');

%select data to be analyzed 
if iscell(f)
    % select set of data
    idx = size(r,1);
    r = r{idx};
    f = f{idx};
    phi = phi{idx};
    eigval = eigval{idx};
end

%--- ELABORATING RESULTS ------------------------------------------------

%%
%print eigenvalue
fprintf("eigenvalue: %.12e\n", eigval)

%%
%compute mass(r) and velocity(r), plots velocity(r)
iffig = 'y';
savefig = 'y';
[mass,vc]=massVelCurv(r,f,epsPhi,k,iffig,savefig,outputPath);

%%
%plots f(r), phi(r), log(f(r))
savefig_f = 'y';
savefig_phi = 'y';
savefig_logf ='y';
plot1dRes(r,f,phi,k,savefig_f,savefig_phi,savefig_logf,outputPath);

%%
%plot slices
if k<30 %can be extended (nr of points along r may be increased)
    ftype = "square"; %available: "abs" (for |f|), "square" (for f^2)
    viewtype = "3d"; %available: "3d" or "2d"
    savefig_f = 'y';
    plot3dRes_slices(r,f,k,ftype,viewtype,savefig_f,outputPath);
end

%%
%fit peaks and shows plot (original scale and loglog)
if k>4 %below, too few peaks to fit
    iffig = 'y';
    iffig_log = 'y';
    savefig = 'y';
    savefig_log = 'y';
    rejend = round(k*0.2); %nr of points to be rejected at the end
    [rmax,fmax,m,q,r2,qbis,r2bis,c,cbis]=fitPeaks(r,abs(f),k,epsPhi,iffig,iffig_log,savefig,savefig_log,rejend,outputPath);
end

%%

%fit velocity curve
if k>4 %below, too few peaks to fit
    iffig = 'y';
    savefig = 'y';
    rejbeg = round(k*0.05); %nr of points to be rejected at the beginning
    rejend = round(k*0.05); %nr of points to be rejected at the end
    [vc,ia,ib,m,q,r2]=fitVc(r,f,epsPhi,k,iffig,savefig,rejbeg,rejend,outputPath);
end










