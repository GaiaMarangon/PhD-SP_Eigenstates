clearvars, clc, close all

%--- PARAMETERS ---------------------------------------------------
inputFolder  = "../TestSensibility/results_NormOne"; % "../ComputeStationary/results_NormOne";
outputFolder = pwd+"\NormOne";                   %created, if not existing
textComment = "normalized to 1 in radial coord"; %included in header of output files

%--- PRINTS -------------------------------------------------------
%create subfolders for data, if needed
createSubfolder(outputFolder);
createSubfolder(outputFolder+"\rfNumData");
createSubfolder(outputFolder+"\rmNumData");
createSubfolder(outputFolder+"\rvNumData");
createSubfolder(outputFolder+"\rpNumData");
createSubfolder(outputFolder+"\neNumData");
%run through files with results and:
files=dir(fullfile(inputFolder,'*.mat'));
eigPrint = zeros(numel(files),1);
nPrint   = zeros(numel(files),1);
for i=1:numel(files)
    %load data from file 
    loadName=string(files(i).folder) + "\" +string(files(i).name);
    load(loadName,"k","r","f","phi","eigval","epsPhi");

    %extract number of nodes k for that file, with zeros
    kstr=extractBetween(files(i).name,strlength("results_k")+1,strlength("results_k")+3);
    %set name for file to write to
    writeName_rf= outputFolder+"\rfNumData\rfNumData_n"+kstr+".dat"; %f(r)
    writeName_rm= outputFolder+"\rmNumData\rmNumData_n"+kstr+".dat"; %mass(r)
    writeName_rv= outputFolder+"\rvNumData\rvNumData_n"+kstr+".dat"; %vel(r)
    writeName_rp= outputFolder+"\rpNumData\rpNumData_n"+kstr+".dat"; %phi(r)
    writeName_ne= outputFolder+"\neNumData\neNumData.dat";        %eigval(n)

    %extract variables to print, according to the type of results 
    if iscell(eigval)
        nPrint(i)   = k;
        eigPrint(i) = eigval{length(eigval)};
        %
        rPrint   = r{length(r)};
        fPrint   = f{length(f)};
        phiPrint = phi{length(phi)};
        %
        [mass,vel]=massVelCurv(r{length(r)},f{length(f)},epsPhi,k,'n','n');
        massPrint = mass;
        velPrint  = vel;
    else
        nPrint(i)   = k;
        eigPrint(i) = eigval;
        %
        rPrint   = r;
        fPrint   = f;
        phiPrint = phi;
        %
        [mass,vel]=massVelCurv(r,f,epsPhi,k,'n','n');
        massPrint = mass;
        velPrint  = vel;
    end
    
    %write results to file - f(r)
    fileID = fopen(writeName_rf,'w');
    fprintf(fileID,"#eigenfunctions f(r) from MATLAB code, "+textComment+"\n");
    fprintf(fileID,'#r=radius, f(r)= eigenfunction\n');
    fprintf(fileID,'#r\tf(r)\n');
    fclose(fileID);
    %
    writematrix([rPrint,fPrint], writeName_rf,'Delimiter',"\t",'WriteMode','append');

    %write results to file - m(r)
    fileID = fopen(writeName_rm,'w');
    fprintf(fileID,"#masses from MATLAB code, "+textComment+"\n");
    fprintf(fileID,'#r=radius, m(r)= mass\n');
    fprintf(fileID,'#r\tm(r)\n');
    fclose(fileID);
    %
    writematrix([rPrint,massPrint], writeName_rm,'Delimiter',"\t",'WriteMode','append');

    %write results to file - v(r)
    fileID = fopen(writeName_rv,'w');
    fprintf(fileID,"#velocities from MATLAB code, "+textComment+"\n");
    fprintf(fileID,'#r=radius, v(r)= velocity, rotation curve\n');
    fprintf(fileID,'#r\tv(r)\n');
    fclose(fileID);
    %
    writematrix([rPrint,velPrint], writeName_rv,'Delimiter',"\t",'WriteMode','append');

    %write results to file - phi(r)
    fileID = fopen(writeName_rp,'w');
    fprintf(fileID,"#potentials from MATLAB code"+textComment+"\n");
    fprintf(fileID,'#r=radius, phi(r)= gravitational potential\n');
    fprintf(fileID,'#r\tphi(r)\n');
    fclose(fileID);
    %
    writematrix([rPrint,phiPrint], writeName_rp,'Delimiter',"\t",'WriteMode','append');
end

%write results to file - phi(r)
fileID = fopen(writeName_ne,'w');
fprintf(fileID,"#eigenvalues from MATLAB code,"+textComment+"\n");
fprintf(fileID,'#n=nr nodes, eig(n) = eigenvalues\n');
fprintf(fileID,'#n\teig(n)\n');
fclose(fileID);
%
writematrix([nPrint,eigPrint], writeName_ne,'Delimiter',"\t",'WriteMode','append');
