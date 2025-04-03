clearvars, clc, close all

%--- PARAMETERS ----------------------------------------------------------
%nodes 
k = 80;  
%tolerances:
epsPhi = 1e-015;       %for regularizing 1/r at r=0
tolInt = 1e-08;        %for finding const eigval, in internal iter
tolExt = 1e-10;        %for finding const eigval, in external iter
maxiterInt = 100;    %maximum number of internal iterations
maxiterExt = 100;   %maximum number of external iterations
%saving parameters
pathfile = "results_NormOne/";  %output folder
nametest = "results";           %file name basis
%grid stepsizes, to be tested
expDom95 = 10*k^2+17*k+8; %approximate domain, based on nr of nodes
if k>19
    h = expDom95*[0.00015,0.00013,0.00011,0.0001,0.00009]; %k=400,350,300,250,200,150,100,90,80,70,60,50,40,30,25,20
elseif k>15
    h = expDom95*[0.00016,0.00014,0.00012,0.00010,0.00009];  %k=19,18,17,16
elseif k>12 
    h = expDom95*[0.00018,0.00015,0.00013,0.00012,0.00011];  %k=15,14,13
elseif k>8
    h = expDom95*[0.00021,0.00018,0.00015,0.00013,0.00012];  %k=12,11,10,9
elseif k>4
    h = expDom95*[0.00025,0.00021,0.00018,0.00015,0.00013];  %k=5,6,7,8
elseif k>3
    h = expDom95*[0.00026,0.00022,0.00019,0.00016,0.00014];  %k=4
elseif k>2
    h = expDom95*[0.00028,0.00024,0.00021,0.00018,0.00015];  %k=3
elseif k>1
    h = expDom95*[0.0003,0.00026,0.00023,0.00020,0.00017];  %k=2
elseif k>0
    h = expDom95*[0.00035,0.0003,0.00025,0.00021,0.00018];  %k=1
else
    h = expDom95*[0.0004,0.00035,0.0003,0.00025,0.00020];  %k=0
end

%--- INITIAL SETTING -----------------------------------------------------
%nr of sensibility tests, for each k
nSens = length(h);
%setting output
eigval = cell(nSens,1);
r = cell(nSens,1);
u1 = cell(nSens,1);
phi = cell(nSens,1);
f = cell(nSens,1);
elapsedTime = cell(nSens,1);
%create output folder, if not existing
createSubfolder(pathfile);
%name file setting
if k<10
    namek = "_k00";
elseif k<100
    namek = "_k0";
else
    namek = "_k";
end
%opening file for log
fileID = fopen(pathfile+nametest+"_log"+namek+k+".txt",'w');

%--- SOLVING PROBLEM -----------------------------------------------------
for i=1:nSens
    %solving the coupled problem
    tic
    [eigval{i},r{i},u1{i},phi{i},f{i}] = kthSolver(k,h(i),epsPhi,tolInt,tolExt,maxiterInt,maxiterExt,fileID);
    elapsedTime{i} = toc;
    %scale for 1 normalization (comment to stay with 4pi normalization)
    [eigval{i}, f{i}, phi{i}, r{i}, u1{i}]=Norm4piToNormOne(eigval{i}, f{i}, phi{i}, r{i}, u1{i});    
end

%--- CHECKS --------------------------------------------------------------
%visual comparison (plot function f(r) for last two tested cases)
figure()
hold on
plot(r{nSens-1},f{nSens-1},'-r')
plot(r{nSens},f{nSens},':k')

%--- SAVING RESULTS ------------------------------------------------------
%save data
save(pathfile+nametest+namek+k+".mat");
%close saving data
fclose(fileID);


