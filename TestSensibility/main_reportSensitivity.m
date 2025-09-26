clearvars, clc,close all

%--- PARAMETERS ----------------------------------------------------------
%nodes
k=5;
%loading and saving parameters
pathfile = "results_NormOne/";  %output folder
nametest = "results";           %file name basis

%--- GRAPHICS
colorVec_default = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];
% figure settings
set(groot,'defaultAxesFontSize',20)     % figures font size
set(groot,'DefaultTextFontSize',20)     % figures font size
set(groot, 'defaultAxesTickLabelInterpreter','latex');  % latex labels
set(groot, 'defaultLegendInterpreter','latex');         % latex labels
set(groot, 'defaultTextInterpreter','latex');           % latex labels

%--- DATA LOADING ------------------------------------------------
%setting file name
if k<10
    namek = "_k00";
elseif k<100
    namek = "_k0";
else
    namek = "_k";
end
%loading data
load(pathfile+nametest+namek+k+".mat",'nSens','eigval','r','f','u1','epsPhi')
%opening file for report 
fileID = fopen(pathfile+nametest+"_reportSens"+namek+k+".txt",'w');

% %--- OPTIONAL (PLOTS)--------------------
% %plot for reference - f(r)
% figure()
% plot(r{nSens},f{nSens})
% %plot for reference - v(r)
% [mass,vel]=massVelCurv(r{nSens},f{nSens},epsPhi,k,'n','n');
% figure()
% plot(r{nSens},vel)


%%

%--- eigval ------------------
fprintf(fileID,'eigenvalues:\n');
figure()
hold on
%nSens = 1;
for i=1:nSens
    h = r{i}(end) - r{i}(end-1);
    fprintf(fileID,'\th: %f\teigval: %.10e\n',h,eigval{i});
    scatter(h,eigval{i},[],colorVec_default(1,:),'o','filled')
end
xlabel('h')
ylabel('eigval')
hold off

%%
%--- relative error eigenval --------------------
fprintf(fileID,'\nRelative error on eigenvalues:\n');
fprintf(fileID,'\th\tRelDiff\n');
figure()
hold on
for i=2:nSens
    h = r{i}(end) - r{i}(end-1);
    relErr = (eigval{i}-eigval{i-1}) / eigval{nSens};
    fprintf(fileID,'\t%f\t%e\n',h,relErr);
    scatter(h,relErr,[],colorVec_default(1,:),'o','filled')
end
xlabel('h')
ylabel('$\Delta \varepsilon_n$')
hold off
%

%%
%--- Outermost local extremum f(r) -----------
fprintf(fileID,'\nOutermost local extremum |f|(r):\n');
fprintf(fileID,'\th\tr_Out\t|f_Out|\n');

rmax_vec = zeros(nSens-1,1);
fmax_vec = zeros(nSens-1,1);
for i=1:nSens
    h = r{i}(end) - r{i}(end-1);
    [rmax,fmax,~]=findLocMaxima(r{i},abs(f{i}),k);
    [rmax_vec(i),idxmax] = max(rmax);
    fmax_vec(i) = fmax(idxmax);
    fprintf(fileID,'\t%f\t%e\t%e\n', h,rmax_vec(i),fmax_vec(i));
end

fprintf(fileID,'\n\th\tRelDiff_r\tRelDiff_|f|\n');
figure(10)
hold on
figure(11)
hold on
for i=2:nSens
        h = r{i}(end) - r{i}(end-1);
        relErr_r = (rmax_vec(i)-rmax_vec(i-1)) / rmax_vec(end);
        relErr_f = (fmax_vec(i)-fmax_vec(i-1)) / fmax_vec(end);
        fprintf(fileID,'\t%f\t%e\t%e\n', h,relErr_r,relErr_f);
        figure(10)
        scatter(h,relErr_r,[],colorVec_default(1,:),'o','filled')
        figure(11)
        scatter(h,relErr_f,[],colorVec_default(1,:),'o','filled')
end
figure(10)
xlabel('h')
ylabel('$\Delta r_{Out}$')
hold off
figure(11)
xlabel('h')
ylabel('$\Delta f_{Out}$')
hold off

%%
%--- Outermost local extremum v(r) -----------
fprintf(fileID,'\nOutermost local extremum v(r):\n');
fprintf(fileID,'\th\tr_Out\tv_Out\n');

rmax_vec = zeros(nSens-1,1);
vmax_vec = zeros(nSens-1,1);
for i=1:nSens
    h = r{i}(end) - r{i}(end-1);
    [mass,vel]=massVelCurv(r{i},f{i},epsPhi,k,'n','n');
    [rmax,vmax,~]=findLocMaxima(r{i},vel,2*k+1);
    [rmax_vec(i),idxmax] = max(rmax);
    vmax_vec(i) = vmax(idxmax);
    fprintf(fileID,'\t%f\t%e\t%e\n', h,rmax_vec(i),vmax_vec(i));
end

fprintf(fileID,'\n\th\tRelDiff_r\tRelDiff_v\n');
figure(20)
hold on
figure(21)
hold on
for i=2:nSens
        h = r{i}(end) - r{i}(end-1);
        relErr_r = (rmax_vec(i)-rmax_vec(i-1)) / rmax_vec(end);
        relErr_v = (vmax_vec(i)-vmax_vec(i-1)) / vmax_vec(end);
        fprintf(fileID,'\t%f\t%e\t%e\n', h,relErr_r,relErr_v);
        figure(20)
        scatter(h,relErr_r,[],colorVec_default(1,:),'o','filled')
        figure(21)
        scatter(h,relErr_v,[],colorVec_default(1,:),'o','filled')
end
figure(20)
xlabel('h')
ylabel('$\Delta r_{Out}$')
hold off
figure(21)
xlabel('h')
ylabel('$\Delta v_{Out}$')
hold off

%%
%--- Innermost local extremum v(r) -----------
fprintf(fileID,'\nInnermost local extremum v(r):\n');
fprintf(fileID,'\th\tr_Inn\tv_Inn\n');

rmin_vec = zeros(nSens-1,1);
vmin_vec = zeros(nSens-1,1);
for i=1:nSens
    h = r{i}(end) - r{i}(end-1);
    [mass,vel]=massVelCurv(r{i},f{i},epsPhi,k,'n','n');
    [rmax,vmax,~]=findLocMaxima(r{i},vel,k+1);
    [rmin_vec(i),idxmin] = min(rmax(2:end));
    vmin_vec(i) = vmax(idxmin+1);
    fprintf(fileID,'\t%f\t%e\t%e\n', h,rmin_vec(i),vmin_vec(i));
end

fprintf(fileID,'\n\th\tRelDiff_r\tRelDiff_v\n');
figure(30)
hold on
figure(31)
hold on
for i=2:nSens
        h = r{i}(end) - r{i}(end-1);
        relErr_r = (rmin_vec(i)-rmin_vec(i-1)) / rmin_vec(end);
        relErr_v = (vmin_vec(i)-vmin_vec(i-1)) / vmin_vec(end);
        fprintf(fileID,'\t%f\t%e\t%e\n', h,relErr_r,relErr_v);
        figure(30)
        scatter(h,relErr_r,[],colorVec_default(1,:),'o','filled')
        figure(31)
        scatter(h,relErr_v,[],colorVec_default(1,:),'o','filled')
end
figure(30)
xlabel('h')
ylabel('$\Delta r_{inn}$')
hold off
figure(31)
xlabel('h')
ylabel('$\Delta v_{inn}$')
hold off

%%
%--- maxima position -----------
fprintf(fileID,'\nRelative Difference of Local Maxima of Eigenfunction (average):\n');
fprintf(fileID,'\th\tAvgDiff_r\tAvgDiff_|f|\n');
figure(40)
hold on
figure(41)
hold on

rmax = cell(nSens,1);
fmax = cell(nSens,1);
for i=1:nSens
    [rmax{i},fmax{i},~]=findLocMaxima(r{i},abs(f{i}),k);
    if i>1
        h = r{i}(end) - r{i}(end-1);
        h_prev = r{i-1}(end) - r{i-1}(end-1);
        avgRel_r = sum( abs(rmax{i}-rmax{i-1}) ) / length(rmax{i})  /rmax{i}(end);
        avgRel_f = sum( abs(fmax{i}-fmax{i-1}) ) / length(fmax{i})  /f{i}(1);
        fprintf(fileID,'\t%f\t%e\t%e\n', h,avgRel_r, avgRel_f);
        figure(40)
        scatter(h,avgRel_r,[],colorVec_default(1,:),'o','filled')
        figure(41)
        scatter(h,avgRel_f,[],colorVec_default(1,:),'o','filled')
    end
end
figure(40)
xlabel('h')
ylabel('$\frac{\sum_j |r_{max,j}(h_i)-r_{max,j}(h_{i-1})|}{ n_{max} }  \frac{1}{ r_{Out}}$')
hold off
figure(41)
xlabel('h')
ylabel('$\frac{\sum_j |f_{max,j}(h_i)-f_{max,j}(h_{i-1})|}{ n_{max}  } \frac{1}{f_{0}}$')
hold off


%%
%--- elapsed time --------
% fprintf(fileID,'\nelapsed times:\n');
% figure(5)
% hold on
% for i=1:nSens
%     h = r{i}(end) - r{i}(end-1);
%     fprintf(fileID,'\th: %d\ttime (s): %.3f\n',h,elapsedTime{i});
%     plot(h,elapsedTime{i},'o')
% end
% xlabel('h')
% ylabel('elapsed time (s)')
% hold off

%%

fclose(fileID);
%-------------------------------------------------------------------------
