clear all, close all, clc

% DESCRIPTION:
% Given a density, performs a series of tests on the associated stationary state, varying the density parameters.
% - (computes) loads the numerical curve with no sources
% - plots the density and the potential for the reference parameters
% - varies r0, (computing) loading numerical curves; and plotting: (1)densities (2)eigenfunction (3)potential (4)velocity (5)total velocity;
% - varies a,  (computing) loading numerical curves; and plotting: (1)densities (2)eigenfunction (3)potential (4)velocity (5)total velocity;


%-- PARAMETERS ---------------------------------------------------

% density name
rhoName = "Exp";
% rhoName = "TruncatedPlummer";
% rhoName = "HardBall";

% varied parameters
rScaleVec = exp(linspace(log(0.1),log(0.3),6));
aScaleVec = linspace(0.03,0.13,6);
% fixed parameters
rScaleRef = 0.18;
aScaleRef = 0.03;
% number of nodes
n = 5;
% tolerances:
epsPhi = 1e-015;        %for regularizing 1/r at r=0
tolInt = 1e-08;         %for finding const eigval, in internal iter
tolExt = 1e-10;         %for finding const eigval, in external iter
maxiterInt = 100;       %maximum number of internal iterations
maxiterExt = 100;       %maximum number of external iterations
h = (10*n^2+17*n+8)*0.002; %discretization step
% saving parameters
outDir = "output/"+rhoName+"/";
%create (sub)folders if not present
createSubfolder(outDir+"/noSource")
createSubfolder(outDir+"/Figures")
createSubfolder(outDir+"/rVary")
createSubfolder(outDir+"/aVary")

% colors
colorVec_default = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];

% figgure settings
set(groot,'defaultAxesFontSize',20)     % figures font size
set(groot,'DefaultTextFontSize',20)     % figures font size
set(groot, 'defaultAxesTickLabelInterpreter','latex');  % latex labels
set(groot, 'defaultLegendInterpreter','latex');         % latex labels
set(groot, 'defaultTextInterpreter','latex');         % latex labels
%%
 
% %-- (COMPUTE) NO SOURCE --------------------------------------------------
% nameFile = outDir+"noSource/noSource";
% % compute and save stationary state
%     computeStationaryWithSource(1,0,nameFile,...
%     n,h,rhoName,epsPhi,tolInt,tolExt,maxiterInt,maxiterExt);

%%

%-- LOAD NO SOURCE -----------------------------------------------------------

% load data
temp = load(outDir+"noSource/noSource.mat", 'f_Num', 'phi_Num', 'rDM_Num', 'vDM_Num');
f_noSource   = temp.f_Num;
phi_noSource = temp.phi_Num;
r_noSource   = temp.rDM_Num;
vDM_noSource = temp.vDM_Num;
clear('temp')
% get r_noSource_out 
[rmax,~,~]=findLocMaxima(r_noSource,abs(f_noSource),n);
r_noSource_out = rmax(n+1);

%%

%-- PLOT DENSITY AND POTENTIAL FOR REFERENCE PARAMETERS ------------------
 
% Compute density and potential (expressed in norm one)
params = [rScaleRef*r_noSource_out, aScaleRef*f_noSource(1)^2, epsPhi];
rho = rhoExtSource(r_noSource,params,rhoName);
phi = phiOneSource(r_noSource,params,rhoName);

% plot density
figure()
plot(r_noSource,rho,'LineWidth',1.5)
xlabel("$r$")
ylabel("$\rho$")
saveas(gcf,outDir+"Figures/rhoOneRef_"+rhoName+".png");

% plot potential
figure()
plot(r_noSource,phi,'LineWidth',1.5)
xlabel("$r$")
ylabel("$\phi_{\rho}$")
saveas(gcf,outDir+"Figures/phiOneRef_"+rhoName+".png");


%%

% %--- (COMPUTE) R VARY ------------------------------------------------------
% % scales
% rScale = rScaleVec;
% aScale = aScaleRef;
% % fixed parameter
% a = f_noSource(1)^2 *aScale;
% 
% for i=1:length(rScale)
%     % varying parameter
%     r0 = r_noSource_out*rScale(i);
%     nameFile = outDir+"rVary/aScale"+string(round(aScale,2)*100)+"rScale"+string(round(rScale(i),2)*100);
%     % switch back to 4pi norm for solver computations
%     r0_4pi = r0 / (4*pi);
%     a_4pi  = a  * (4*pi)^4;
%     % compute and save stationary state
%     computeStationaryWithSource(r0_4pi,a_4pi,nameFile,...
%     n,h,rhoName,epsPhi,tolInt,tolExt,maxiterInt,maxiterExt);
% end

%%

close all

%-- PLOTS VARYING R ----------------------------------------------

% scales
rScale = rScaleVec;
aScale = aScaleRef;

% figure setting - square root of density
figure(1)
hLegend = legend('Location','northeast','Interpreter','latex');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$\sqrt{\rho}$")
hold on
plot(r_noSource,f_noSource,'-','Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName','$f_n$ no source')
% figure setting - eigenfunction
figure(2)
hLegend = legend('Location','northeast','Interpreter','latex');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
set(gcf,'Position',[100 100 560*1 420])
xlim([0,6000])
xlabel("$r$")
ylabel("$f_"+string(n)+"$")
hold on
plot(r_noSource,f_noSource,'-','Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName','no source')
% figure setting - potential
figure(3)
hLegend = legend('Location','southeast');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$\phi_"+string(n)+"$")
hold on
plot(r_noSource,phi_noSource,'-','Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName','no source')
% figure setting - velocity
figure(4)
hLegend = legend('Location','southeast','Interpreter','latex');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
set(gcf,'Position',[100 100 560*1 420])
xlim([0,6000])
xlabel("$r$")
ylabel("$v_"+string(n)+"$")
hold on
plot(r_noSource,vDM_noSource,'-','Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName','no source')
% figure setting - composing velocity
figure(5)
hLegend = legend('Location','northeastoutside','Interpreter','latex');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
set(gcf,'Position',[100 100 560*1.4 420])
xlabel("$r$")
ylabel("$v_{Tot}$")
xlim([0,6000])
hold on
plot(nan,nan,'-','Color',[0.7,0.7,0.7],'LineWidth',1.5,'DisplayName', "$v_{Tot}$")
plot(nan,nan,'--','Color',[0.7,0.7,0.7],'LineWidth',1.5,'DisplayName', "$v_{\rho}$")
plot(nan,nan,':','Color',[0.7,0.7,0.7],'LineWidth',1.5,'DisplayName', "$v_"+string(n)+"$")
scatter(nan,nan,'o','MarkerFaceColor',[0.7,0.7,0.7],'MarkerEdgeColor','none','DisplayName', "$r_{0}$")

% add plots as r varies
for i=1:length(rScale)
    % load computed data
    nameFile = outDir+"rVary/aScale"+string(round(aScale,2)*100)+"rScale"+string(round(rScale(i),2)*100);
    load(nameFile, "r0","a","rDM_Num","phi_Num","f_Num","vDM_Num","rhoExt","vSource_Num","vTot_Num")
    % common setting
    label = "$s_{r_0}$="+string(round(rScale(i),2));
    transp = i/(length(rScale));
    % add to figure - density
    figure(1)
    plot(rDM_Num,sqrt(rhoExt),'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,'DisplayName', label)
    % add to figure - eigenfunction
    figure(2)
    plot(rDM_Num,f_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,'DisplayName',label)
    % add to figure - potential
    figure(3)
    plot(rDM_Num,phi_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,'DisplayName',label)
    % add to figure - velocity
    figure(4)
    plot(rDM_Num,vDM_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,'DisplayName', label)
    % add to figure - composing velocity
    figure(5)
    plot(rDM_Num,vDM_Num,    ':', 'Color',[colorVec_default(1,:),transp],'LineWidth',1.5,'HandleVisibility','off');%,'DisplayName', label + ", $V_{DM}$")
    plot(rDM_Num,vSource_Num,'--','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,'HandleVisibility','off');%'DisplayName', label + ", $V_{Source}$")
    plot(rDM_Num,vTot_Num,   '-', 'Color',[colorVec_default(1,:),transp],'LineWidth',1.5,'DisplayName', label)
    scatter(r0,0,'o','MarkerFaceColor',colorVec_default(1,:),'MarkerFaceAlpha',transp,'MarkerEdgeColor','none','HandleVisibility','off');%'DisplayName', label + ", $r_{0}$")
end

% save figure - density
figure(1)
saveas(gcf,outDir+"Figures/varyR_density.png")
% save figure - eigenfunction
figure(2)
saveas(gcf,outDir+"Figures/varyR_eigenfunction.png")
% save figure - potential
figure(3)
saveas(gcf,outDir+"Figures/varyR_potential.png")
% save figure - velocity
figure(4)
saveas(gcf,outDir+"Figures/varyR_velocity.png")
% save figure - composing velocity
figure(5)
saveas(gcf,outDir+"Figures/varyR_composingVelocity.png")

%%

% %--- (COMPUTE) A VARY ------------------------------------------------------
% % scales
% rScale = rScaleRef;
% aScale = aScaleVec;
% % fixed parameter
% r0 = r_noSource_out*rScale;
% 
% for i=1:length(aScale)
%     % varying parameter
%     a = f_noSource(1)^2 *aScale(i);
%     nameFile = outDir+"aVary/aScale"+string(round(aScale(i),2)*100)+"rScale"+string(round(rScale,2)*100);
%     % switch back to 4pi norm for solver computations
%     r0_4pi = r0 / (4*pi);
%     a_4pi  = a  * (4*pi)^4;
%     % compute and save stationary state
%     computeStationaryWithSource(r0_4pi,a_4pi,nameFile,...
%     n,h,rhoName,epsPhi,tolInt,tolExt,maxiterInt,maxiterExt);
% end

%%

close all

%-- PLOTS VARYING A ----------------------------------------------

% scales
rScale = rScaleRef;
aScale = aScaleVec;

% figure setting - density
figure(1)
hLegend = legend('Location','northeast','Interpreter','latex');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$\sqrt{\rho}$")
hold on
plot(r_noSource,f_noSource,'-','Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName','$f_n$ no source')
% figure setting - eigenfunction
figure(2)
hLegend = legend('Location','northeast','Interpreter','latex');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
set(gcf,'Position',[100 100 560*1 420])
xlim([0,6000])
xlabel("$r$")
ylabel("$f_"+string(n)+"$")
hold on
plot(r_noSource,f_noSource,'-','Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName','no source')
% figure setting - potential
figure(3)
hLegend = legend('Location','southeast','Interpreter','latex');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$\phi_"+string(n)+"$")
hold on
plot(r_noSource,phi_noSource,'-','Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName','no source')
% figure setting - velocity
figure(4)
hLegend = legend('Location','southeast','Interpreter','latex');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
set(gcf,'Position',[100 100 560*1 420])
xlim([0,6000])
xlabel("$r$")
ylabel("$v_"+string(n)+"$")
hold on
plot(r_noSource,vDM_noSource,'-','Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName','no source')
% figure setting - composing velocity
figure(5)
hLegend = legend('Location','northeastoutside','Interpreter','latex');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
set(gcf,'Position',[100 100 560*1.4 420])
xlabel("$r$")
ylabel("$v_{Tot}$")
xlim([0,6000])
hold on
plot(nan,nan,'-','Color',[0.7,0.7,0.7],'LineWidth',1.5,'DisplayName', "$v_{Tot}$")
plot(nan,nan,'--','Color',[0.7,0.7,0.7],'LineWidth',1.5,'DisplayName', "$v_{\rho}$")
plot(nan,nan,':','Color',[0.7,0.7,0.7],'LineWidth',1.5,'DisplayName', "$v_"+string(n)+"$")

% add plot as a varies
for i=1:length(aScale)
    % load computed data
    nameFile = outDir+"aVary/aScale"+string(round(aScale(i),2)*100)+"rScale"+string(round(rScale,2)*100);
    load(nameFile)
    % common setting
    label = "$s_a=$"+string(round(aScale(i),2));
    transp = i/(length(aScale));
    % add to figure - velocity
    figure(1)
    plot(rDM_Num,sqrt(rhoExt),'-','Color',[colorVec_default(2,:),transp],'LineWidth',1.5,'DisplayName',label)
    % add to figure - eigenfunction
    figure(2)
    plot(rDM_Num,f_Num,'-','Color',[colorVec_default(2,:),transp],'LineWidth',1.5,'DisplayName',label)
    % add to figure - potential
    figure(3)
    plot(rDM_Num,phi_Num,'-','Color',[colorVec_default(2,:),transp],'LineWidth',1.5,'DisplayName',label)
    % add to figure - velocity
    figure(4)
    plot(rDM_Num,vDM_Num,'-','Color',[colorVec_default(2,:),transp],'LineWidth',1.5,'DisplayName', label)
    % add to figure - composing velocity
    figure(5)
    plot(rDM_Num,vDM_Num,    ':', 'Color',[colorVec_default(2,:),transp],'LineWidth',1.5,'HandleVisibility','off');%'DisplayName', label + ", $V_{DM}$")
    plot(rDM_Num,vSource_Num,'--','Color',[colorVec_default(2,:),transp],'LineWidth',1.5,'HandleVisibility','off');%,'DisplayName', label + ", $V_{Source}$")
    plot(rDM_Num,vTot_Num,   '-', 'Color',[colorVec_default(2,:),transp],'LineWidth',1.5,'DisplayName', label)
end

% save figure - velocity
figure(1)
saveas(gcf,outDir+"Figures/varyA_density.png")
% save figure - eigenfunction
figure(2)
saveas(gcf,outDir+"Figures/varyA_eigenfunction.png")
% save figure - potential
figure(3)
saveas(gcf,outDir+"Figures/varyA_potential.png")
% save figure - velocity
figure(4)
saveas(gcf,outDir+"Figures/varyA_velocity.png")
% save figure - composing velocity
figure(5)
saveas(gcf,outDir+"Figures/varyA_composingVelocity.png")

