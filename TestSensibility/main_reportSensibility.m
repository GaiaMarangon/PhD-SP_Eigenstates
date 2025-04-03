clearvars, clc,close all

%--- PARAMETERS ----------------------------------------------------------
%nodes
k=80;
%loading and saving parameters
pathfile = "results_NormOne/";  %output folder
nametest = "results";           %file name basis

%--- SENSIBILITY ANALYSIS ------------------------------------------------
%settingifle name
if k<10
    namek = "_k00";
elseif k<100
    namek = "_k0";
else
    namek = "_k";
end
%loading data
load(pathfile+nametest+namek+k+".mat",'tolInt','tolExt','h','nSens','eigval','r','f','u1','elapsedTime')
%opening file for report 
fileID = fopen(pathfile+nametest+"_reportSens"+namek+k+".txt",'w');
    %--- eigval ------------------
    fprintf(fileID,'eigenvalues:\n');
    figure(1)
    hold on
    %nSens = 1;
    for i=1:nSens
        fprintf(fileID,'\th: %f\teigval: %.10e\n',h(i),eigval{i});
        plot(h(i),eigval{i},'o')
    end
    xlabel('h')
    ylabel('eigval(h)')
    hold off
    %--- error eigenval --------------------
    fprintf(fileID,'err eigenvalues:\n');
    figure(2)
    hold on
    for i=2:nSens
        fprintf(fileID,'\th: %f\tabsolute err: %.10e\n',h(i),abs(eigval{i}-eigval{i-1}));
        plot(h(i),abs(eigval{i}-eigval{i-1}),'o')
    end
    xlabel('h')
    ylabel('|eigval(h)-eigval(h_{prev})|')
    hold off
    %
    fprintf(fileID,'\n');
    for i=2:nSens
        fprintf(fileID,'\th: %f\tpercent  err: %.10e\n',h(i),abs(eigval{i}-eigval{i-1})/eigval{i});
    end
    %--- maxima position -----------
    rmax = cell(nSens,1);
    fprintf(fileID,'\nloc maxima difference (2-norm):\n');
    figure(3)
    hold on
    for i=1:nSens
        [rmax{i},~,~]=findLocMaxima(r{i},abs(f{i}),k);
        if i>1
            res = sqrt(sum( (rmax{i}-rmax{i-1}).^2 ));
            avgerr = sum( abs(rmax{i}-rmax{i-1}) ) / length(rmax{i});
            fprintf(fileID,'\ths: %f,\t%f\tnorm percent: %.5e\taverage percent local err: %.5e\n', ...
                h(i-1),h(i),res/r{i}(end), avgerr/r{i}(end));
            plot(h(i),avgerr/r{i}(end),'o')
        end
    end
    xlabel('h')
    ylabel('|r_{max}(h)-r_{max}(h_{prev})|_1 / n_{max}')
    hold off
    %--- r95 -----------------
    r95 = zeros(nSens,1);
    fprintf(fileID,'\nradius:\n');
    for i=1:nSens
        partInt = 0;
        for j=2:length(r{i})
            if partInt < 0.95
                rule1 = 'c';
                rule2 = 'c';
                if j==2
                    rule1='f';
                elseif j== length(r{i})
                    rule2='b';
                end
                partInt = partInt+CTintegrate(j-1,j,u1{i}.^2,r{i}(j)-r{i}(j-1),rule1,rule2);
            else
                r95(i) = r{i}(j);
                break;
            end
        end
        fprintf(fileID,'\tr95: %.5e\trEnd: %.5e\tpercent dist: %.5e\n', ...
            r95(i),r{i}(end),abs(r95(i)-r{i}(end))/r{i}(end));
    end
    %--- elapsed time --------
    fprintf(fileID,'\nelapsed times:\n');
    figure(5)
    hold on
    for i=1:nSens
        fprintf(fileID,'\th: %d\ttime (s): %.3f\n',h(i),elapsedTime{i});
        plot(h(i),elapsedTime{i},'o')
    end
    xlabel('h')
    ylabel('elapsed time (s)')
    hold off
fclose(fileID);
%-------------------------------------------------------------------------
