function plot3dRes_slices(r,f,k,ftype,viewtype,savef,outputPath)
    % Create 3D visualization of spherically symmetric functions
    % Inputs:
    %   r: vector of radial positions
    %   f, phi: vectors of function values at positions r
    %   k: excitation index
    %   ftype: string, 'abs'|'square' - how to process f values
    %   viewtype: string, '2d'|'3d' - type of visualization
    
    % Process function values based on ftype
    switch ftype
        case 'abs'
            fplot = abs(f);
        case 'square'
            fplot = f.^2;
    end
    
    % Find last extremum and set plot radius
    [rmax,fmax,~] = findLocMaxima(r,abs(fplot),k);
    last_extremum = 1.4 * rmax(k+1); % plotted radius
    plot_radius = 1.05 * last_extremum;   % plotted axis
    
    % Sort maxima to find second highest for colorbar
    sorted_maxima = sort(fmax, 'descend');
    color_max = sorted_maxima(2);  % Second highest maximum
    
    % Create 3D grid
    [X, Y, Z] = sphere(150);
    
    % Create interpolation function for radial data
    fr = griddedInterpolant(r, fplot, 'linear', 'nearest');
    
    % Calculate radial distance for each point
    R = sqrt(X.^2 + Y.^2 + Z.^2);
    
    % Scale sphere points and interpolate function values
    V = zeros(size(R));
    for i = 1:size(R,1)
        for j = 1:size(R,2)
            if R(i,j) <= last_extremum
                V(i,j) = fr(last_extremum);
            end
        end
    end
    
    figure(); 
    colormap jet
    
    if strcmp(viewtype, '2d')
        % 2D slice visualization
        [X2D, Y2D] = meshgrid(linspace(-plot_radius, plot_radius, 400));
        Z2D = zeros(size(X2D));
        R2D = sqrt(X2D.^2 + Y2D.^2);
        V2D = fr(R2D);
        V2D(R2D > last_extremum) = NaN;
        
        surf(X2D, Y2D, Z2D, V2D);
        view(2);
        shading interp;
        c = colorbar ;
        clim([0, color_max]);  % Set colorbar range
        
        % Add text in top right corner
        text(plot_radius*0.9, plot_radius*0.9, 0, ['n = ' num2str(k)], ...
             'HorizontalAlignment', 'right');
        
    else % 3D visualization with cutout
        % Create 7/8 of sphere by masking out one octant
        mask = ~(X >= 0 & Y >= 0 & Z >= 0);
        X(~mask) = NaN;
        Y(~mask) = NaN;
        Z(~mask) = NaN;
        V(~mask) = NaN;
        
        % Plot 7/8 of sphere
        s = surf(X.*plot_radius, Y.*plot_radius, Z.*plot_radius, V);
        hold on;
        
        % Create three slices for the cutout region
        [Xs, Ys, Zs] = meshgrid(linspace(0, plot_radius, 50));
        Vs = zeros(size(Xs));
        R3D = sqrt(Xs.^2 + Ys.^2 + Zs.^2);
        
        for i = 1:size(Xs,1)
            for j = 1:size(Xs,2)
                for m = 1:size(Xs,3)
                    if R3D(i,j,m) <= last_extremum
                        Vs(i,j,m) = fr(R3D(i,j,m));
                    else
                        Vs(i,j,m) = NaN;
                    end
                end
            end
        end
        
        % Plot slices in the cutout region
        slice(Xs, Ys, Zs, Vs, plot_radius, 0, 0);
        slice(Xs, Ys, Zs, Vs, 0, plot_radius, 0);
        slice(Xs, Ys, Zs, Vs, 0, 0, plot_radius);
        
        view(142.5, 30);
        shading interp;
        c = colorbar;
        clim([0, color_max]);  % Set colorbar range
        
        % Add text in visible corner
        text(plot_radius*0.6, -plot_radius*0.7, plot_radius*0.8,['n = ' num2str(k)], ...
             'HorizontalAlignment', 'right');
    end
    
    % Set common properties
    axis equal;
    ax = gca;
    % ax.Position = [0.1 0.1 0.8 0.8];  % Maximize plot size in figure
    xlabel('x');
    ylabel('y');
    zlabel('z');
    % title(['Function visualization: ' ftype]);
    
    % Set axis limits to fill figure
    axlim = plot_radius;  % Extend limits a bit beyond the sphere
    xlim([-axlim axlim]);
    ylim([-axlim axlim]);
    zlim([-axlim axlim]);

    % Save f figure
    if strcmp(savef, 'y')
        namefile = outputPath+"/section_"+ftype+viewtype+"_";
        if(k<10)
          saveas(gcf,namefile+"_k00"+k+".png");
        elseif k<100
          saveas(gcf,namefile+"_k0"+k+".png");
        else
          saveas(gcf,namefile+"_k"+k+".png");
        end
    end



end
