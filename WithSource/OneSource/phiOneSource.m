function phi = phiOneSource(r,params,rhoName)
    % Computes the contribution to phi(r) due to the external density rho(r), specified by rhoName and by parameters params, on domain r.
    % Available density types: "Exp", "TruncatedPlummer", "HardBall".
    % For some densities, the analitical expression for phi(r) is used to speed up computations ("Exp", "HardBall").
    % In other cases, the general expression for a sperically symmetric potential is evaluated numerically.
    % For numerical evaluation, regularization parameter eps is expected as the last parameter in params. 

    switch rhoName
        case "Exp"
            % Get parameters
            r0 = params(1);
            a  = params(2);
            % Compute potential
            phi = a*r0^2./r .*( -2*r0 + exp(-r/r0).*(r+2*r0) );

        case "HardBall"
            % Get parameters
            r0 = params(1);
            a  = params(2);
            % Compute potential
            [~,idxr]  = min(abs(r-r0));
            phi = zeros(size(r));
            phi(1:idxr) = -a*(r(1:idxr).^2/3 + (r0^2-r(1:idxr).^2)/2);
            phi(idxr+1:end) = -a*r0^3/3 ./r(idxr+1:end);

        case "TruncatedPlummer"
            % Get density
            rho = rhoExtSource(r,params,rhoName);
            % Compute potential
            phi = phiGenericRho(r,rho,params(end));

        otherwise
            phi = zeros(size(r));
    end
end
    
    