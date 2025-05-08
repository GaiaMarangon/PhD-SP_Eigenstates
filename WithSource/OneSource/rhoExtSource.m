function rho = rhoExtSource(r,params,rhoName)
    % Compute density specified by rhoName, with parameters params, on domain r.
    % Available density types: "Exp", "TruncatedPlummer", "HardBall":
    %  - rhoExp =  a e^(-r/r0)
    %  - rhoTruncatedPlummer =  a / ( (r/r0)^2 +1 )^2.5     for r<=r0
    %                           rho(r0) * exp(r0-r)         for r> r0
    %  - rhoHardBall =  a      for r<=r0
    %                   0      for r> r0
    

    switch rhoName
        case "Exp"
            % define parameters
            r0 = params(1);
            a  = params(2);
            % define density
            rho = a*exp(-r/r0);
        case "TruncatedPlummer"
            % define parameters
            r0 = params(1);
            a  = params(2);
            % define density
            rho = zeros(size(r));
            [~,idxr0]  = min(abs(r-r0));
            rho(1:idxr0) = a ./ ( (r(1:idxr0)/r0).^2 +1 ).^2.5;
            rhor0 = rho(idxr0);
            rho(idxr0:end) =  rhor0 .* exp(r0-r(idxr0:end));
        case "HardBall"
            % define parameters
            r0 = params(1);
            a  = params(2);
            % define density
            rho = zeros(size(r));
            [~,idxr]  = min(abs(r-r0));
            rho(1:idxr) = ones(idxr,1)*a;     
        otherwise
            rho = zeros(size(r));
    end
end