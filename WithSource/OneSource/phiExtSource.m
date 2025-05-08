function phi = phiExtSource(r,params,rhoNames,eps)
    % Computes the contribution to phi due to external sources, on domain r.
    % - params is a matrix, with each row containing the parameters for a single source.
    % - rhoNames is a vector, with the names of the density of each source.
    % - eps is the regularizing parameter used at r=0 in numerical evaluations.
    
    phi = zeros(size(r));
    % run through sources
    for i=1:length(rhoNames)
        % Get source data
        param = [params(i,:), eps];
        rhoName = rhoNames(i);
        % Compute phi contribution
        phi = phi + phiOneSource(r,param,rhoName);
    end
end