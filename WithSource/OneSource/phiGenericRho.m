function phi = phiGenericRho(r,rho,eps)
    % Compute the potential phi(r) fro a generic source rho(r), on domain r.    
    % Parameter eps is used to regularize at r=0.

    phi =  zeros(size(r));
    for i=1:length(r)
        % regularize 1/r for r=0
        r(i)= r(i) + eps*(r(i)==0);
        % choose approx rule for first derivative at r
        if i==1
            rule2 = "f";
        elseif i==length(r)
            rule2 = "b";
        else
            rule2 = "c";
        end
        % compute phi
        phi(i) = -( CTintegrate(1,i,        rho.*r.^2,r(2)-r(1),'f',rule2) )/r(i) ...
                 -  CTintegrate(1,length(r),rho.*r,   r(2)-r(1),'f','b') ...
                 +  CTintegrate(1,i,        rho.*r,   r(2)-r(1),'f',rule2);
    end
end