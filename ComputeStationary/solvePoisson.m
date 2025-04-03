function phi = solvePoisson(eps,r,u)
    h=r(2)-r(1);
    ngrid = length(r);
    uquad = u.*u;
    uquads =  uquad./(r+eps*(r==0));
    phi = zeros(ngrid,1);
    for i=1:ngrid
        %regularize 1/r for r=0
        r(i)= r(i) + eps*(r(i)==0);
        %choose approx rule for first derivative at r
        if i==1
            rule2 = "f";
        elseif i==ngrid
            rule2 = "b";
        else
            rule2 = "c";
        end
        %compute phi
        phi(i) = -( CTintegrate(1,i,uquad,h,'f',rule2) )/r(i)  ...
                 +( CTintegrate(1,i,uquads,h,'f',rule2) )  ...
                 -( CTintegrate(1,ngrid,uquads,h,'f','b') ); 
    end
end