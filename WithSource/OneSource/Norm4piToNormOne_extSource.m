function [eigval, f, phi, r, rho]=Norm4piToNormOne_extSource(eigval4pi, f4pi, phi4pi, r4pi, rho4pi)
    % from 4pi normalization to one normalization
    eigval = eigval4pi /(4*pi)^2;
    f      = f4pi      /(4*pi)^2;
    phi    = phi4pi    /(4*pi)^2;
    r      = r4pi      *(4*pi);
    rho    = rho4pi    /(4*pi)^4;
end