function [eigval, f, phi, r, u1]=Norm4piToNormOne(eigval4pi, f4pi, phi4pi, r4pi, u14pi)
    % from 4pi normalization to one normalization
    eigval = eigval4pi /(4*pi)^2;
    f      = f4pi      /(4*pi)^2;
    phi    = phi4pi    /(4*pi)^2;
    r      = r4pi      *(4*pi);
    u1     = u14pi     /(4*pi);
end