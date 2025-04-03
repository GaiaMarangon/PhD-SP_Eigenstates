function deriv = fdFirstOrdDeriv(f,j,h,rule)
    %given as domain  a grid of points r1,...,rn with spacing h
    %computes the finite difference approximation of the first derivative
    %of the function f=f(r1),...,f(rn) at a point rj 
    %according to the selected rule: forward(f), backward(b), centered)
    switch(rule)
        case 'f'
            if (j> length(f)-2) 
                fprintf("forward rule cannot be applied here\n")
                deriv = NaN;
            else
                deriv = ( -3*f(j) +4*f(j+1) -f(j+2) ) /(2*h);
            end
        case 'b'
            if (j< 3) 
                fprintf("backward rule cannot be applied here\n")
                deriv = NaN;
            else
                deriv = ( 3*f(j) -4*f(j-1) +f(j-2) ) /(2*h);
            end
        case 'c'
            if or( j>length(f)-1, j<2) 
                fprintf("centered rule cannot be applied here\n")
                deriv = NaN;
            else
                deriv = ( f(j+1) -f(j-1) ) /(2*h);
            end
        otherwise 
            fprintf("invalid rule for first order derivative\n");
    end
end