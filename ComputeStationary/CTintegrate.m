function integr = CTintegrate(i,j,f,h,rule1,rule2)
%given as domain  a grid of points r1,...,rn
%here we integrate from a=ri to b=rj the function f(r),
%which is given as f(r1),...,f(rn) on the grid points
%with corrected trapezoidal integration rule, with
%first order deriv at the extremals approximated according 
%to the specified fd rule (rule1 for a=ri, rule2 for b=rj) 
    Th = sum( ( f(i:j-1) + f(i+1:j) ) *h/2 );
    deriv1 = fdFirstOrdDeriv(f,i,h,rule1);
    deriv2 = fdFirstOrdDeriv(f,j,h,rule2);
    integr = Th -( deriv2-deriv1 )*h^2/12;
end