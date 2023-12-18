function [rvec, V]  = chebvander(p)
% p: integer
%
% Returns rx1 first kind Chebyshev nodes r_i
% and pxp pseudo-Vandermonde matrix 
% V_{ij} = T_{j-1}(r_i)
% such that Vc=f(rvec) where c are the Chebyshev coefficients of f
    rvec = chebpts(p, 1);
    V = approx.chebevalmat(rvec, p);
end
