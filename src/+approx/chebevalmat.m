function E = chebevalmat(xvec, p)
% xvec: Nx1 vector of points in [-1, 1]
% p: integer
%
% Return Nxp matrix 
% E_{ij} = T_{j-1}(x_i)
% such that E*c evaluates Chebyshev expansion with coefficients c
% at points xvec
    [N, m] = size(xvec);
    assert(m==1);
    E = zeros(N, p);
    E(:, 1) = 1;
    E(:, 2) = xvec;
    for n=3:p
        E(:, n) = 2*xvec.*E(:, n-1) - E(:, n-2);
    end
end
