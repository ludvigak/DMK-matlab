function est = estimate_interp_error(p_list, sigma_0)
% Estimate DMK interpolation error (relative) by measuring
% error when interpolating difference kernel along a 1D line
    Nsrc  = 100;
    Neval =  25;
    rs = RandStream('mt19937ar', seed=1);
    xsrc = rs.rand(Nsrc, 1) - 1/2;
    charges = rs.rand(Nsrc, 1) - 1/2;
    sources = [xsrc zeros(Nsrc, 2)];   
    xeval = linspace(-1/2, 1/2, Neval)';
    K = @(x) laplace_diffkernel([x(:) 0*x(:), 0*x(:)], sources, charges, sigma_0);
    Kref = K(xeval);
    est = 0*p_list;
    for i=1:numel(p_list)
        p = p_list(i);
        [rvec, V] = approx.chebvander(p);
        xp   = rvec/2;
        Kval = K(xp);
        E = approx.chebevalmat(2*xeval, p);
        Ki = E*(V\Kval);
        est(i) = norm(Kref-Ki, inf) / norm(Kref, inf);
    end
end

