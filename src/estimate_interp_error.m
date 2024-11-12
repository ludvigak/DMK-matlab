function est = estimate_interp_error(p_list, kernel)
% Estimate DMK interpolation error (relative) by measuring
% error when interpolating difference kernel along a 1D line at level 0
    Nsrc  = 100;
    Neval =  25;
    rs = RandStream('mt19937ar', seed=1);
    charges = rs.rand(Nsrc, kernel.dim_in) - 1/2;
    sources = rs.rand(Nsrc, 3) - 1/2;
    xeval = linspace(-1/2, 1/2, Neval)';
    level = 0;
    K = @(x) kernel.diffkernel([x(:) 0*x(:), 0*x(:)], sources, charges, level);
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

