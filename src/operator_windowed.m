function Tfar = operator_windowed(p, hf, nf, Ctrunc, sigma, kernel)
% Operator that converts proxy_charges to local expansion of windowed potential
    k = hf*(-nf:nf);
    [k1, k2, k3] = ndgrid(k);
    k1 = k1(:); k2 = k2(:); k3 = k3(:);
    [rvec, V] = approx.chebvander(p);
    kdotr = exp(-1i*k(:).*rvec'/2);
    M = V\(kdotr');
    W0hat = kernel.winkernel_fourier(k1, k2, k3, sigma, Ctrunc);
    w = hf^3 / (2*pi)^3;
    Nf = numel(k1);
    dim_in  = kernel.dim_in;
    dim_out = kernel.dim_out;
    function ufar_expa = apply(proxy_charges)
        ghat = zeros(Nf, dim_in, like=1+1i);
        for d=1:dim_in
            ghat(:, d) = approx.kronmat_apply(kdotr, proxy_charges(:, d), 3);
        end
        uhat = w.*W0hat(ghat);
        ufar_expa = zeros(p^3, dim_out);
        for d=1:dim_out
            ufar_expa(:,d) = real( approx.kronmat_apply(M, uhat(:, d), 3) );
        end
    end
    Tfar = @apply;
end
