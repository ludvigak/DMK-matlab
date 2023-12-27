function Tfar = operator_far(p, hf, nf, Ctrunc, sigma)
% Operator that converts proxy_charges to local exapnsion of far field potential
    k = hf*(-nf:nf);
    [k1, k2, k3] = ndgrid(k);
    k1 = k1(:); k2 = k2(:); k3 = k3(:);
    [rvec, V] = approx.chebvander(p);
    kdotr = exp(-1i*k(:).*rvec'/2);
    M = V\(kdotr');
    W0hat = laplace_winkernel_fourier(k1, k2, k3, sigma, Ctrunc);
    w = W0hat  * hf^3 / (2*pi)^3;
    function ufar_expa = apply(proxy_charges)
        ghat = approx.kronmat_apply(kdotr, proxy_charges, 3);
        uhat = w.*ghat;
        ufar_expa = real( approx.kronmat_apply(M, uhat, 3) );
    end
    Tfar = @apply;
end
