function Tfar = operator_periodic(p, nf, kernel)
    % Operator that converts proxy_charges to local expansion of periodic potential
    k = 2*pi*(-nf:nf);
    [k1, k2, k3] = ndgrid(k);
    k1 = k1(:); k2 = k2(:); k3 = k3(:);
    [rvec, V] = approx.chebvander(p);
    kdotr = exp(-1i*k(:).*rvec'/2);
    M = V\(kdotr');
    M0hat = kernel.mollkernel_fourier(k1, k2, k3, 0);
    Nf = numel(k1);
    dim_out = kernel.dim_out;
    function [ufar_expa, uhat] = apply(proxy_charges)
        dim_in = size(proxy_charges, 2);
        ghat = zeros(Nf, dim_in, like=1+1i);
        for d=1:dim_in
            ghat(:, d) = approx.kronmat_apply(kdotr, proxy_charges(:, d), 3);
        end
        uhat = M0hat(ghat);
        ufar_expa = zeros(p^3, dim_out);
        for d=1:dim_out
            ufar_expa(:,d) = real( approx.kronmat_apply(M, uhat(:, d), 3) );
        end
    end
    Tfar = @apply;
end
