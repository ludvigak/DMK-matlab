function Tpw2poly = operator_planewave2local(p, h0, nf, max_level, sigma_0, kernel)
    ViMinc = cell(1, max_level+1);
    w      = cell(1, max_level+1);
    [rvec, V] = approx.chebvander(p);
    for l=0:max_level
        rl = 1/2.^l;
        hl = h0 / rl;
        sigma_l = sigma_0 / 2^l;
        k = hl*(-nf:nf)';
        Minc = exp(1i*rvec.*k'/2*rl);
        ViMinc{l+1} = V\Minc;
        [k1, k2, k3] = ndgrid(k, k, k);
        % Todo: Dlhat needs to be a closure that is applied to Psi at runtime,
        % to account for vectors
        Dlhat = kernel.diffkernel_fourier(k1(:), k2(:), k3(:), sigma_l);
        % w_l factor 
        w{l+1} = 1/(2*pi)^3 * hl^3 * Dlhat;
    end
    function Lambda = apply(Psi, level)
    % Convert incoming expansion Psi to local expansion Lambda
        wl = w{level+1}; % Fourier kernel scaling
        uhat = wl.*Psi;
        % Evaluate incoming at proxy points and convert to local expansion
        Lambda = real( approx.kronmat_apply(ViMinc{level+1}, uhat, 3) );
    end
    Tpw2poly = @apply;
end

