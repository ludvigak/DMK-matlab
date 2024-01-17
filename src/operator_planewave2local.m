function Tpw2poly = operator_planewave2local(p, h0, nf, max_level, kernel)
    ViMinc = cell(1, max_level+1);
    Dlhat  = cell(1, max_level+1);
    w      = cell(1, max_level+1);
    [rvec, V] = approx.chebvander(p);
    for l=0:max_level
        rl = 1/2.^l;
        hl = h0 / rl;
        k = hl*(-nf:nf)';
        Minc = exp(1i*rvec.*k'/2*rl);
        ViMinc{l+1} = V\Minc;
        [k1, k2, k3] = ndgrid(k, k, k);
        % Dlhat is a closure that is applied to Psi at runtime,
        % to account for vectors
        Dlhat{l+1} = kernel.diffkernel_fourier(k1(:), k2(:), k3(:), l);
        % w_l factor 
        w{l+1} = 1/(2*pi)^3 * hl^3;
    end
    dim = kernel.dim_in;
    function Lambda = apply(Psi, level)
    % Convert incoming expansion Psi to local expansion Lambda
        wl = w{level+1}; % Fourier kernel scaling
        op = Dlhat{level+1};
        uhat = wl.*op(Psi);
        % Evaluate incoming at proxy points and convert to local expansion
        Lambda = zeros(p^3, dim);
        for d=1:dim
            Lambda(:, d) = real( approx.kronmat_apply(ViMinc{level+1}, uhat(:, d), 3) );
        end
    end
    Tpw2poly = @apply;
end

