function Tpw2poly = operator_planewave2local(p, h0, nf, max_level)
    ViMinc = cell(1, max_level+1);
    [rvec, V] = approx.chebvander(p);
    for l=0:max_level
        rl = 1/2.^l;
        hl = h0 / rl;
        k = hl*(-nf:nf)';
        Minc = exp(1i*rvec.*k'/2*rl);
        ViMinc{l+1} = V\Minc;
    end
    function Lambda = apply(Psi, level)
    % Convert incoming expansion Psi to local expansion Lambda
    % Evaluate incoming at proxy points and convert to local expansion
        dim = size(Psi, 2);
        Lambda = zeros(p^3, dim);
        for d=1:dim
            Lambda(:, d) = real( approx.kronmat_apply(ViMinc{level+1}, Psi(:, d), 3) );
        end
    end
    Tpw2poly = @apply;
end

