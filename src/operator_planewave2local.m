function Tpw2poly = operator_planewave2local(p, h0, nf, max_level)
    Minc = cell(1, max_level+1);
    [rvec, V] = approx.chebvander(p);
    Vi = inv(V);
    for l=0:max_level
        rl = 1/2.^l;
        hl = h0 / rl;
        % Convert proxy charges to outgoing
        k = hl*(-nf:nf)';
        Minc{l+1} = exp(1i*rvec.*k'/2*rl);
    end
    function Lambda = apply(Psi, level)
    % Convert incoming expansion Psi to local expansion Lambda
        % Evaluate incoming at proxy points
        proxy_vals = approx.kronmat_apply(Minc{level+1}, Psi, 3);
        % Convert to local expansion
        Lambda = approx.kronmat_apply(Vi, proxy_vals, 3);
    end
    Tpw2poly = @apply;
end
