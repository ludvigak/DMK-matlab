function Tprox2pw = operator_proxy2planewave(p, h0, nf, max_level)
    Mout = cell(1, max_level+1);
    [rvec, ~] = approx.chebvander(p);
    for l=0:max_level
        rl = 1/2.^l;
        hl = h0 / rl;
        % Convert proxy charges to outgoing
        k = hl*(-nf:nf)';
        Mout{l+1} = exp(-1i*k.*rvec'/2*rl);
    end
    Nf = (2*nf+1)^3;
    function Phi = apply(proxy_charges, level)
        dim = size(proxy_charges, 2);
        Phi = zeros(Nf, dim, like=(1+1i));
        for d=1:dim
            Phi(:, d) = approx.kronmat_apply(Mout{level+1}, proxy_charges(:, d), 3);
        end
    end
    Tprox2pw = @apply;
end

