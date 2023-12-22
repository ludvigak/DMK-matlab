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
    function Phi = apply(proxy_charges, level)
        Phi = approx.kronmat_apply(Mout{level+1}, proxy_charges, 3);
    end
    Tprox2pw = @apply;
end

