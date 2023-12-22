function Tpwshift = operator_planewave_shift(h0, nf)
    r0 = 1;
    k = h0*(-nf:nf)';
    S = cell(1,3);
    for i=-1:1
        S{i+2} = exp(-1i*k*r0*i);
    end
    function shifted = apply(expa, i, j, k)
    % (i,j,k) represent location of outgoing box vs incoming box in 3x3x3 grid
        shifted = expa .* kron(S{k+2}, kron(S{j+2}, S{i+2}));
    end
    Tpwshift = @apply;
end

