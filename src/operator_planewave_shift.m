function Tpwshift = operator_planewave_shift(h0, nf)
    r0 = 1;
    k = h0*(-nf:nf)';
    S = cell(1,3);
    for i=-1:1
        S{i+2} = exp(-1i*k*r0*i);
    end
    K = cell(3,3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                K{i,j,k} = kron(S{k}, kron(S{j}, S{i}));
            end
        end
    end
    
    function shifted = apply(expa, i, j, k)
    % (i,j,k) represent location of outgoing box vs incoming box in 3x3x3 grid
        shifted = expa .* K{i+2,j+2,k+2}; % TODO: speed up somehow
    end
    Tpwshift = @apply;
end

