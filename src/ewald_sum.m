function [u, u_moll, u_res, u_self] = ewald_sum(points, charges, kernel)
% Direct Ewald sum, no optimizations
% Assumes unit box

    N = size(points, 1);

    %% Residual sum
    % Collect 3^3 nearest neighbor boxes
    u_res = 0;
    for p1=-1:1
        for p2 = -1:1
            for p3 = -1:1
                shift = [p1 p2 p3];
                src = points+shift;
                trg = points;
                u_res = u_res + kernel.reskernel(trg, src, charges, 0);
            end
        end
    end

    %% Fourier sum
    nf = ceil(kernel.Kmax/(2*pi));
    kvec = 2*pi*(-nf:nf);
    [k1, k2, k3] = ndgrid(kvec);
    k1 = k1(:); k2 = k2(:); k3 = k3(:);
    kr = k1*(points(:,1)).' + k2*(points(:,2)).' + k3*(points(:,3)).';
    if isa(kernel, 'kernels.StressletBase')
        f_charges = kernel.input_product(charges);
    else
        f_charges = charges;
    end
    f_dim = size(f_charges, 2);
    fhat = zeros(numel(k1), f_dim);
    for d=1:f_dim
        fhat(:, d) = sum(exp(-1i*kr) .* (f_charges(:, d)).', 2);
    end
    moll_op = kernel.mollkernel_fourier(k1, k2, k3, 0);
    GF = moll_op(fhat);
    u_moll = zeros(N, kernel.dim_out);
    for d=1:kernel.dim_out
        u_moll(:, d) = (sum( real(GF(:, d).*exp(1i*kr)), 1))';
    end

    %% Self interaction
    u_self = kernel.self_interaction(charges, 0);

    %% Return sum
    u = u_moll + u_res + u_self;
end
