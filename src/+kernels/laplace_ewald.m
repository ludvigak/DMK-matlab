classdef laplace_ewald < kernels.SplitKernelInterface
% Laplace kernel split using the orignal Ewald decomposition

    properties (Constant)
        name    = "Laplace Ewald";
        dim_in  = 1;
        dim_out = 1;
    end
    
    methods (Static)
        function u = direct(targets, points, charges)
            N = numel(charges);
            assert(size(points, 1)==N)
            assert(size(points, 2)==3)
            if size(targets, 1) ~= 3
                targets = targets.';
            end
            assert(size(targets, 1)==3);
            u = zeros(size(targets, 2), 1);
            batch_size = 64;
            % Vectorize over batch size
            for idx_begin=1:batch_size:N
                idx_end = min(N, idx_begin+batch_size-1);
                range =  idx_begin:idx_end;
                r = sqrt( (points(range, 1)-targets(1, :)).^2 + ...
                          (points(range, 2)-targets(2, :)).^2 + ...
                          (points(range, 3)-targets(3, :)).^2 );
                K = 1./r;
                K(r==0) = 0;
                u = u + (charges(range).' * K).';
            end
        end

        function ures = reskernel(targets, points, charges, sigma_l)
        % Laplace residual kernel R_l(r)
            if size(targets, 1) ~= 3
                targets = targets.';
            end
            assert(size(targets, 1)==3);
            r = sqrt( (points(:, 1)-targets(1, :)).^2 + ...
                      (points(:, 2)-targets(2, :)).^2 + ...
                      (points(:, 3)-targets(3, :)).^2 );
            Rl = erfc(r/sigma_l)./r;
            Rl(r==0) = 0;
            ures = (charges.' * Rl).';
        end

        function udiff = diffkernel(targets, points, charges, sigma_l)
        % Laplace difference kernel D_l(r)
            if size(targets, 1) ~= 3
                targets = targets.';
            end
            assert(size(targets, 1)==3);
            sigma_lp1 = sigma_l/2; % sigma_{l+1} = sigma_l / 2
            r = sqrt( (points(:, 1)-targets(1, :)).^2 + ...
                      (points(:, 2)-targets(2, :)).^2 + ...
                      (points(:, 3)-targets(3, :)).^2 );
            Dl = (erf(r/sigma_lp1) - erf(r/sigma_l))./r;
            udiff = Dl.' * charges;
        end
        
        function Dlhat_fun = diffkernel_fourier(k1, k2, k3, sigma_l)
        % Fourier transform of Laplace difference kernel
        % Returns operator that applies kernel
            sigma_lp1 = sigma_l/2; % sigma_{l+1} = sigma_l / 2
            ksq = k1.^2 + k2.^2 + k3.^2;
            Dlhat = 4*pi*(exp(-ksq*sigma_lp1^2/4) - exp(-ksq*sigma_l^2/4))./ksq;
            Dlhat(ksq==0) = pi*(sigma_l^2-sigma_lp1^2);
            function uhat=apply(Psi)
                uhat = Dlhat.*Psi;
            end
            Dlhat_fun = @apply;
        end

        function W0hat_fun = winkernel_fourier(k1, k2, k3, sigma_0, Ctrunc)
        % Fourier transform of  windowed mollified Laplace kernel (i.e. far field)
            ksq = k1.^2 + k2.^2 + k3.^2;
            k = sqrt(ksq);
            W0hat = 8*pi*(sin(Ctrunc*k/2)./k).^2 .* exp(-ksq * sigma_0^2/4);
            W0hat(ksq==0) = 8*pi * Ctrunc^2/4;
            function uhat=apply(fhat)
                uhat = W0hat .* fhat;
            end
            W0hat_fun = @apply;
        end

        function uself = self_interaction(charges, sigma_l)
            uself = -charges*2/(sqrt(pi)*sigma_l);
        end
    end
end
