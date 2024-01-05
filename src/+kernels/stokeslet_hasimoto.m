classdef stokeslet_hasimoto < kernels.SplitKernelInterface
% Stokes kernel split using the Hasimoto decomposition

    properties (Constant)
        name    = "Stokes Hasimoto";
        dim_in  = 3;
        dim_out = 3;
    end
    
    methods (Static)
        function [Sdiag, Soffd] = real_decay(r, xi)
            Sdiag = -2*xi/sqrt(pi)*exp(-xi^2*r.^2).*r + erfc(xi*r);
            Soffd =  2*xi/sqrt(pi)*exp(-xi^2*r.^2).*r + erfc(xi*r);
        end

        function B = fourier_scaling(ksq, xi)
            B = 8*pi./ksq.^2 .* (1 + ksq/(4*xi^2)) .*  exp(-ksq/4/xi^2);
            B(ksq==0) = 0;
        end
        
        function u = direct(targets, sources, f)
        % Stokeslet kernel S(r)
            assert(size(sources, 2)==3)
            if size(targets, 1) ~= 3
                targets = targets.';
            end
            assert(size(targets, 1)==3);
            Ntrg = size(targets, 2);
            u = zeros(Ntrg, 3);
            batch_size = 64;
            % Vectorize over batch size
            for idx_begin=1:batch_size:Ntrg
                idx_end = min(Ntrg, idx_begin+batch_size-1);
                range =  idx_begin:idx_end;
                r1 = targets(1, range)-sources(:, 1);
                r2 = targets(2, range)-sources(:, 2);
                r3 = targets(3, range)-sources(:, 3);
                
                r = sqrt(r1.^2 + r2.^2 + r3.^2);
                udiag = f' * (1./r);
                fdotrd3 = (f(:, 1).*r1 + f(:, 2).*r2  + f(:, 3).*r3)./r.^3;
                uoffd = [sum(r1.*fdotrd3, 1); sum(r2.*fdotrd3, 1); sum(r3.*fdotrd3, 1)];
                u(range, :) = (udiag + uoffd).';
            end
        end
    end

    methods
        function ures = reskernel(self, targets, sources, f, sigma_l)
        % Stokes residual kernel R_l(r)
            assert(size(sources, 2)==3)
            if size(targets, 1) ~= 3
                targets = targets.';
            end
            assert(size(targets, 1)==3);
            r1 = targets(1, :)-sources(:, 1);
            r2 = targets(2, :)-sources(:, 2);
            r3 = targets(3, :)-sources(:, 3);
            r = sqrt(r1.^2 + r2.^2 + r3.^2);
            xi = 1/sigma_l;
            [Sdiag, Soffd] = self.real_decay(r, xi);
            udiag = f' * (Sdiag./r);
            fdotrd3 = (f(:, 1).*r1 + f(:, 2).*r2  + f(:, 3).*r3).*Soffd./r.^3;
            uoffd = [sum(r1.*fdotrd3, 1); sum(r2.*fdotrd3, 1); sum(r3.*fdotrd3, 1)];
            ures = (udiag + uoffd).';
        end

        function umoll = mollkernel(self, targets, sources, f, sigma_l)
        % Mollified kernel M_l = S(r) - R_l(r)
            u = self.direct(targets, sources, f);
            ures = self.reskernel(targets, sources, f, sigma_l);
            umoll = u - ures;
        end

        function udiff = diffkernel(self, targets, sources, f, sigma_l)
        % Laplace difference kernel D_l(r) = M_{l+1}(r) - M_l(r)
            sigma_lp1 = sigma_l/2; % sigma_{l+1} = sigma_l / 2
            udiff = mollkernel(self, targets, sources, f, sigma_lp1) - ...
                    mollkernel(self, targets, sources, f, sigma_l);
        end
        
        function Dlhat = diffkernel_fourier(self, k1, k2, k3, sigma_l)
        % Fourier transform of Laplace difference kernel
            % sigma_lp1 = sigma_l/2; % sigma_{l+1} = sigma_l / 2
            % ksq = k1.^2 + k2.^2 + k3.^2;
            % Dlhat = 4*pi*(exp(-ksq*sigma_lp1^2/4) - exp(-ksq*sigma_l^2/4))./ksq;
            % Dlhat(ksq==0) = pi*(sigma_l^2-sigma_lp1^2);
        end

        function W0hat = winkernel_fourier(self, k1, k2, k3, sigma_0, Ctrunc)
        % Fourier transform of  windowed mollified Laplace kernel (i.e. far field)
            % ksq = k1.^2 + k2.^2 + k3.^2;
            % k = sqrt(ksq);
            % W0hat = 8*pi*(sin(Ctrunc*k/2)./k).^2 .* exp(-ksq * sigma_0^2/4);
            % W0hat(ksq==0) = 8*pi * Ctrunc^2/4;
        end

        function uself = self_interaction(self, charges, sigma_l)
        % uself = -charges*2/(sqrt(pi)*sigma_l);
        end
    end
end
