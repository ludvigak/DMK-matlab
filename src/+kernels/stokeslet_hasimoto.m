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
            B(ksq==0) = 0; % Not actually zero, just avoiding inf
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
                self_mask = (r==0);
                rinv = 1./r;
                rinv(self_mask) = 0;
                udiag = f' * rinv;
                fdotrd3 = (f(:, 1).*r1 + f(:, 2).*r2  + f(:, 3).*r3)./r.^3;
                fdotrd3(self_mask) = 0;
                uoffd = [sum(r1.*fdotrd3, 1); sum(r2.*fdotrd3, 1); sum(r3.*fdotrd3, 1)];
                u(range, :) = (udiag + uoffd).';
            end
        end

        function uself = self_interaction(charges, sigma_l)
            uself = -4/(sigma_l*sqrt(pi)) * charges;
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
            Sdrinv = Sdiag./r;
            self_mask = (r==0);
            Sdrinv(self_mask) = 0;
            udiag = f' * Sdrinv;
            fdotrd3 = (f(:, 1).*r1 + f(:, 2).*r2  + f(:, 3).*r3).*Soffd./r.^3;
            fdotrd3(self_mask) = 0;
            uoffd = [sum(r1.*fdotrd3, 1); sum(r2.*fdotrd3, 1); sum(r3.*fdotrd3, 1)];
            ures = (udiag + uoffd).';
        end

        function umoll = mollkernel(self, targets, sources, f, sigma_l)
        % Mollified kernel M_l = S(r) - R_l(r)
            u = self.direct(targets, sources, f);
            ures = self.reskernel(targets, sources, f, sigma_l);
            umoll = u - ures;
        end

        function Mlhat_fun = mollkernel_fourier(self, k1, k2, k3, sigma_l)
            ksq = k1.^2 + k2.^2 + k3.^2;
            xi = 1/sigma_l;
            B = self.fourier_scaling(ksq, xi);
            Nf = numel(k1);
            function uhat=apply(fhat)
                assert(all(size(fhat)==[Nf 3]));
                kdotf = k1.*fhat(:, 1) + k2.*fhat(:, 2) + k3.*fhat(:, 3);
                uhat = B .* (ksq.*fhat - [k1 k2 k3].*kdotf);
                % Note k==0 artificially set to zero
            end
            Mlhat_fun = @apply;
        end
        
        function udiff = diffkernel(self, targets, sources, f, sigma_l)
        % Difference kernel D_l(r) = M_{l+1}(r) - M_l(r)
            sigma_lp1 = sigma_l/2; % sigma_{l+1} = sigma_l / 2
            udiff = mollkernel(self, targets, sources, f, sigma_lp1) - ...
                    mollkernel(self, targets, sources, f, sigma_l);
        end
        
        function Dlhat_fun = diffkernel_fourier(self, k1, k2, k3, sigma_l)
        % Fourier transform of difference kernel
            ksq = k1.^2 + k2.^2 + k3.^2;
            sigma_lp1 = sigma_l/2; % sigma_{l+1} = sigma_l / 2;
            xi_l = 1/sigma_l;
            xi_lp1 = 1/sigma_lp1;
            Bl = self.fourier_scaling(ksq, xi_l);
            Blp1 = self.fourier_scaling(ksq, xi_lp1);
            Bdiff = Blp1 - Bl;
            % Difference kernel is zero at k=0, so artificial zeroing of scaling is ok
            function uhat=apply(fhat)
                kdotf = k1.*fhat(:, 1) + k2.*fhat(:, 2) + k3.*fhat(:, 3);
                uhat = zeros(size(fhat), like=1+1i);
                uhat(:, 1) = Bdiff .* (ksq.*fhat(:, 1) - k1.*kdotf);
                uhat(:, 2) = Bdiff .* (ksq.*fhat(:, 2) - k2.*kdotf);
                uhat(:, 3) = Bdiff .* (ksq.*fhat(:, 3) - k3.*kdotf);
            end
            Dlhat_fun = @apply;
        end
        
        function W0hat_fun = winkernel_fourier(self, k1, k2, k3, sigma_0, Ctrunc)
        % Fourier transform of windowed mollified kernel (i.e. far field)
            ksq = k1.^2 + k2.^2 + k3.^2;
            hf = k1(2)-k1(1); % Should be safe way to get hf
            k = sqrt(ksq);
            xi = 1/sigma_0;
            B = self.fourier_scaling(ksq, xi);
            Bwin = B .* (1 + 1/2*cos(Ctrunc*k) - 3/2*sin(Ctrunc*k)./(Ctrunc*k));
            Bwin(k==0) = 0; % FIX?
            function uhat=apply(fhat)
                kdotf = k1.*fhat(:, 1) + k2.*fhat(:, 2) + k3.*fhat(:, 3);
                uhat = zeros(size(fhat), like=1+1i);
                uhat(:, 1) = Bwin .* (ksq.*fhat(:, 1) - k1.*kdotf);
                uhat(:, 2) = Bwin .* (ksq.*fhat(:, 2) - k2.*kdotf);
                uhat(:, 3) = Bwin .* (ksq.*fhat(:, 3) - k3.*kdotf);
                % Adjust k=0 scaling for windowed B
                % Correcting for weighting in operator here, need to
                % derive what is going on
                uhat(ksq==0, :) = fhat(ksq==0, :) * 2/Ctrunc * (2*pi)^3 / hf^3;
            end
            W0hat_fun = @apply;
        end
    end
end
