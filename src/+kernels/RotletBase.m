classdef RotletBase < kernels.SplitKernelInterface
% Base class for rotlet splits, contains everything that is not split-specific

% Interface functions to be defined in subclass
    methods (Abstract)
        fourier_decay = fourier_scaling(self, ksq, level)
        % Fourier scaling, corresponding to Fourier transform of screening function
        [Rdiag, Roffd] = real_decay(self, r, level)
        % Diagonal and offdiagonal decay of residual rotlet, scaled by r
    end

    properties (Constant)
        dim_in  = 3;
        dim_out = 3;
    end

    methods (Static)
        function u = direct(targets, sources, f)
        % Rotlet kernel Omega(r)
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
                % r-vectors are Nsrc x Ntrg (of batch)
                r1 = targets(1, range)-sources(:, 1);
                r2 = targets(2, range)-sources(:, 2);
                r3 = targets(3, range)-sources(:, 3);
                r = sqrt(r1.^2 + r2.^2 + r3.^2);
                self_mask = (r==0);
                rinv3 = 1./r.^3;
                rinv3(self_mask) = 0;
                % Source vectors are Nsrc x 1
                f1 = f(:, 1);
                f2 = f(:, 2);
                f3 = f(:, 3);
                % Cross product
                c1 = (f2.*r3 - f3.*r2).*rinv3;
                c2 = (f3.*r1 - f1.*r3).*rinv3;
                c3 = (f1.*r2 - f2.*r1).*rinv3;
                % Sum over targets
                u(range, 1) = sum(c1, 1);
                u(range, 2) = sum(c2, 1);
                u(range, 3) = sum(c3, 1);
            end
        end

        function uhat = fourier_composition(fhat, k1, k2, k3)
            f1 = fhat(:, 1);
            f2 = fhat(:, 2);
            f3 = fhat(:, 3);
            % Cross product
            c1 = (f2.*k3 - f3.*k2);
            c2 = (f3.*k1 - f1.*k3);
            c3 = (f1.*k2 - f2.*k1);
            uhat = [c1 c2 c3];
        end
    end

    methods
        function ures = reskernel(self, targets, sources, f, level)
        % Rotlet residual kernel R_l(r)
            assert(size(sources, 2)==3)
            if size(targets, 1) ~= 3
                targets = targets.';
            end
            assert(size(targets, 1)==3);
            % r-vectors are Nsrc x Ntrg (of batch)
            r1 = targets(1, :)-sources(:, 1);
            r2 = targets(2, :)-sources(:, 2);
            r3 = targets(3, :)-sources(:, 3);
            r = sqrt(r1.^2 + r2.^2 + r3.^2);
            self_mask = (r==0);
            rinv3 = 1./r.^3;
            rinv3(self_mask) = 0;
            % Source vectors are Nsrc x 1
            f1 = f(:, 1);
            f2 = f(:, 2);
            f3 = f(:, 3);
            % Decay
            R = self.real_decay(r, level);
            % Scaled cross product
            c1 = R.*(f2.*r3 - f3.*r2).*rinv3;
            c2 = R.*(f3.*r1 - f1.*r3).*rinv3;
            c3 = R.*(f1.*r2 - f2.*r1).*rinv3;
            % Sum over targets
            ures = zeros(size(targets))';
            ures(:, 1) = sum(c1, 1);
            ures(:, 2) = sum(c2, 1);
            ures(:, 3) = sum(c3, 1);
        end

        function umoll = mollkernel(self, targets, sources, f, level)
        % Mollified kernel M_l = S(r) - R_l(r)
            u = self.direct(targets, sources, f);
            ures = self.reskernel(targets, sources, f, level);
            umoll = u - ures;
        end

        function Mlhat_fun = mollkernel_fourier(self, k1, k2, k3, level)
            ksq = k1.^2 + k2.^2 + k3.^2;
            H = 4*pi./ksq.*self.fourier_scaling(ksq, level);
            H(ksq==0)=0;
            Nf = numel(k1);
            function uhat=apply(fhat)
                assert(all(size(fhat)==[Nf 3]));
                uhat = -1i * H .* self.fourier_composition(fhat, k1, k2, k3);
                % Note k==0 artificially set to zero
            end
            Mlhat_fun = @apply;
        end

        function udiff = diffkernel(self, targets, sources, f, level)
        % Difference kernel D_l(r) = M_{l+1}(r) - M_l(r)
            udiff = mollkernel(self, targets, sources, f, level+1) - ...
                    mollkernel(self, targets, sources, f, level);
        end

        function Dlhat_fun = diffkernel_fourier(self, k1, k2, k3, level)
        % Fourier transform of difference kernel
            ksq = k1.^2 + k2.^2 + k3.^2;
            Hl = self.fourier_scaling(ksq, level);
            Hlp1 = self.fourier_scaling(ksq, level+1);
            Hdiff = 4*pi./ksq.*(Hlp1 - Hl);
            Hdiff(ksq==0) = 0; % Difference kernel is zero at k=0
            function uhat=apply(fhat)
                uhat = -1i * Hdiff .* self.fourier_composition(fhat, k1, k2, k3);
            end
            Dlhat_fun = @apply;
        end

        function W0hat_fun = winkernel_fourier(self, k1, k2, k3, Ctrunc)
        % Fourier transform of windowed mollified kernel (i.e. far field)
            ksq = k1.^2 + k2.^2 + k3.^2;
            k = sqrt(ksq);
            Hwin = windowed_harmonic(k, Ctrunc);
            W0hat = Hwin .* self.fourier_scaling(ksq, 0);
            function [uhat, const]=apply(fhat)
                const = 0;
                uhat = -1i * W0hat .* self.fourier_composition(fhat, k1, k2, k3);
            end
            W0hat_fun = @apply;
        end
    end % End methods
end

