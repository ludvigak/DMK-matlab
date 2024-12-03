classdef StressletBase < kernels.SplitKernelInterface
% Base class for Stresslet splits, contains everything that is not split-specific
% Input f is packed as [q1 q2 q3 n1 n2 n3]

% Interface functions to be defined in subclass
    methods (Abstract)        
        fourier_decay = fourier_scaling(self, ksq, level)
        % Fourier scaling, corresponding to Fourier transform of screening function
        [Rdiag, Roffd] = real_decay(self, r, level)
        % Diagonal and offdiagonal decay of residual stresslet, scaled by r
    end

    properties (Constant)
        dim_in  = 6;
        dim_out = 3;
    end

    methods (Static)
        function u = direct(targets, sources, f)
        % Stresslet kernel T(r)
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
                rinv5 = r.^(-5);
                rinv5(self_mask) = 0;
                % Source vectors are Nsrc x 1
                q1 = f(:, 1);
                q2 = f(:, 2);
                q3 = f(:, 3);
                n1 = f(:, 4);
                n2 = f(:, 5);
                n3 = f(:, 6);
                rdotq = q1.*r1 + q2.*r2  + q3.*r3; % Nsrc x Ntrg
                rdotn = n1.*r1 + n2.*r2  + n3.*r3; % Nsrc x Ntrg
                fact = -6*rdotq.*rdotn.*rinv5; % Nsrc x Ntrg
                u(range, :) = [sum(r1.*fact, 1) % 1 x Ntrg
                               sum(r2.*fact, 1)
                               sum(r3.*fact, 1)].'; % Ntrg x 3
            end
        end

        function uhat = fourier_composition(fhat, k1, k2, k3, ksq)
            N = size(fhat, 1);
            if ndims(fhat)==2 && all(size(fhat)==[N 9])
                fhat = reshape(fhat, N, 3, 3);
            end
            assert(ndims(fhat)==3 && all(size(fhat)==[N 3 3]))
            uhat = zeros(N, 3, like=1+1i);

            if false
                % Explicit code:
                kvec = {k1, k2, k3};
                for l=1:3
                    for m=1:3
                        for n=1:3
                            uhat(:, l) = uhat(:, l) + 1i*( ...
                                ksq.*(kvec{l}*(m==n) + kvec{m}*(n==l) + kvec{n}*(l==m)) + ...
                                -2*kvec{l}.*kvec{m}.*kvec{n} ...
                                                         ) .* fhat(:, m, n);
                        end
                    end
                end
            else
                % Slightly operation-optimized code:
                double_product = ...
                    + k1.*k1.*fhat(:, 1, 1) ...
                    + k1.*k2.*fhat(:, 1, 2) ...
                    + k1.*k3.*fhat(:, 1, 3) ...
                    + k2.*k1.*fhat(:, 2, 1) ...
                    + k2.*k2.*fhat(:, 2, 2) ...
                    + k2.*k3.*fhat(:, 2, 3) ...
                    + k3.*k1.*fhat(:, 3, 1) ...
                    + k3.*k2.*fhat(:, 3, 2) ...
                    + k3.*k3.*fhat(:, 3, 3);
                fhat_diag = fhat(:, 1, 1) + fhat(:, 2, 2) + fhat(:, 3, 3);
                product_1 = fhat(:, :, 1).*k1 + fhat(:, :, 2).*k2 + fhat(:, :, 3).*k3;
                product_2 = fhat(:, 1, :).*k1 + fhat(:, 2, :).*k2 + fhat(:, 3, :).*k3;
                product_2 = squeeze(product_2);
                products = product_1 + product_2;
                uhat = 1i*([k1 k2 k3].*(ksq.*fhat_diag - 2*double_product) + ksq.*products);
            end
        end

        function fprod = input_product(f)
            N = size(f, 1);
            fprod = zeros(N, 3, 3);
            for i1=1:3
                for i2=1:3
                    fprod(:, i1, i2) = f(:, i1) .* f(:,3+i2);
                end
            end
            fprod = reshape(fprod, N, 9);
        end

    end

    methods
        function ures = reskernel(self, targets, sources, f, level)
        % Stresslet residual kernel R_l(r)
            assert(size(sources, 2)==3)
            if size(targets, 1) ~= 3
                targets = targets.';
            end
            assert(size(targets, 1)==3);
            % r-vectors are Nsrc x Ntrg
            r1 = targets(1, :)-sources(:, 1);
            r2 = targets(2, :)-sources(:, 2);
            r3 = targets(3, :)-sources(:, 3);
            r = sqrt(r1.^2 + r2.^2 + r3.^2);
            self_mask = (r==0);
            rinv = 1./r;
            rinv(self_mask) = 0;
            rinv3 = rinv.^3;
            rinv5 = rinv.^5;
            % Source vectors are Nsrc x 1
            q1 = f(:, 1);
            q2 = f(:, 2);
            q3 = f(:, 3);
            n1 = f(:, 4);
            n2 = f(:, 5);
            n3 = f(:, 6);
            rdotq = q1.*r1 + q2.*r2  + q3.*r3; % Nsrc x Ntrg
            rdotn = n1.*r1 + n2.*r2  + n3.*r3; % Nsrc x Ntrg
            qdotn = n1.*q1 + n2.*q2  + n3.*q3; % Nsrc x 1
            [Rdiag, Roffd] = self.real_decay(r, level);
            % Diagonal term
            diag_fact  = Rdiag.*rinv3; % common factor
            udiag = [ sum( (r1.*qdotn + q1.*rdotn + n1.*rdotq).*diag_fact, 1)    % 1 x Ntrg
                      sum( (r2.*qdotn + q2.*rdotn + n2.*rdotq).*diag_fact, 1)
                      sum( (r3.*qdotn + q3.*rdotn + n3.*rdotq).*diag_fact, 1) ]; % 3 x Ntrg
            % Off-diagonal term
            offd_fact = -6*rdotq.*rdotn.*rinv5.*Roffd; % common factor
            uoffd = [sum(r1.*offd_fact, 1)   % 1 x Ntrg
                     sum(r2.*offd_fact, 1)
                     sum(r3.*offd_fact, 1)]; % 3 x Ntrg
            ures = (udiag + uoffd).';
        end

        function umoll = mollkernel(self, targets, sources, f, level)
        % Mollified kernel M_l = S(r) - R_l(r)
            u = self.direct(targets, sources, f);
            ures = self.reskernel(targets, sources, f, level);
            umoll = u - ures;
        end

        function Mlhat_fun = mollkernel_fourier(self, k1, k2, k3, level)
            ksq = k1.^2 + k2.^2 + k3.^2;
            B = -8*pi./ksq.^2.*self.fourier_scaling(ksq, level);
            B(ksq==0)=0;
            Nf = numel(k1);
            function uhat=apply(fhat)
                assert(all(size(fhat)==[Nf 3 3]));
                uhat = -B .* self.fourier_composition(fhat, k1, k2, k3, ksq);
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
            Bl = self.fourier_scaling(ksq, level);
            Blp1 = self.fourier_scaling(ksq, level+1);
            Bdiff = 8*pi./ksq.^2.*(Blp1 - Bl);
            Bdiff(ksq==0) = 0; % Difference kernel is zero at k=0
            function uhat=apply(fhat)
                uhat = Bdiff .* self.fourier_composition(fhat, k1, k2, k3, ksq);
            end
            Dlhat_fun = @apply;
        end

        function W0hat_fun = winkernel_fourier(self, k1, k2, k3, Ctrunc)
        % Fourier transform of windowed mollified kernel (i.e. far field)
            ksq = k1.^2 + k2.^2 + k3.^2;
            k = sqrt(ksq);
            Bwin = windowed_biharmonic(k, Ctrunc);
            Bwin = -Bwin .* self.fourier_scaling(ksq, 0);
            function [uhat, const]=apply(fhat)
            % const = constant term to be added to solution
                uhat = self.fourier_composition(fhat, k1, k2, k3, ksq);
                uhat = Bwin.*uhat;
                % No constant term needed for stresslet
                const = 0;
            end
            W0hat_fun = @apply;
        end
    end % End methods
end

