classdef laplace_pswf < kernels.SplitKernelInterface
% Laplace kernel split using the orignal Ewald decomposition

    properties (Constant)
        name    = "Laplace PSWF";
        dim_in  = 1;
        dim_out = 1;
    end

    properties
        c_pswf
        gamma_hat
        d2gamma_hat
        Phi
        dPhi
    end

    methods
        function obj = laplace_pswf(args)
        % Constructor, takes either tolerance or sigma_0 as input
            arguments
                args.tolerance (1,1) double = 0 % tolerance default=0 (unset)
                args.c_pswf (1,1) double = 0   % PSWF constant   default=0 (unset)
            end
            obj = obj@kernels.SplitKernelInterface(tolerance=args.tolerance);
            if args.tolerance==0
                if args.c_pswf==0
                    error("Need either tolerance or c_pswf")
                end
                obj.c_pswf = args.c_pswf;
            else
                % TODO: need something better for bandwidth selection
                for c_pswf=1:0.5:60
                    psi = pswf(0, c_pswf);
                    if psi(1) < args.tolerance
                        break
                    end
                end
                obj.c_pswf = c_pswf - 3; % Heuristic
                %obj.c_pswf = -log(args.tolerance / 1.2) + 1; % Heuristic from DKM paper values
                fprintf('[laplace_pswf] auto-selected c_pswf=%g\n', obj.c_pswf);
            end
            % Init PSWF
            obj.Kmax = obj.c_pswf;
            h = harmonic_pswf_split(obj.c_pswf);
            obj.gamma_hat = h.gamma_hat;
            obj.d2gamma_hat = h.d2gamma_hat;
            obj.Phi = h.Phi;
            obj.dPhi = h.dPhi;
        end

        function gamma_hat = fourier_scaling(self, ksq, level)
            rl = 1/2^level;
            kabs = sqrt(ksq);
            mask = (kabs*rl/self.c_pswf) < 1;
            gamma_hat = zeros(size(ksq));
            gamma_hat(mask) = self.gamma_hat(kabs(mask)*rl/self.c_pswf);
        end

        function R = real_decay(self, r, level)
            rl = 1/2^level;
            mask = (r/rl <= 1);
            R = zeros(size(r));
            R(mask) = self.Phi(r(mask)/rl);
        end

    end
    
    methods (Static)
        function u = direct(targets, points, charges)
            N = numel(charges);
            assert(size(points, 1)==N)
            assert(size(points, 2)==3)
            targets = targets.';
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
    end
    
    methods
        function ures = reskernel(self, targets, points, charges, level)
        % Laplace residual kernel R_l(r)
            targets = targets.';
            assert(size(targets, 1)==3);
            r = sqrt( (points(:, 1)-targets(1, :)).^2 + ...
                      (points(:, 2)-targets(2, :)).^2 + ...
                      (points(:, 3)-targets(3, :)).^2 );
            Rl = self.real_decay(r, level)./r;
            Rl(r==0) = 0;
            ures = (charges.' * Rl).';
        end

        function udiff = diffkernel(self, targets, points, charges, level)
        % Laplace difference kernel D_l(r)
            targets = targets.';
            assert(size(targets, 1)==3);
            rl = 1/2^level;
            rlp1 = rl/2;
            r = sqrt( (points(:, 1)-targets(1, :)).^2 + ...
                      (points(:, 2)-targets(2, :)).^2 + ...
                      (points(:, 3)-targets(3, :)).^2 );
            Dl = (-self.real_decay(r, level+1) + self.real_decay(r, level))./r;
            udiff = Dl.' * charges;
        end

        function Mlhat_fun = mollkernel_fourier(self, k1, k2, k3, level)
            ksq = k1.^2 + k2.^2 + k3.^2;
            Mlhat = 4*pi./ksq.*self.fourier_scaling(ksq, level);
            Mlhat(ksq==0) = 0;
            function uhat=apply(fhat)
                uhat = Mlhat .* fhat;
            end
            Mlhat_fun = @apply;
        end

        function Dlhat_fun = diffkernel_fourier(self, k1, k2, k3, level)
        % Fourier transform of Laplace difference kernel
        % Returns operator that applies kernel
            rl = 1/2^level;
            rlp1 = rl/2;
            ksq = k1.^2 + k2.^2 + k3.^2;
            kabs = sqrt(ksq);
            Dlhat = 4*pi./ksq.*(self.fourier_scaling(ksq, level+1) - ...
                                self.fourier_scaling(ksq, level) );
            Dlhat(ksq==0) = 2*pi*self.d2gamma_hat(0)/self.c_pswf^2*(rlp1^2-rl^2);
            function uhat=apply(Psi)
                uhat = Dlhat.*Psi;
            end
            Dlhat_fun = @apply;
        end

        function W0hat_fun = winkernel_fourier(self, k1, k2, k3, Ctrunc)
        % Fourier transform of  windowed mollified Laplace kernel (i.e. far field)
            ksq = k1.^2 + k2.^2 + k3.^2;
            k = sqrt(ksq);
            Hwin = windowed_harmonic(k, Ctrunc);
            W0hat = Hwin .* self.fourier_scaling(ksq, 0);
            function [uhat, const] = apply(fhat)
                const = 0;
                uhat = W0hat .* fhat;
            end
            W0hat_fun = @apply;
        end

        function uself = self_interaction(self, charges, level)
            rl = 1/2^level;
            const = self.dPhi(0)/rl;
            uself = charges*const;
        end
    end
end
