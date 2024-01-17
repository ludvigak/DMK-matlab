classdef laplace_pswf < kernels.SplitKernelInterface
% Laplace kernel split using the orignal Ewald decomposition

    properties (Constant)
        name    = "Laplace PSWF";
        dim_in  = 1;
        dim_out = 1;
    end

    properties
        c_pswf
        pswf_cheb
        pswf_c0
        pswf_c2
        pswf_erf_cheb
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
                obj.c_pswf = -log(args.tolerance / 1.2) + 1; % Heuristic from DKM paper values
            end
            % Init PSWF
            obj.Kmax = obj.c_pswf;
            obj.pswf_cheb = pswf(0, obj.c_pswf);
            obj.pswf_c0 = integral(obj.pswf_cheb, 0, 1);
            x = chebfun(@(x) x);
            obj.pswf_c2 = integral(x*x*obj.pswf_cheb, 0, 1);
            pswf_erf = @(z) integral(obj.pswf_cheb, 0, z) / obj.pswf_c0;
            obj.pswf_erf_cheb = chebfun(pswf_erf, [0 1]);
        end

        function y = psi(self, x)
            y = zeros(size(x));
            y(x <= 1) = self.pswf_cheb(x(x <= 1));
        end

        function y = pswf_erf(self, x)
            y = ones(size(x));
            y(x < 1) = self.pswf_erf_cheb(x(x < 1));
        end

        function y = pswf_erfc(self, x)
            y = 1 - self.pswf_erf(x);
        end
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
    end
    
    methods
        function ures = reskernel(self, targets, points, charges, level)
        % Laplace residual kernel R_l(r)
            rl = 1/2^level;
            if size(targets, 1) ~= 3
                targets = targets.';
            end
            assert(size(targets, 1)==3);
            r = sqrt( (points(:, 1)-targets(1, :)).^2 + ...
                      (points(:, 2)-targets(2, :)).^2 + ...
                      (points(:, 3)-targets(3, :)).^2 );
            Rl = self.pswf_erfc(r/rl)./r;
            Rl(r==0) = 0;
            ures = (charges.' * Rl).';
        end

        function udiff = diffkernel(self, targets, points, charges, level)
        % Laplace difference kernel D_l(r)
            if size(targets, 1) ~= 3
                targets = targets.';
            end
            assert(size(targets, 1)==3);
            rl = 1/2^level;
            rlp1 = rl/2;
            r = sqrt( (points(:, 1)-targets(1, :)).^2 + ...
                      (points(:, 2)-targets(2, :)).^2 + ...
                      (points(:, 3)-targets(3, :)).^2 );
            Dl = (self.pswf_erf(r/rlp1) - self.pswf_erf(r/rl))./r;
            udiff = Dl.' * charges;
        end
        
        function Dlhat_fun = diffkernel_fourier(self, k1, k2, k3, level)
        % Fourier transform of Laplace difference kernel
        % Returns operator that applies kernel
            rl = 1/2^level;
            rlp1 = rl/2;
            ksq = k1.^2 + k2.^2 + k3.^2;
            kabs = sqrt(ksq);
            Dlhat = 4*pi*(self.psi(kabs*rlp1/self.c_pswf) ...
                          - self.psi(kabs*rl/self.c_pswf))./ksq/self.psi(0);
            Dlhat(ksq==0) = 2*pi*self.pswf_c2/self.pswf_c0*(rl^2-rlp1^2);
            function uhat=apply(Psi)
                uhat = Dlhat.*Psi;
            end
            Dlhat_fun = @apply;
        end

        function W0hat_fun = winkernel_fourier(self, k1, k2, k3, Ctrunc)
        % Fourier transform of  windowed mollified Laplace kernel (i.e. far field)
            ksq = k1.^2 + k2.^2 + k3.^2;
            k = sqrt(ksq);
            W0hat = 8*pi*(sin(Ctrunc*k/2)./k).^2 .* self.psi(k/self.c_pswf)/self.psi(0);
            W0hat(ksq==0) = 8*pi * Ctrunc^2/4;
            function uhat=apply(fhat)
                uhat = W0hat .* fhat;
            end
            W0hat_fun = @apply;
        end

        function uself = self_interaction(self, charges, level)
            rl = 1/2^level;
            const = -self.psi(0)/(self.pswf_c0*rl);
            uself = charges*const;
        end
    end
end
