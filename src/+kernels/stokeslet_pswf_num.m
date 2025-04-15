classdef stokeslet_pswf_num < kernels.StokesletFourierSplit
% Stokes kernel split using the Hasimoto-like PSWF decomposition (the bad one)
% and numerically computed residual kernel.
    
    properties (Constant)
        name    = "Stokes PSWF (num)";
        dim_in  = 3;
        dim_out = 3;
    end

    properties
        c_pswf
        c_self;
        pswf_cheb
    end

    methods
        function obj = stokeslet_pswf_num(args)
        % Constructor, takes either tolerance or sigma_0 as input
            arguments
                args.tolerance (1,1) double = 0 % tolerance default=0 (unset)
                args.c_pswf (1,1) double = 0   % sigma_0   default=0 (unset)
            end
            obj = obj@kernels.StokesletFourierSplit(tolerance=args.tolerance);
            if args.tolerance==0
                if args.c_pswf==0
                    error("Need either tolerance or c_pswf")
                end
                obj.c_pswf = args.c_pswf;
            else
                % Autoselect PSWF bandwidth
                % TODO: need something better
                for c_pswf=1:0.5:60
                    psi = pswf(0, c_pswf);
                    if psi(1) < args.tolerance
                        break
                    end
                end
                obj.c_pswf = c_pswf + 4; % Heuristic
                fprintf('[stokeslet_pswf_num] auto-selected c_pswf=%g\n', obj.c_pswf);
            end
            % Init PSWF
            psi = pswf(0, obj.c_pswf);
            obj.Kmax = obj.c_pswf;           
            obj.pswf_cheb = psi;
            % Precompute residual decay
            % (enough to do at zero-level since scale-invariant)
            [obj.Rdiag_cheb, obj.Roffd_cheb, obj.c_self] = obj.precompute_real_decay(0);
        end

        function gamma_hat = fourier_scaling(self, ksq, level)
            rl = 1/2^level;
            psi = self.pswf_cheb;
            psi = psi/psi(0);
            dpsi = diff(psi);
            d2psi = diff(dpsi);
            alpha = -d2psi(0)*rl^2/self.c_pswf^2/2;
            k = sqrt(ksq);
            psi_arg = k*rl/self.c_pswf;
            gamma_hat = zeros(size(psi_arg));
            supp_mask = psi_arg <= 1; % Truncate to PSWF support
            gamma_hat(supp_mask) = psi(psi_arg(supp_mask)) .* (1 + alpha*ksq(supp_mask));
        end
    end
end
