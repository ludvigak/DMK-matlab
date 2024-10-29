classdef stokeslet_pswf2 < kernels.StokesletFourierSplit
% Stokes kernel split using the Hasimoto-like PSWF decomposition

    properties (Constant)
        name    = "Stokes PSWF v2";
        dim_in  = 3;
        dim_out = 3;
    end

    properties
        c_pswf
        c_self;
        pswf_cheb
    end

    methods
        function obj = stokeslet_pswf2(args)
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
                fprintf('[stokeslet_pswf2] auto-selected c_pswf=%g\n', obj.c_pswf);
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
            k = sqrt(ksq);
            psi_arg = k*rl/self.c_pswf;
            gamma_hat = zeros(size(psi_arg));
            supp_mask = psi_arg <= 1; % Truncate to PSWF support
            psi_arg_mask = psi_arg(supp_mask);
            gamma_hat(supp_mask) = psi(psi_arg_mask) - 1/2*psi_arg_mask.*dpsi(psi_arg_mask);
        end
    end
end
