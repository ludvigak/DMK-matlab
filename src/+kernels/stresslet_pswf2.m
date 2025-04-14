classdef stresslet_pswf2 < kernels.StressletFourierSplit
% Stresslet split using the Hasimoto-like PSWF decomposition

    properties (Constant)
        name    = "Stresslet PSWF v2";
    end

    properties
        c_pswf
        c_self;
        gamma_hat;
    end

    methods
        function obj = stresslet_pswf2(args)
        % Constructor, takes either tolerance or sigma_0 as input
            arguments
                args.tolerance (1,1) double = 0 % tolerance default=0 (unset)
                args.c_pswf (1,1) double = 0   % sigma_0   default=0 (unset)
            end
            obj = obj@kernels.StressletFourierSplit(tolerance=args.tolerance);
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
                    if abs(psi(1)) < args.tolerance
                        break
                    end
                end
                obj.c_pswf = c_pswf + 7; % Heuristic (probably too much, but needed for tests)
                                         % c_pswf + 5 is enough to get interpolation working
                fprintf('[stresslet_pswf2] auto-selected c_pswf=%g\n', obj.c_pswf);
            end
            obj.Kmax = obj.c_pswf;
            b = biharmonic_pswf_split(obj.c_pswf);
            obj.gamma_hat = b.gamma_hat;
            obj.Rdiag_cheb = b.r.^2 .* b.d3Bres;
            obj.Roffd_cheb = b.dBres - b.r .* b.d2Bres + 1/3*b.r.^2 .* b.d3Bres;
            obj.c_self = 0;
        end

        function gamma_hat = fourier_scaling(self, ksq, level)
            rl = 1/2^level;
            k = sqrt(ksq);
            psi_arg = k*rl/self.c_pswf;
            supp_mask = psi_arg <= 1; % Truncate to PSWF support
            gamma_hat = zeros(size(psi_arg));
            gamma_hat(supp_mask) = self.gamma_hat(psi_arg(supp_mask));
        end
    end
end
