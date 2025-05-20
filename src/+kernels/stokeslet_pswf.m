classdef stokeslet_pswf < kernels.StokesletFourierSplit
% Stokes kernel split using the Hasimoto-like PSWF decomposition

    properties (Constant)
        name    = "Stokeslet PSWF";
        dim_in  = 3;
        dim_out = 3;
    end

    properties
        p
        c_pswf
        c_self;
        gamma_hat
    end

    methods
        function obj = stokeslet_pswf(args)
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
                obj.c_pswf = pi/3*ceil(3/pi*(log10(args.tolerance)-1.11) / -0.41);
                fprintf('[stokeslet_pswf] auto-selected c_pswf=%g\n', obj.c_pswf);
            end
            obj.p = ceil(1.43*obj.c_pswf - 2.76);
            obj.Kmax = obj.c_pswf;
            b = biharmonic_pswf_split(obj.c_pswf);
            obj.Rdiag_cheb =  b.r .* b.d2Bres + b.dBres;
            obj.Roffd_cheb = -b.r .* b.d2Bres + b.dBres;
            obj.c_self = 2*b.d2Bmoll(0);
            obj.gamma_hat = b.gamma_hat;
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
