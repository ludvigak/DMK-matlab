classdef stokeslet_exp4 < kernels.StokesletFourierSplit
% Stokes kernel split using the exp(-k^4) decay

    properties (Constant)
        name    = "Stokes Exp(-k^4)";
        dim_in  = 3;
        dim_out = 3;
    end

    properties
        sigma_0;
        c_self;
    end

    methods
        function obj = stokeslet_exp4(args)
        % Constructor, takes either tolerance or sigma_0 as input
            arguments
                args.tolerance (1,1) double = 0 % tolerance default=0 (unset)
            end
            obj = obj@kernels.StokesletFourierSplit(tolerance=args.tolerance);
            if args.tolerance==0
                error("Need tolerance")
            end
            % Parameter selection = guesswork
            obj.sigma_0 = 0.9/log(1/obj.tolerance);
            obj.Kmax = log(1/obj.tolerance)^(1/4)/obj.sigma_0;
            % Precompute residual decay
            % (enough to do at zero-level since scale-invariant)
            [obj.Rdiag_cheb, obj.Roffd_cheb, obj.c_self] = obj.precompute_real_decay(0);
        end
        
        function gamma_hat = fourier_scaling(self, ksq, level)
            sigma_l = self.sigma_0 / 2^level;
            gamma_hat = exp(-ksq.^2 * sigma_l^4);
        end
    end
end
