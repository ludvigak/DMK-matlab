classdef stokeslet_hasimoto < kernels.StokesletBase
% Stokes kernel split using the Hasimoto decomposition

    properties (Constant)
        name    = "Stokes Hasimoto";
        dim_in  = 3;
        dim_out = 3;
    end

    properties
        sigma_0
    end

    methods
        function obj = stokeslet_hasimoto(args)
        % Constructor, takes either tolerance or sigma_0 as input
            arguments
                args.tolerance (1,1) double = 0 % tolerance default=0 (unset)
                args.sigma_0 (1,1) double = 0   % sigma_0   default=0 (unset)
            end
            obj = obj@kernels.StokesletBase(tolerance=args.tolerance);
            if args.tolerance==0
                if args.sigma_0==0
                    error("Need either tolerance or sigma_0")
                end
                obj.sigma_0 = args.sigma_0;
            else
                obj.sigma_0 = 1/sqrt(log(1/obj.tolerance)); % TODO: Good choice also here?
            end
            obj.Kmax = 2/obj.sigma_0^2;
        end

        function sigma_l = sigma_level(self, level)
            sigma_l = self.sigma_0 / 2^level;
        end
    end

    methods
        function B = fourier_scaling(self, ksq, level)
            xi = 1/self.sigma_level(level);
            B = (1 + ksq/(4*xi^2)) .*  exp(-ksq/4/xi^2);
        end

        function [Sdiag, Soffd] = real_decay(self, r, level)
            sigma_l = self.sigma_level(level);
            xi = 1/sigma_l;
            Sdiag = -2*xi/sqrt(pi)*exp(-xi^2*r.^2).*r + erfc(xi*r);
            Soffd =  2*xi/sqrt(pi)*exp(-xi^2*r.^2).*r + erfc(xi*r);
        end

        function uself = self_interaction(self, charges, level)
            sigma_l = self.sigma_level(level);
            uself = -4/(sigma_l*sqrt(pi)) * charges;
        end        
    end
end
