classdef stresslet_hasimoto < kernels.StressletBase
% Stresslet kernel split using the Hasimoto decomposition

    properties (Constant)
        name    = "Stresslet Hasimoto";
    end

    properties
        sigma_0
    end

    methods
        function obj = stresslet_hasimoto(args)
        % Constructor, takes either tolerance or sigma_0 as input
            arguments
                args.tolerance (1,1) double = 0 % tolerance default=0 (unset)
                args.sigma_0 (1,1) double = 0   % sigma_0   default=0 (unset)
            end
            obj = obj@kernels.StressletBase(tolerance=args.tolerance);
            if args.tolerance==0
                if args.sigma_0==0
                    error("Need either tolerance or sigma_0")
                end
                obj.sigma_0 = args.sigma_0;
            else
                obj.sigma_0 = 1/sqrt(log(1/obj.tolerance))*0.9; % TODO
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

        function [Rdiag, Roffd] = real_decay(self, r, level)
            sigma_l = self.sigma_level(level);
            xi = 1/sigma_l;
            % Hasimoto
            r2 = r.^2;
            c = xi^2*r2;
            Roffd = (2/6)*( 3.0*erfc(xi*r) + 2.0*xi/sqrt(pi)*r.*(3.0+2.0*c).*exp(-c) );
            Rdiag = 4/sqrt(pi)*xi^3.*exp(-c);
        end

        function uself = self_interaction(self, charges, level)
            sigma_l = self.sigma_level(level);
            uself = 0; % TODO
        end
    end
end
