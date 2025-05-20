classdef rotlet_pswf < kernels.RotletBase
% Stokes kernel split using the Pswf decomposition

    properties (Constant)
        name    = "Rotlet PSWF";
    end

    properties
        p
        sigma_0
        c_pswf
        gamma_hat
        R_cheb
    end

    methods
        function obj = rotlet_pswf(args)
        % Constructor, takes either tolerance or sigma_0 as input
            arguments
                args.tolerance (1,1) double = 0 % tolerance default=0 (unset)
                args.c_pswf (1,1) double = 0    % PSWF constant   default=0 (unset)
            end
            obj = obj@kernels.RotletBase(tolerance=args.tolerance);
            if args.tolerance==0
                if args.c_pswf==0
                    error("Need either tolerance or c_pswf")
                end
                obj.c_pswf = args.c_pswf;
            else
                % Autoselect PSWF bandwidth
                obj.c_pswf = pi/3*ceil(3/pi*(log10(args.tolerance)+0.40) / -0.42);
                fprintf('[stokeslet_pswf] auto-selected c_pswf=%g\n', obj.c_pswf);
            end
            obj.p = ceil(1.44*obj.c_pswf - 1.95);
            % Init PSWF
            obj.Kmax = obj.c_pswf;
            h = harmonic_pswf_split(obj.c_pswf);
            obj.R_cheb = h.Phi - h.r*h.dPhi;
            obj.gamma_hat = h.gamma_hat;
        end
    end

    methods
        function H = fourier_scaling(self, ksq, level)
            rl = 1/2^level;
            kabs = sqrt(ksq);
            H = zeros(size(kabs));
            mask = (kabs*rl/self.c_pswf <= 1);
            H(mask) = self.gamma_hat(kabs(mask)*rl/self.c_pswf);
        end

        function R = real_decay(self, r, level)
            rl = 1/2^level;
            mask = (r/rl <= 1);
            R = zeros(size(r));
            R(mask) = self.R_cheb(r(mask)/rl);
        end

        function uself = self_interaction(self, charges, level)
            uself = zeros(size(charges));
        end
    end
end
