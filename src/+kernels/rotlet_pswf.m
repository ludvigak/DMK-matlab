classdef rotlet_pswf < kernels.RotletBase
% Stokes kernel split using the Pswf decomposition

    properties (Constant)
        name    = "Rotlet PSWF";
    end

    properties
        sigma_0
        c_pswf
        pswf_cheb
        pswf_c0
        pswf_erfc_cheb
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
                % TODO: need something better for bandwidth selection
                for c_pswf=1:0.5:60
                    psi = pswf(0, c_pswf);
                    if psi(1) < args.tolerance
                        break
                    end
                end
                obj.c_pswf = c_pswf + 1; % Heuristic
                fprintf('[rotlet_pswf] auto-selected c_pswf=%g\n', obj.c_pswf);
            end
            % Init PSWF
            obj.Kmax = obj.c_pswf;
            obj.pswf_cheb = pswf(0, obj.c_pswf);
            obj.pswf_c0 = integral(obj.pswf_cheb, 0, 1);
            pswf_erfc = @(z) integral(obj.pswf_cheb, z, 1) / obj.pswf_c0;
            obj.pswf_erfc_cheb = chebfun(pswf_erfc, [0 1]);
        end

        function y = psi(self, x)
            y = zeros(size(x));
            y(x <= 1) = self.pswf_cheb(x(x <= 1));
        end

        function y = pswf_erfc(self, x)
            y = zeros(size(x));
            mask = (x<1);
            xm = x(mask);
            y(mask) = self.pswf_erfc_cheb(xm);
        end

    end

    methods
        function B = fourier_scaling(self, ksq, level)
            rl = 1/2^level;
            kabs = sqrt(ksq);
            B = self.psi(kabs*rl/self.c_pswf)/self.psi(0);
        end

        function R = real_decay(self, r, level)
            rl = 1/2^level;
            R =  r.*self.psi(r/rl)/(rl*self.pswf_c0) + self.pswf_erfc(r/rl);
        end

        function uself = self_interaction(self, charges, level)
            uself = zeros(size(charges));
        end
    end
end
