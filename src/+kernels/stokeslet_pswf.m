classdef stokeslet_pswf < kernels.StokesletBase
% Stokes kernel split using the Hasimoto-like PSWF decomposition

    properties (Constant)
        name    = "Stokes PSWF";
        dim_in  = 3;
        dim_out = 3;
    end

    properties
        c_pswf
        c_self;
        pswf_cheb
        Rdiag_cheb
        Roffd_cheb
    end

    methods
        function obj = stokeslet_pswf(args)
        % Constructor, takes either tolerance or sigma_0 as input
            arguments
                args.tolerance (1,1) double = 0 % tolerance default=0 (unset)
                args.c_pswf (1,1) double = 0   % sigma_0   default=0 (unset)
            end
            obj = obj@kernels.StokesletBase(tolerance=args.tolerance);
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
            end
            % Init PSWF
            psi = pswf(0, obj.c_pswf);
            obj.Kmax = obj.c_pswf;           
            obj.pswf_cheb = psi;
            % Precompute residual decay
            % (enough to do at zero-level since scale-invariant)
            [obj.Rdiag_cheb, obj.Roffd_cheb, obj.c_self] = obj.precompute_real_decay(0);
        end

        function [Sdiag_cheb, Soffd_cheb, c_self] = precompute_real_decay(self, level)
        % Precompute decay of residual function
            rl = 1/2^level;
            % Setup Fourier modes for windowed kernel
            Nk = ceil(self.Kmax/rl);
            k = linspace(0, self.Kmax/rl, Nk).';
            hk = k(2)-k(1);
            R = sqrt(3)+1;
            % Explicitly duplicate the windowing here, since corrections needed
            Bwin = 8*pi*((1 + 1/2*cos(R*k) - 3/2*sin(R*k)./(R*k)))./k.^4;
            Bwin(k==0) = pi*R^4/15;
            gamma_hat = fourier_scaling(self, k.^2, level);
            Bmollhat = Bwin .* gamma_hat;
            % Setup the integrands that correspond to the radial Fourier transforms,
            % and their derivatives w.r.t. radius r
            %f = @(k,r)  k./r   .*(sin(k.*r)); % Not needed, kept for completeness
            function y=df(k,r)
                assert(size(r, 1)==1)
                y = -k./r.^2.*(sin(k.*r)-k.*r.*cos(k.*r));
                if any(r==0)
                    y(:,r==0) = 0;
                end
            end
            function y=d2f(k,r)
                assert(size(r, 1)==1)
                y = -k./r.^3.*(sin(k.*r).*-2.0+k.^2.*r.^2.*sin(k.*r)+k.*r.*cos(k.*r).*2.0);
                if any(r==0)
                    y(:, r==0) = -k.^4/3;
                end
            end
            % Fourier integrals for evaluating mollified biharmonic
            %Bmoll   = @(r) -4*pi*hk*sum(  f(k,r) .* Bmollhat, 1) / (2*pi)^3;
            dBmoll  = @(r) -4*pi*hk*sum( df(k,r) .* Bmollhat, 1) / (2*pi)^3;
            d2Bmoll = @(r) -4*pi*hk*sum(d2f(k,r) .* Bmollhat, 1) / (2*pi)^3;
            %  Residual biharmonic, including corrections
            %Bres   = @(r) r - r.^2/2/R - R/2 - Bmoll(r);
            dBres  = @(r) 1 - r/R - dBmoll(r);
            d2Bres = @(r) - 1/R - d2Bmoll(r);
            % Construct chebfuns of diagonal and offdiagonal
            Sdiag_cheb = chebfun(@(r)  r .* d2Bres(r) + dBres(r), [0 rl]);
            Soffd_cheb = chebfun(@(r) -r .* d2Bres(r) + dBres(r), [0 rl]);
            % Self interaction
            yself = -(2*k.^4)/3; % Limit d2f(r)-df(r)/r as r->0;
            c_self = -4*pi*hk*sum( yself .* Bmollhat, 1) / (2*pi)^3 + 2/R;
        end        
    end

    methods
        function gamma_hat = fourier_scaling(self, ksq, level)
            rl = 1/2^level;
            psi = self.pswf_cheb;
            dpsi = diff(psi);
            d2psi = diff(dpsi);
            alpha = -d2psi(0)/psi(0)*rl^2/self.c_pswf^2/2;
            k = sqrt(ksq);
            psi_arg = k*rl/self.c_pswf;
            gamma_hat = zeros(size(psi_arg));
            supp_mask = psi_arg <= 1; % Truncate to PSWF support
            gamma_hat(supp_mask) = psi(psi_arg(supp_mask))/psi(0) .* (1 + alpha*ksq(supp_mask));
        end

        function [Rdiag, Roffd] = real_decay(self, r, level)
            rl = 1/2^level;
            % Rescale to 0-level, where decay is precomputed
            r = r/rl;
            [Rdiag, Roffd] = deal(zeros(size(r)));
            supp_mask = r <= 1; % Truncate to PSWF support
            Rdiag(supp_mask) = self.Rdiag_cheb(r(supp_mask));
            Roffd(supp_mask) = self.Roffd_cheb(r(supp_mask));
        end

        function uself = self_interaction(self, charges, level)
            rl = 1/2^level;
            % Scale correct since c_self computed at 0-level
            uself = -self.c_self/rl * charges;
        end        
    end
end
