classdef StokesletFourierSplit < kernels.StokesletBase
% Base class for Stokeslet split using Fourier decomposition
%
% This contains routines for precomputing and evaluting the real space
% and self-interaction components, given a Fourier split
    
    properties
        Rdiag_cheb
        Roffd_cheb
        % chebfuns of the diagonal and off-diagonal components of the residual kernel
    end
    
    methods
        function [Rdiag, Roffd] = real_decay(self, r, level)
            rl = 1/2^level;
            % Rescale to 0-level, where decay is precomputed
            r = r/rl;
            [Rdiag, Roffd] = deal(zeros(size(r)));
            supp_mask = r <= 1; % Truncate to chebfun support
            Rdiag(supp_mask) = self.Rdiag_cheb(r(supp_mask));
            Roffd(supp_mask) = self.Roffd_cheb(r(supp_mask));
        end

        function uself = self_interaction(self, charges, level)
            rl = 1/2^level;
            % Scale correct since c_self computed at 0-level
            uself = -self.c_self/rl * charges;
        end        

        function [Rdiag_cheb, Roffd_cheb, c_self] = precompute_real_decay(self, level)
        % Precompute decay of residual function, given a Fourier split
            rl = 1/2^level;
            % Setup Fourier modes for windowed kernel
            Nk = ceil(self.Kmax/rl);
            k = linspace(0, self.Kmax/rl, Nk).';
            hk = k(2)-k(1);
            R = sqrt(3)+1;
            % Explicitly duplicate the windowing here, since corrections needed
            Bwin = -8*pi*((1 + 1/2*cos(R*k) - 3/2*sin(R*k)./(R*k)))./k.^4;
            Bwin(k==0) = -pi*R^4/15;
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
            %Bmoll   = @(r) 4*pi*hk*sum(  f(k,r) .* Bmollhat, 1) / (2*pi)^3;
            dBmoll  = @(r) 4*pi*hk*sum( df(k,r) .* Bmollhat, 1) / (2*pi)^3;
            d2Bmoll = @(r) 4*pi*hk*sum(d2f(k,r) .* Bmollhat, 1) / (2*pi)^3;
            %  Residual biharmonic, including corrections
            %Bres   = @(r) r - r.^2/2/R - R/2 - Bmoll(r);
            dBres  = @(r) 1 - r/R - dBmoll(r);
            d2Bres = @(r) - 1/R - d2Bmoll(r);
            % Construct chebfuns of diagonal and offdiagonal
            Rdiag_cheb = chebfun(@(r)  r .* d2Bres(r) + dBres(r), [0 rl]);
            Roffd_cheb = chebfun(@(r) -r .* d2Bres(r) + dBres(r), [0 rl]);
            % Self interaction
            yself = -(2*k.^4)/3; % Limit df(r)/r+d2f(r) as r->0;
            c_self = 4*pi*hk*sum( yself .* Bmollhat, 1) / (2*pi)^3 + 2/R;
        end
    end
end

