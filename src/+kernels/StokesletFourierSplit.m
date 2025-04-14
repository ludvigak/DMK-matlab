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
            [Bwin, alpha, beta] = windowed_biharmonic(k, R);
            gamma_hat = fourier_scaling(self, k.^2, level);
            Bmollhat = Bwin .* gamma_hat;
            % Fourier integrals for evaluating mollified biharmonic
            if true
                % See comment in StressletFourierSplit.m
                [f, df, d2f] = radial_fourier_kernels_gen();
                %Bmoll   = @(r) 4*pi*hk*sum(  f(k,r) .* Bmollhat, 1) / (2*pi)^3;
                dBmoll  = @(r) 4*pi*hk*sum( df(k,r) .* Bmollhat, 1) / (2*pi)^3;
                d2Bmoll = @(r) 4*pi*hk*sum(d2f(k,r) .* Bmollhat, 1) / (2*pi)^3;
                %  Residual biharmonic, including corrections
                %Bres   = @(r) r + alpha + beta*r.^2 - Bmoll(r);
                dBres  = @(r) 1 + 2*beta*r - dBmoll(r);
                d2Bres = @(r) 2*beta - d2Bmoll(r);
                % Construct chebfuns of diagonal and offdiagonal
                Rdiag_cheb = chebfun(@(r)  r .* d2Bres(r) + dBres(r), [0 rl]);
                Roffd_cheb = chebfun(@(r) -r .* d2Bres(r) + dBres(r), [0 rl]);
            else
                f = @(k,r)  k./r   .*(sin(k.*r));
                Bmoll = @(r) 4*pi*hk*sum(  f(k,r) .* Bmollhat, 1) / (2*pi)^3;
                Bres = @(r) r + alpha + beta*r.^2 - Bmoll(r);
                Bres_cheb = chebfun(Bres, [0 rl], 'turbo');
                r_cheb = chebfun(@(r) r, [0 rl]);
                dBres_cheb = diff(Bres_cheb);
                d2Bres_cheb = diff(Bres_cheb, 2);
                d3Bres_cheb = diff(Bres_cheb, 3);
                Rdiag_cheb = r_cheb .* d2Bres_cheb + dBres_cheb;
                Roffd_cheb = -r_cheb .* d2Bres_cheb + dBres_cheb;
            end
            % Self interaction
            yself = -(2*k.^4)/3; % Limit df(r)/r+d2f(r) as r->0;
            c_self = 4*pi*hk*sum( yself .* Bmollhat, 1) / (2*pi)^3 - 4*beta;
        end
    end
end

