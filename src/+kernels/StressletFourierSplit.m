classdef StressletFourierSplit < kernels.StressletBase
% Base class for stresslet split using Fourier decomposition
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
            Ncharge = size(charges, 1);
            % Unecessary code, this is all zero for all splits (?)
            uself = -self.c_self/rl * ones(Ncharge, 3);
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
            % Setup the integrands that correspond to the radial Fourier transforms,
            % and their derivatives w.r.t. radius r
            % === Generating code:
            % syms k r; f = k/r*sin(k*r);
            % for d=0:3
            %     fprintf('d%df     = %s\n', d, func2str(matlabFunction(simplify(diff(f,d)))))
            %     fprintf('d%df_lim = %s\n', d, func2str(matlabFunction(limit(diff(f,d)), r, 0)))
            % end
            % ==================
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
            function y=d3f(k,r)
                assert(size(r, 1)==1)
                y = -k.*1.0./r.^4.*(sin(k.*r).*6.0+k.^3.*r.^3.*cos(k.*r)-k.^2.*r.^2.*sin(k.*r).*3.0-k.*r.*cos(k.*r).*6.0);
                if any(r==0)
                    y(:, r==0) = 0;
                end
            end
            % Fourier integrals for evaluating mollified biharmonic
            %Bmoll   = @(r) -4*pi*hk*sum(  f(k,r) .* Bmollhat, 1) / (2*pi)^3;
            dBmoll  = @(r) 4*pi*hk*sum( df(k,r) .* Bmollhat, 1) / (2*pi)^3;
            d2Bmoll = @(r) 4*pi*hk*sum(d2f(k,r) .* Bmollhat, 1) / (2*pi)^3;
            d3Bmoll = @(r) 4*pi*hk*sum(d3f(k,r) .* Bmollhat, 1) / (2*pi)^3;
            %  Residual biharmonic, including corrections
            %Bres   = @(r) r + alpha + beta*r.^2 - Bmoll(r);
            dBres  = @(r) 1 + 2*beta*r - dBmoll(r);
            d2Bres = @(r) 2*beta - d2Bmoll(r);
            d3Bres = @(r) - d3Bmoll(r);
            % Construct chebfuns of diagonal and offdiagonal
            Rdiag_cheb = chebfun(@(r) r.^2.*d3Bres(r), [0 rl]);
            Roffd_cheb = chebfun(@(r) dBres(r) - r .* d2Bres(r) + 1/3*r.^2 .* d3Bres(r), [0 rl]);
            % Self interaciton is zero
            c_self = 0;
        end 
    end
end

