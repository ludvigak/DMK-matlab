%% Main function to generate tests
function tests = test_planewave
tests = functiontests(localfunctions);
end

%% Test Functions


function test_outgoing(testCase)
    tol = 1e-9;
    p = 45;
    max_level = 3;
    sigma_0 = 1/sqrt(log(1/tol));
    r0 = 1;
    D = 3*r0;
    h0 = 2*pi/D;
    K0 = 4/1 * log(1/tol);
    nf = ceil(K0/h0);
    Tprox2pw = operator_proxy2planewave(p, h0, nf, max_level);
    Tpw2poly = operator_planewave2local(p, h0, nf, max_level);
    for l = 0:max_level
        rl = 1/2^l;
        sigma_l = sigma_0 / 2^l;
        sigma_lp1 = sigma_l / 2;
        hl = h0 / rl;
        Kl = 4/rl * log(1/tol);
        nf = ceil(Kl/hl);
        % Set up random sources in box [-1/2, 1/2]^3
        N = 16;
        points =  rl*(rand(N, 3)-1/2);
        charges = rl*(rand(N, 1)-1/2);
        % Test sample point
        x = rl*0.1; y = rl*0.2; z = rl*0.3;
        % Setup Fourier vectors
        [m1, m2, m3] = ndgrid(-nf:nf, -nf:nf, -nf:nf);
        k1 = hl*m1(:);
        k2 = hl*m2(:);
        k3 = hl*m3(:);    
        D0hat = laplace_diffkernel_fourier(k1, k2, k3, sigma_l);
        % w_l factor 
        wl = 1/(2*pi)^3 * hl^3 * D0hat;
        % Source expansion
        kdoty = k1.*points(:, 1)' + k2.*points(:, 2)' + k3.*points(:, 3)';
        ghat = exp(-1i*kdoty) * charges;
        kdotx = k1.*x + k2.*y + k3.*z;
        % Compute udiff
        udiff_fourier = real( (wl .* ghat).' * exp(1i*kdotx) );
        % Compare to reference
        r = sqrt( (points(:, 1)-x).^2 + ...
                  (points(:, 2)-y).^2 + ...
                  (points(:, 3)-z).^2 );
        D0 = (erf(r/sigma_lp1) - erf(r/sigma_l))./r;
        % u_diff computed directly
        udiff_direct = charges'*D0;
        % Compare direct and Fourier
        testCase.verifyEqual(udiff_fourier, udiff_direct, 'abstol', tol);
        % Form proxy charges:
        proxy_charges = points2proxy(points*2/rl, charges, p);
        % Convert proxy charges to outgoing
        Phi = Tprox2pw(proxy_charges, l);
        % Convert to incoming (same box)
        Psi = wl.*Phi;
        % Compare to directly computed incoming
        testCase.verifyLessThan(norm(wl.*ghat-Psi, inf), tol);
        % Convert to local
        lambda = Tpw2poly(Psi, l);
        % Evaluate expansion at sample point
        Ex = approx.chebevalmat(x*2/rl, p);
        Ey = approx.chebevalmat(y*2/rl, p);
        Ez = approx.chebevalmat(z*2/rl, p);
        udiff_expa = real( kron(Ez, kron(Ey, Ex))* lambda );
        testCase.verifyEqual(udiff_expa, udiff_direct, 'abstol', tol);
    end
end
