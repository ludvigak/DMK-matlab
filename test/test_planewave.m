%% Main function to generate tests
function tests = test_planewave
tests = functiontests(localfunctions);
end

%% Test Functions


function test_outgoing(testCase)
    tol = 1e-9;
    p = 45;
    l = 0;
    rl = 1/2^l;
    sigma_0 = 1/sqrt(log(1/tol));
    sigma_1 = sigma_0 / 2;
    D = 3;
    hl = 2*pi/D;
    Kl = 4/rl * log(1/tol);
    nf = ceil(Kl/hl);
    % Set up random sources in box [-1/2, 1/2]^3
    N = 100;
    points = rand(N, 3)-1/2;
    charges = rand(N, 1)-1/2;
    % Random sample point
    x = 0.1; y = 0.2; z = 0.3;
    % Setup Fourier vectors
    [m1, m2, m3] = ndgrid(-nf:nf, -nf:nf, -nf:nf);
    k1 = hl*m1(:);
    k2 = hl*m2(:);
    k3 = hl*m3(:);
    ksq = k1.^2 + k2.^2 + k3.^2;
    % Difference kernel Fourier transform
    D0hat = 4*pi*(exp(-ksq*sigma_1^2/4) - exp(-ksq*sigma_0^2/4))./ksq;
    D0hat(k1==0 & k2==0 & k3==0) = pi*(sigma_0^2-sigma_1^2);
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
    D0 = (erf(r/sigma_1) - erf(r/sigma_0))./r;
    % u_diff computed directly
    udiff_direct = charges'*D0;
    % Compare direct and Fourier
    testCase.verifyEqual(udiff_fourier, udiff_direct, 'abstol', tol);
    % Form proxy charges:
    [rvec, V] = approx.chebvander(p);
    ViT = transpose(inv(V));
    scaled_points = points*2;
    proj = approx.chebevalmat3_trans_apply(...
        scaled_points(:, 1), ...
        scaled_points(:, 2), ...
        scaled_points(:, 3), ...
        p, charges);
    proxy_charges = approx.kronmat_apply(ViT, proj, 3);
    % Convert proxy charges to outgoing
    k = hl*(-nf:nf)';
    Mout = exp(-1i*k.*rvec'/2); % note scaling
    Phi = approx.kronmat3_apply(Mout, Mout, Mout, proxy_charges);
    % Convert to incoming (same box)
    Psi = wl.*Phi;
    % Compare to directly computed incoming
    testCase.verifyLessThan(norm(wl.*ghat-Psi, inf), tol);
    % Evaluate at proxy points
    Minc = Mout';
    udiff_proxy = approx.kronmat3_apply(Minc, Minc, Minc, Psi);
    % Expand at proxy points
    lambda = approx.kronmat_apply(inv(V), udiff_proxy, 3);
    % Evaluate expansion at sample point
    Ex = approx.chebevalmat(2*x, p);
    Ey = approx.chebevalmat(2*y, p);
    Ez = approx.chebevalmat(2*z, p);
    udiff_expa = real( kron(Ez, kron(Ey, Ex))* lambda );
    testCase.verifyEqual(udiff_expa, udiff_direct, 'abstol', tol);
end
