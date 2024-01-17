%% Main function to generate tests
function tests = test_planewave
tests = functiontests(localfunctions);
end

%% Test Functions


function run_planeswaves(testCase, kernel_ref)
    rng(1);
    tol = 1e-9;
    kernel = kernel_ref(tolerance=tol);
    p = 45;
    max_level = 2;
    r0 = 1;
    D = 3*r0;
    h0 = 2*pi/D;
    K0 = 2*kernel.Kmax;
    nf = ceil(K0/h0);
    dim = kernel.dim_in;
    Tprox2pw = operator_proxy2planewave(p, h0, nf, max_level);
    Tpw2poly = operator_planewave2local(p, h0, nf, max_level, kernel);
    Tpwshift = operator_planewave_shift(h0, nf);
    for l = 0:max_level
        rl = 1/2^l;
        hl = h0 / rl;
        Kl = 2/rl * kernel.Kmax;
        nf = ceil(Kl/hl);
        % Set up random sources in box [-1/2, 1/2]^3
        N = 20;
        points =  rl*(rand(N, 3)-1/2);
        charges = rl*(rand(N, dim)-1/2);
        % Test sample point
        target = rl*[0.4 0.2 0.3];
        % Setup Fourier vectors
        [m1, m2, m3] = ndgrid(-nf:nf, -nf:nf, -nf:nf);
        k1 = hl*m1(:);
        k2 = hl*m2(:);
        k3 = hl*m3(:);    
        D0hat = kernel.diffkernel_fourier(k1, k2, k3, l);
        % w_l factor 
        wl = 1/(2*pi)^3 * hl^3;
        % Source expansion
        kdoty = k1.*points(:, 1)' + k2.*points(:, 2)' + k3.*points(:, 3)';
        ghat = exp(-1i*kdoty) * charges;
        kdotx = k1.*target(1) + k2.*target(2) + k3.*target(3);
        % Compute udiff
        udiff_fourier = real( (wl .* D0hat(ghat)).' * exp(1i*kdotx) ).';
        % Compare to direct reference
        udiff_direct = kernel.diffkernel(target, points, charges, l);
        % Compare direct and Fourier
        testCase.verifyEqual(udiff_fourier, udiff_direct, 'abstol', tol);
        % Form proxy charges:
        proxy_charges = zeros(p^3, dim);
        for d=1:dim
            proxy_charges(:, d) = points2proxy(points*2/rl, charges(:, d), p);
        end
        % Convert proxy charges to outgoing
        Phi = Tprox2pw(proxy_charges, l);
        % Convert to incoming (same box)
        Psi = Phi;
        % Compare to directly computed incoming
        testCase.verifyLessThan(norm(wl.*D0hat(ghat)-wl.*D0hat(Psi), inf), tol);
        % Convert to local
        lambda = Tpw2poly(Psi, l);
        % Evaluate expansion at sample point
        Ex = approx.chebevalmat(target(1)*2/rl, p);
        Ey = approx.chebevalmat(target(2)*2/rl, p);
        Ez = approx.chebevalmat(target(3)*2/rl, p);
        udiff_expa = real( kron(Ez, kron(Ey, Ex))* lambda );
        testCase.verifyEqual(udiff_expa, udiff_direct, 'abstol', tol);
        % Shift sample point around together with expansions
        for i=-1:1
            for j=-1:1
                for k=-1:1
                    shift = [i j k];
                    newCenter = rl*shift;                   
                    Psi_shifted = Tpwshift(Psi, i, j, k);
                    lambda_shifted = Tpw2poly(Psi_shifted, l);
                    udiff_shifted = real( kron(Ez, kron(Ey, Ex))* lambda_shifted );
                    udiff_shifted_direct = kernel.diffkernel(target - newCenter, points, charges, l);
                    testCase.verifyEqual(udiff_shifted, ...
                                         udiff_shifted_direct, ...
                                         'abstol', tol);
                end
            end
        end
    end
end

function test_planeswaves_laplace(testCase)
    run_planeswaves(testCase, @kernels.laplace_ewald);
end

function test_planeswaves_laplace_pswf(testCase)
    run_planeswaves(testCase, @kernels.laplace_pswf);
end

function test_planeswaves_stokeslet(testCase)
    run_planeswaves(testCase, @kernels.stokeslet_hasimoto);
end
