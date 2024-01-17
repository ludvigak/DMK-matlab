%% Main function to generate tests
function tests = test_stokeslet
tests = functiontests(localfunctions);
end

%% Test Functions


function test_kernel(testCase)
% Basic checks that the Stokeslet kernel does what we want
    tol = 1e-13;
    kernel = kernels.stokeslet_hasimoto(tolerance=tol);
    testCase.verifyEqual(kernel.dim_in, 3);
    testCase.verifyEqual(kernel.dim_out, 3);

    Nsrc = 5;
    Ntrg = 2;
    rng(1);
    sources = rand(Nsrc, 3);
    forces  = rand(Nsrc, 3);
    targets = rand(Ntrg, 3)+[2 0 0]; % Well-separated targets

    uref = zeros(Ntrg, 3);
    for idx_src=1:Nsrc
        for idx_trg=1:Ntrg
            rvec = targets(idx_trg, :) - sources(idx_src, :);
            r = sqrt(sum(rvec.^2));
            f = forces(idx_src, :);
            for i=1:3
                for j=1:3
                    uref(idx_trg, i) = uref(idx_trg, i) + ...
                        (i==j)*f(j)/r + rvec(i)*rvec(j)*f(j)/r^3;
                end
            end
        end
    end
    u = kernel.direct(targets, sources, forces);
    testCase.verifyEqual(u, uref, abstol=1e-14);
    % Test residual
    ures = kernel.reskernel(targets, sources, forces, 0);
    % Should be properly decayed
    testCase.verifyLessThan(norm(ures(:), inf), tol);
    % Should be ~ same with very large sigma
    usame = kernel.reskernel(targets, sources, forces, -20); % Artifically make sigma large
    testCase.verifyEqual(uref, usame, reltol=1e-4)
end

function test_residual_decay(testCase)
    tol = 1e-13;
    kernel = kernels.stokeslet_hasimoto(tolerance=tol);
    level = 0;
    hasimoto_trunc_err = norm(kernel.diffkernel([1 0 0], [0 0 0], [1 1 1], level), inf);
    testCase.verifyLessThan(hasimoto_trunc_err, 1e-12);
end

function test_moll(testCase)
    tol = 1e-8;
    kernel = kernels.stokeslet_hasimoto(tolerance=tol);
    Nsrc = 15;
    Ntrg = 5;
    rng(1);
    sources = rand(Nsrc, 3);
    forces  = rand(Nsrc, 3);
    targets = rand(Ntrg, 3);
    umoll = kernel.mollkernel(targets, sources, forces, 0);
    ures = kernel.reskernel(targets, sources, forces, 0);
    u = kernel.direct(targets, sources, forces);
    testCase.verifyEqual(u, umoll+ures, reltol=1e-15)
end

function test_diff(testCase)
% Simple check that diffkernel is defined correctly in relation to mollkernel
    tol = 1e-8;
    kernel = kernels.stokeslet_hasimoto(tolerance=tol);
    Nsrc = 15;
    Ntrg = 5;
    rng(1);
    sources = rand(Nsrc, 3);
    forces  = rand(Nsrc, 3);
    targets = rand(Ntrg, 3);
    umoll0 = kernel.mollkernel(targets, sources, forces, 0);
    umoll1 = kernel.mollkernel(targets, sources, forces, 1);
    udiff0 = kernel.diffkernel(targets, sources, forces, 0);
    testCase.verifyEqual(udiff0, umoll1-umoll0, reltol=1e-15)
    % Fourier modes
    [k1, k2, k3] = ndgrid(-5:5);
    k1 = k1(:); k2 = k2(:); k3 = k3(:);
    Nf = numel(k1);
    Fmoll0 = kernel.mollkernel_fourier(k1, k2, k3, 0);
    Fmoll1 = kernel.mollkernel_fourier(k1, k2, k3, 1);
    Fdiff0 = kernel.diffkernel_fourier(k1, k2, k3, 0);
    fhat = rand(Nf, 3) + 1i*rand(Nf, 3);
    testCase.verifyEqual(Fdiff0(fhat), (Fmoll1(fhat)-Fmoll0(fhat)), abstol=1e-14);
end
