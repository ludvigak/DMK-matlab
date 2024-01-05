%% Main function to generate tests
function tests = test_stokeslet
tests = functiontests(localfunctions);
end

%% Test Functions


function test_kernel(testCase)
% Basic checks that the Stokeslet kernel does what we want
    kernel = kernels.stokeslet_hasimoto();
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
    tol = 1e-13;
    sigma_0 = 1/sqrt(log(1/tol));
    ures = kernel.reskernel(targets, sources, forces, sigma_0);
    % Should be properly decayed
    testCase.verifyLessThan(norm(ures(:), inf), tol);
    % Should be ~ same with very large sigma_0
    usame = kernel.reskernel(targets, sources, forces, 1e5);
    testCase.verifyEqual(uref, usame, reltol=1e-4)
end

function test_residual_decay(testCase)
    kernel = kernels.stokeslet_hasimoto();
    tol = 1e-13;
    sigma_0 = 1/sqrt(log(1/tol));    
    hasimoto_trunc_err = norm(kernel.diffkernel([1 0 0], [0 0 0], [1 1 1], sigma_0), inf);
    testCase.verifyLessThan(hasimoto_trunc_err, 1e-12);
end

function test_moll(testCase)
    kernel = kernels.stokeslet_hasimoto();
    Nsrc = 15;
    Ntrg = 5;
    rng(1);
    sources = rand(Nsrc, 3);
    forces  = rand(Nsrc, 3);
    targets = rand(Ntrg, 3);
    tol = 1e-8;
    sigma_0 = 1/sqrt(log(1/tol));    
    umoll = kernel.mollkernel(targets, sources, forces, sigma_0);
    ures = kernel.reskernel(targets, sources, forces, sigma_0);
    u = kernel.direct(targets, sources, forces);
    testCase.verifyEqual(u, umoll+ures, reltol=1e-15)
end
