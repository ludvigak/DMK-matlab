%% Main function to generate tests
function tests = test_periodic()
tests = functiontests(localfunctions);
end

%% Test Functions
function run_periodic_dmk(testCase, kernel, dmk_tol, test_tol)
    rng(1);
    N = 100;
    max_level = 1;
    dmk_opt   = dmk_default_opts(tolerance=dmk_tol, kernel=kernel, periodic=true);
    points = rand(N, 3)-1/2;
    charges = rand(N, dmk_opt.kernel.dim_in)-1/2;
    charges = charges - sum(charges, 1)/N;
    charges(end, :) = charges(end, :)-sum(charges, 1);    
    dmk_state = dmk_init(points, max_level, dmk_opt);
    u_dmk     = dmk_apply(charges, dmk_state);
    % Stricter tolerance for Ewald sum: also avoids using exactly same kernel params
    u_ref = ewald_sum(points, charges, kernel(tolerance=dmk_tol/10));
    max_rel_err = norm(u_ref - u_dmk, inf) / norm(u_ref, inf);
    testCase.verifyLessThan(max_rel_err, test_tol);
end

function test_periodic_dmk_laplace_ewald(testCase)
    run_periodic_dmk(testCase, @kernels.laplace_ewald, 1e-10, 1e-10);
end

function test_periodic_dmk_laplace_pswf(testCase)
    run_periodic_dmk(testCase, @kernels.laplace_pswf, 1e-10, 1e-10);
end

function test_periodic_dmk_stokeslet_hasimoto(testCase)
    run_periodic_dmk(testCase, @kernels.stokeslet_hasimoto, 1e-10, 1e-10);
end

function test_periodic_dmk_stokeslet_pswf(testCase)
    run_periodic_dmk(testCase, @kernels.stokeslet_pswf, 1e-10, 1e-10);
end

function test_periodic_dmk_stresslet_hasimoto(testCase)
    run_periodic_dmk(testCase, @kernels.stresslet_hasimoto, 1e-10, 1e-10);
end

function test_periodic_dmk_stresslet_pswf(testCase)
    run_periodic_dmk(testCase, @kernels.stresslet_pswf, 1e-10, 1e-10);
end

function test_periodic_dmk_rotlet_hasimoto(testCase)
    run_periodic_dmk(testCase, @kernels.rotlet_ewald, 1e-10, 1e-10);
end

function test_periodic_dmk_rotlet_pswf(testCase)
    run_periodic_dmk(testCase, @kernels.rotlet_pswf, 1e-10, 1e-10);
end
