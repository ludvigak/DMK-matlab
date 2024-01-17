%% Main function to generate tests
function tests = test_dmk()
tests = functiontests(localfunctions);
end

%% Test Functions
function test_point2point_laplace(testCase)
% Test complete point to point sum
    rng(1);
    tol = 1e-9;
    N = 1e3;
    max_level = 1;
    points = rand(N, 3)-1/2;
    charges = rand(N, 1)-1/2;
    dmk_opt   = dmk_default_opts(tolerance=tol);
    dmk_state = dmk_init(points, max_level, dmk_opt);
    u_dmk     = dmk_apply(charges, dmk_state);
    u_ref = dmk_opt.kernel.direct(points, points, charges);
    max_rel_err = norm(u_ref - u_dmk, inf) / norm(u_ref, inf);
    testCase.verifyLessThan(max_rel_err, tol);
end

function test_point2point_laplace_pswf(testCase)
% Test complete point to point sum
    rng(1);
    tol = 1e-9;
    N = 1e3;
    max_level = 1;
    points = rand(N, 3)-1/2;
    charges = rand(N, 1)-1/2;
    dmk_opt   = dmk_default_opts(tolerance=tol, kernel=@kernels.laplace_pswf);
    dmk_state = dmk_init(points, max_level, dmk_opt);
    u_dmk     = dmk_apply(charges, dmk_state);
    u_ref = dmk_opt.kernel.direct(points, points, charges);
    max_rel_err = norm(u_ref - u_dmk, inf) / norm(u_ref, inf);
    testCase.verifyLessThan(max_rel_err, tol);
end

function test_point2point_stokeslet(testCase)
% Test complete point to point sum
    rng(1);
    tol = 1e-9;
    N = 1e3;
    max_level = 1;
    points = rand(N, 3)-1/2;
    charges = rand(N, 3)-1/2;
    dmk_opt   = dmk_default_opts(tolerance=tol, kernel=@kernels.stokeslet_hasimoto);
    dmk_state = dmk_init(points, max_level, dmk_opt);
    u_dmk     = dmk_apply(charges, dmk_state);
    u_ref = dmk_opt.kernel.direct(points, points, charges);
    max_rel_err = norm(u_ref(:) - u_dmk(:), inf) / norm(u_ref(:), inf);
    testCase.verifyLessThan(max_rel_err, tol);
end

