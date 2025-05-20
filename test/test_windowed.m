%% Main function to generate tests
function tests = test_windowed
tests = functiontests(localfunctions);
end

%% Test Functions


function run_windowed_kernel(testCase, kernel_ref, args)
    arguments
        testCase
        kernel_ref
        args.tol = 1e-12
    end
    rng(1);
    p = 40;
    N = 20;
    tol = args.tol;
    tol_factor = 0.1; % The tests here are stricter than what is necessary for total DMK accuracy (TODO: tighten)
    opt = dmk_default_opts(tolerance=tol*tol_factor, kernel=kernel_ref, p=p);
    kernel = opt.kernel;
    % Setup test
    dim = kernel.dim_in;
    points = rand(N, 3) - 1/2;
    charges = rand(N, dim) - 1/2;
    assert(all(abs(sum(charges, 1)) > 1e-3)); % Cannot have charge neutrality
    % Run upward pass and collect root proxy charges and points
    tree = octree(points, 1);
    proxy_charges = init_proxy_charges(tree, charges, p, kernel);
    [rvec, V] = approx.chebvander(p);
    box_proxy_points = tree.box_grid(1, rvec);
    box_proxy_charges = proxy_charges{1};
    % Fourier setup
    % TODO: Move into kernel
    Kmax = ceil( kernel.Kmax );
    nf = opt.nf_win;
    hf = opt.hf_win;
    Ctrunc = opt.Ctrunc;
    Twin = operator_windowed(p, hf, nf, Ctrunc, kernel);
    % Evaluate windowed kernel at proxy points
    far_expa = Twin(box_proxy_charges);
    field_fourier = zeros(p^3, kernel.dim_out);
    for d=1:kernel.dim_out
        field_fourier(:, d) = approx.kronmat_apply(V, far_expa(:, d), 3);
    end
    % Evaluate mollified potential at proxies directly from source points
    targets = box_proxy_points;
    field_direct = kernel.direct(targets, points, charges) - ...
        kernel.reskernel(targets, points, charges, 0);
    % Compare
    Emax = norm(field_direct - field_fourier, inf) / norm(field_direct, inf);
    testCase.verifyLessThan(Emax, tol);
    % Self interaction at one proxy point
    single_charge = box_proxy_charges;
    single_charge(2:end, :) = 0;
    expa = Twin(single_charge);
    field_fourier = zeros(p^3, kernel.dim_out);
    for d=1:kernel.dim_out
        field_fourier(:, d) = approx.kronmat_apply(V, expa(:, d), 3);
    end
    uself = kernel.self_interaction(single_charge, 0);
    testCase.verifyEqual(uself(1, :), -field_fourier(1, :), abstol=tol)
end

function test_windowed_laplace(testCase)
    run_windowed_kernel(testCase, @kernels.laplace_ewald);
end

function test_windowed_pswf(testCase)
    run_windowed_kernel(testCase, @kernels.laplace_pswf);
end

function test_windowed_stokeslet(testCase)
    run_windowed_kernel(testCase, @kernels.stokeslet_hasimoto);
end

function test_windowed_stokeslet_pswf(testCase)
    run_windowed_kernel(testCase, @kernels.stokeslet_pswf);
end

function test_windowed_stokeslet_pswf_num(testCase)
    run_windowed_kernel(testCase, @kernels.stokeslet_pswf_num);
end

function test_windowed_stokeslet_pswf3(testCase)
    run_windowed_kernel(testCase, @kernels.stokeslet_pswf3, tol=1e-11);
end

function test_windowed_stokeslet_pswf_sq(testCase)
    run_windowed_kernel(testCase, @kernels.stokeslet_pswf_sq, tol=1e-7);
end

function test_windowed_stresslet(testCase)
    run_windowed_kernel(testCase, @kernels.stresslet_hasimoto);
end

function test_windowed_stresslet_pswf(testCase)
    run_windowed_kernel(testCase, @kernels.stresslet_pswf);
end

function test_windowed_rotlet(testCase)
    run_windowed_kernel(testCase, @kernels.rotlet_ewald);
end

function test_windowed_rotlet_pswf(testCase)
    run_windowed_kernel(testCase, @kernels.rotlet_pswf);
end
