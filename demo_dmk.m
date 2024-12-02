clear; rng(1);

%kernel = @kernels.laplace_ewald;
%kernel = @kernels.laplace_pswf;
%kernel = @kernels.stokeslet_hasimoto;
%kernel = @kernels.stokeslet_pswf;
%kernel = @kernels.stokeslet_pswf2;
%kernel = @kernels.stresslet_hasimoto;
kernel = @kernels.stresslet_pswf;

tol = 1e-12;
N = 2000;
max_level = 1;

dmk_opt   = dmk_default_opts(tolerance=tol, verbose=true, kernel=kernel);

points = rand(N, 3)-1/2;
charges = rand(N, dmk_opt.kernel.dim_in)-1/2;

disp("* DMK")
atic = tic();
dmk_state = dmk_init(points, max_level, dmk_opt);
u_dmk     = dmk_apply(charges, dmk_state);
toc(atic)

if N <= 1e5
    disp('* Direct eval')
    atic = tic();
    u_ref = dmk_opt.kernel.direct(points, points, charges);
    toc(atic)
    max_rel_err = norm(u_ref(:) - u_dmk(:), inf) / norm(u_ref(:), inf)
end



