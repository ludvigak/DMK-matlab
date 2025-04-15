clear; rng(1);

%kernel = @kernels.laplace_ewald;
%kernel = @kernels.laplace_pswf;
%kernel = @kernels.stokeslet_hasimoto;
%kernel = @kernels.stokeslet_pswf;
%kernel = @kernels.stresslet_hasimoto;
%kernel = @kernels.stresslet_pswf;
%kernel = @kernels.rotlet_ewald;
%kernel = @kernels.rotlet_pswf;

tol       = 1e-10;
N         = 2000;
max_level = 1;
periodic  = false;

dmk_opt = dmk_default_opts(tolerance=tol, verbose=true, kernel=kernel, periodic=periodic);

points = rand(N, 3)-1/2;
charges = rand(N, dmk_opt.kernel.dim_in)-1/2;
if periodic
    charges = charges - sum(charges, 1)/N;
    charges(end, :) = charges(end, :)-sum(charges, 1);
    assert(all(abs(sum(charges)) < 1e-14))
end

disp("* DMK")
atic = tic();
dmk_state = dmk_init(points, max_level, dmk_opt);
u_dmk     = dmk_apply(charges, dmk_state);
toc(atic)

if N <= 1e5
    atic = tic();
    if periodic
        disp('* Direct Ewald eval')
        u_ref = ewald_sum(points, charges, dmk_opt.kernel);
    else
        disp('* Direct eval')
        u_ref = dmk_opt.kernel.direct(points, points, charges);
    end
    toc(atic)
    max_rel_err = norm(u_ref(:) - u_dmk(:), inf) / norm(u_ref(:), inf)
    l2_rel_err = norm(u_ref(:) - u_dmk(:), 2) / norm(u_ref(:), 2)
end
