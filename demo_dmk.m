clear; rng(1);

tol = 1e-6;
N = 3e4;
max_level = 2;

points = rand(N, 3)-1/2;
charges = rand(N, 1)-1/2;

disp("* DMK")
atic = tic();
dmk_opt   = dmk_default_opts(tolerance=tol, verbose=true);
dmk_state = dmk_init(points, max_level, dmk_opt);
u_dmk     = dmk_apply(charges, dmk_state);
toc(atic)

if N <= 1e5
    disp('* Direct eval')
    atic = tic();
    u_ref = dmk_opt.kernel.direct(points, points, charges);
    toc(atic)
    max_rel_err = norm(u_ref - u_dmk, inf) / norm(u_ref, inf)
end
