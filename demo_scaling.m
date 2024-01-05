% Show linear scaling
clear; rng(1);
clf

tol = 1e-6;
ns  = 280;
N_list = [];
t_list = [];
for N=ceil(logspace(2, 5, 30))
    % Set max level to satisfy ns points per box
    max_level = max(ceil(log(N/ns)/log(8)), 1)
    points = rand(N, 3)-1/2;
    charges = rand(N, 1)-1/2;
    disp("* DMK")
    atic = tic();
    dmk_opt   = dmk_default_opts(tolerance=tol, verbose=true);
    dmk_state = dmk_init(points, max_level, dmk_opt);
    u_dmk     = dmk_apply(charges, dmk_state);
    toc(atic)
    wtime = toc(atic);
    N_list(end+1) = N;
    t_list(end+1) = wtime;
end
loglog(N_list, t_list, '.-', displayname='DMK')
hold on
legend
loglog(N_list, N_list * t_list(end)/N_list(end), '--', displayname='O(N)')
xlabel('N')
ylabel('time [s]')
