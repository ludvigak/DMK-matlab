clear all
rng(1);
N = 3e4;
max_level = 2;

points = rand(N, 3)-1/2;
charges = rand(N, 1)-1/2;


tol = 1e-8;
p = 43;

%tol = 1e-3;
%p = 16;


sigma_0 = 1/sqrt(log(1/tol));

% Setup root level windowed kernel
Kmax_win = ceil( 2*log(1/tol) );
nf_win = Kmax_win;
hf_win = Kmax_win/nf_win;
Ctrunc = sqrt(3) + 6*sigma_0;    
Twin = operator_windowed(p, hf_win, nf_win, Ctrunc, sigma_0);
% Setup planewave ops
r0 = 1;
D = 3*r0;
h0 = 2*pi/D;
K0 = 4/1 * log(1/tol);
nf = ceil(K0/h0);
Nf = 2*nf+1;
Tprox2pw = operator_proxy2planewave(p, h0, nf, max_level);
Tpw2poly = operator_planewave2local(p, h0, nf, max_level, sigma_0);
Tpwshift = operator_planewave_shift(h0, nf);
Tp2c = operator_parent2child(p);


disp('=== Error estimates: ===')
[rvec, V] = approx.chebvander(p);
ViT = transpose(inv(V));
r = @(x) abs(x-0.01);
D0 = @(x) (erf(r(x)/(sigma_0/2))-erf(r(x)/sigma_0))./r(x);
xeval = linspace(-0.5, 0.5)';
Di = approx.chebevalmat(2*xeval, p)*(V\D0(rvec/2));
fprintf("interp  = %.2e\n", norm(Di-D0(xeval), inf));
fprintf("trunc_R = %.2e\n", erfc(1/sigma_0));
fprintf("trunc_D = %.2e\n", (erf(1/(sigma_0/2))-erf(1/sigma_0)));
disp('================================')



disp('* Build tree')
tic
tree = octree(points, max_level);
toc
tree

tic_dmk = tic();
tic_far = tic();
disp('* Upward pass')
tic
proxy_charges = init_proxy_charges(tree, charges, p);
toc


disp('* Downward pass')
tic
local_expansions = collect_local_expansions(tree, proxy_charges, ...
                                            Twin, Tprox2pw, Tpw2poly, Tpwshift, Tp2c);
toc

disp('* Eval local expansions')
tic
ufar = evaluate_local_expansions(local_expansions, tree, p);
toc
tfar = toc(tic_far);

disp('* Local interactions')
tic
ures = local_interactions(tree, charges, sigma_0);
toc
tlocal = toc();

uself = self_interaction(tree, charges, sigma_0);
u = ufar + ures + uself;
tdmk = toc(tic_dmk);

if N <= 1e5
    disp('* Direct eval')
    tic
    uref = laplace_kernel(points, points, charges);
    toc
    tdirect = toc();
else
    disp(' * No direct eval, problem too large')
    uref = 0;
    tdirect = 0;
end

disp(' ')
fprintf("| Time DMK: %.2f (far) + %.2f (local) =\t%.2f\n", tfar, tlocal, tdmk)
fprintf("| Time direct: \t\t\t\t%.2f\n", tdirect)
err = norm(u - uref, inf)

return

u = 0;
disp('* Eval')
tic
for box_idx=1:tree.numBoxes
    l = tree.boxLevels(box_idx);
    rl = 1/2^l;
    sigma_l = sigma_0 / 2^l;
    sigma_lp1 = sigma_l / 2;
    if l==tree.maxLevel
        idxs = tree.box_point_idxs(box_idx);
        box_points = points(idxs, :);
        box_charges = charges(idxs);
        d = (target-box_points);
        r = sqrt(sum(d.^2, 2));
        Rl = erfc(r/sigma_l)./r;
        Rl(max(abs(d),[],2) > rl | r==0) = 0;
        uR = sum(Rl .* box_charges);
        u = u + uR;
        continue
    end    
    box_proxy_points = tree.box_grid(box_idx, rvec);
    box_proxy_charges = proxy_charges{box_idx};
    d = target-box_proxy_points;
    r = sqrt(sum(d.^2, 2));
    % K = D_{l}
    K = (erf(r/sigma_lp1) - erf(r/sigma_l))./r;
    if l==0
        % K += W_0
        K = K + erf(r/sigma_0)./r;
    end
    K(max(abs(d),[],2) > rl) = 0;
    uK = sum(K .* box_proxy_charges);
    u = u + uK;
end
toc
u
uref
err = abs(u - uref)

return
u = 0;
disp('* Eval')
tic
for box_idx=1:tree.numBoxes
    l = tree.boxLevels(box_idx);
    rl = 1/2^l;
    sigma_l = sigma_0 / 2^l;
    sigma_lp1 = sigma_l / 2;
    if l==tree.maxLevel
        idxs = tree.box_point_idxs(box_idx);
        box_points = points(idxs, :);
        box_charges = charges(idxs);
        d = (target-box_points);
        r = sqrt(sum(d.^2, 2));
        Rl = erfc(r/sigma_l)./r;
        Rl(max(abs(d),[],2) > rl | r==0) = 0;
        uR = sum(Rl .* box_charges);
        u = u + uR;
        continue
    end    
    box_proxy_points = tree.box_grid(box_idx, rvec);
    box_proxy_charges = proxy_charges{box_idx};
    d = target-box_proxy_points;
    r = sqrt(sum(d.^2, 2));
    % K = D_{l}
    K = (erf(r/sigma_lp1) - erf(r/sigma_l))./r;
    if l==0
        % K += W_0
        K = K + erf(r/sigma_0)./r;
    end
    K(max(abs(d),[],2) > rl) = 0;
    uK = sum(K .* box_proxy_charges);
    u = u + uK;
end
toc
u
uref
err = abs(u - uref)



return

err_proxy = zeros(1, tree.numBoxes);
for box_idx=1:tree.numBoxes
    idxs = tree.box_point_idxs(box_idx);
    box_points = points(idxs, :);
    box_charges = charges(idxs);
    c = tree.box_center(box_idx);
    s = tree.box_size(box_idx);
    box_proxy_points = tree.box_grid(box_idx, rvec);
    box_proxy_charges = proxy_charges{box_idx};
    target = c + 3*s;
    u = laplace_kernel(target, box_points, box_charges);
    u_proxy = laplace_kernel(target, box_proxy_points, box_proxy_charges);
    err_proxy(box_idx) = abs(u-u_proxy);
end
err_proxy
norm(err_proxy, inf)
clf
plot(log10(err_proxy), '.')
