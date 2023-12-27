clear all
rng(1);
N = 10000;
max_level = 3;

points = rand(N, 3)-1/2;
charges = rand(N, 1)-1/2;

targ_idx = 1;
target = points(targ_idx,:)
%charges(targ_idx) = 0;

p = 40;
tol = 1e-8;
sigma_0 = 1/sqrt(log(1/tol));

% Setup root level windowed kernel
Kmax_win = ceil( 2*log(1/tol) );
nf_win = Kmax_win;
hf_win = Kmax_win/nf_win;
Ctrunc = sqrt(3) + 6*sigma_0;    
Tfar = operator_far(p, hf_win, nf_win, Ctrunc, sigma_0);
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

disp('* Upward pass')
tic
proxy_charges = init_proxy_charges(tree, charges, p);
toc


disp('* Downward pass')
tic
local_expansions = collect_local_expansions(tree, proxy_charges, ...
                                            Tfar, Tprox2pw, Tpw2poly, Tpwshift, Tp2c);
toc

uref = laplace_kernel(target, points, charges);



u = 0;
ufar = 0;
udiff = 0;
% Single target eval
home_box = tree.pointBoxes(targ_idx)
box = home_box;

while box > 0
    box
    l = tree.boxLevels(box);
    rl = 1/2^l;
    sigma_l = sigma_0 / 2^l;
    sigma_lp1 = sigma_l / 2;
    clist = tree.boxColleagues{box};
    %incoming = zeros(Nf^3, 1);
    for coll=clist
        if l==tree.maxLevel
            idxs = tree.box_point_idxs(coll);
            box_points = points(idxs, :);
            box_charges = charges(idxs);
            r = sqrt(sum((target-box_points).^2, 2));
            Rl = erfc(r/sigma_l)./r;
            Rl(r==0) = 0;
            uR = sum(Rl .* box_charges);
            u = u + uR;
        else
            if l==0
                % W_0
                box_proxy_points = tree.box_grid(coll, rvec);
                box_proxy_charges = proxy_charges{coll};
                r = sqrt(sum((target-box_proxy_points).^2, 2));
                ufar = ufar + sum(erf(r/sigma_0)./r .* box_proxy_charges);
            end
            % % Add difference kernel
            % %udiff = laplace_diffkernel(target, box_proxy_points, box_proxy_charges, sigma_l);

            % % Compute box expansion
            % pw = Tprox2pw(box_proxy_charges, l);
            % % Shift it
            % shift = (tree.box_center(coll)-tree.box_center(box)) / rl;
            % incoming = incoming + Tpwshift(pw, shift(1), shift(2), shift(3));
        end
    end
    %local_expa = Tpw2poly(incoming, l);
    % if l < tree.maxLevel
    %     local_expa = Tpw2poly(incoming_expansions{box}, l);
    %     scaled_target = (target-(tree.box_center(box)))*2/rl;
    %     Ex = approx.chebevalmat(scaled_target(1), p);
    %     Ey = approx.chebevalmat(scaled_target(2), p);
    %     Ez = approx.chebevalmat(scaled_target(3), p);
    %     udiff_expa = real( kron(Ez, kron(Ey, Ex))* local_expa );
    %     udiff = udiff+udiff_expa;
    % end

    
    if l==tree.maxLevel
        % Self interaction
        uself = - charges(targ_idx)*2/(sqrt(pi)*sigma_l);
    end
    
    % Next iteration: parent
    box = tree.boxParents(box);
end

scaled_target = (target-(tree.box_center(home_box)))*2*2^tree.maxLevel;
Ex = approx.chebevalmat(scaled_target(1), p);
Ey = approx.chebevalmat(scaled_target(2), p);
Ez = approx.chebevalmat(scaled_target(3), p);
udiff = real( kron(Ez, kron(Ey, Ex))* local_expansions{home_box} );


u = u+udiff+uself;



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
