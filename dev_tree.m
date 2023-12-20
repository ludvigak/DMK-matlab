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
[rvec, V] = approx.chebvander(p);
ViT = transpose(inv(V));

disp('* Build tree')
tic
tree = octree(points, max_level);
toc
tree

disp('* Upward pass')
tic
proxy_charges = init_proxy_charges(tree, charges, p);
toc


tol = 1e-8;
sigma_0 = 1/sqrt(log(1/tol));



uref = laplace_kernel(target, points, charges);

disp('=== Error estimates: ===')
r = @(x) abs(x-0.01);
D0 = @(x) (erf(r(x)/(sigma_0/2))-erf(r(x)/sigma_0))./r(x);
xeval = linspace(-0.5, 0.5)';
Di = approx.chebevalmat(2*xeval, p)*(V\D0(rvec/2));
fprintf("interp  = %.2e\n", norm(Di-D0(xeval), inf));
fprintf("trunc_R = %.2e\n", erfc(1/sigma_0));
fprintf("trunc_D = %.2e\n", (erf(1/(sigma_0/2))-erf(1/sigma_0)));
disp('================================')


u = 0;
% Single target eval
home_box = tree.pointBoxes(targ_idx)
box = home_box;

while box > 0
    box
    l = tree.boxLevels(box);
    sigma_l = sigma_0 / 2^l;
    sigma_lp1 = sigma_l / 2;
    clist = tree.boxColleagues{box};
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
            box_proxy_points = tree.box_grid(coll, rvec);
            box_proxy_charges = proxy_charges{coll};
            r = sqrt(sum((target-box_proxy_points).^2, 2));
            % K = D_{l}
            K = (erf(r/sigma_lp1) - erf(r/sigma_l))./r;
            if l==0
                % K += W_0
                K = K + erf(r/sigma_0)./r;
            end
            uK = sum(K .* box_proxy_charges);
            u = u + uK;
        end
    end
    if l==tree.maxLevel
        % Self interaction
        u = u - charges(targ_idx)*2/(sqrt(pi)*sigma_l);
    end
    
    % Next iteration: parent
    box = tree.boxParents(box);
end



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
