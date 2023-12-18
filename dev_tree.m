clear all

N = 10000;
max_level = 3;

points = 2*rand(N, 3)-1;
charges = rand(N, 1);

p = 24;
[rvec, V] = approx.chebvander(p);
ViT = transpose(inv(V));
[xp, yp, zp] = ndgrid(rvec, rvec, rvec);
proxy_points = [xp(:), yp(:), zp(:)];

tree = octree(points, max_level)
Tc2p = precompute_child2parent(p);

tic
proxy_charges = init_proxy_charges(tree, charges, p);
toc

err_proxy = 0;
for leaf_idx=tree.nonempty_leafs()
    idxs = tree.box_point_idxs(leaf_idx);
    box_points = points(idxs, :);
    box_charges = charges(idxs);
    
    c = tree.box_center(leaf_idx);
    s = tree.box_size(leaf_idx);
    box_proxy_points = proxy_points.*s/2 + c;
    box_proxy_charges = proxy_charges{leaf_idx};
    target = c + s;
    u = laplace_kernel(target, box_points, box_charges);
    u_proxy = laplace_kernel(target, box_proxy_points, box_proxy_charges);
    err_proxy = max(err_proxy, abs(u-u_proxy));
end
err_proxy
