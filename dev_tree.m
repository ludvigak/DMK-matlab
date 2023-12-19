clear all
rng(1);
N = 10000;
max_level = 3;

points = 2*rand(N, 3)-1;
charges = rand(N, 1);

p = 24;
[rvec, V] = approx.chebvander(p);
ViT = transpose(inv(V));

tree = octree(points, max_level)
Tc2p = precompute_child2parent(p);

tic
proxy_charges = init_proxy_charges(tree, charges, p);
toc

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
