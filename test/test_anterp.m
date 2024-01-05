%% Main function to generate tests
function tests = test_anterp
tests = functiontests(localfunctions);
end

%% Test Functions
function test_upward_pass(testCase)
    N = 64;
    p = 32;
    max_level = 1;
    points = rand(N, 3)-1/2;
    charges = rand(N, 1)-1/2;
    rvec = chebpts(p, 1);
    tree = octree(points, max_level);
    proxy_charges = init_proxy_charges(tree, charges, p);
    kernel = kernels.laplace_ewald();
    % For each box, check that proxy charges represent the field from
    % all points in box (including children), at a well-separated distance
    for box_idx=1:tree.numBoxes
        idxs = tree.box_point_idxs(box_idx);
        box_points = points(idxs, :);
        box_charges = charges(idxs);
        c = tree.box_center(box_idx);
        s = tree.box_size(box_idx);
        box_proxy_points = tree.box_grid(box_idx, rvec);
        box_proxy_charges = proxy_charges{box_idx};
        target = c + s; % half-box separation from center
        u = kernel.direct(target, box_points, box_charges);
        u_proxy = kernel.direct(target, box_proxy_points, box_proxy_charges);
        testCase.verifyEqual(u, u_proxy, 'reltol', eps(1000));
    end
end

