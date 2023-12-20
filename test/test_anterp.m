%% Main function to generate tests
function tests = test_anterp
tests = functiontests(localfunctions);
end

%% Test Functions
function test_c2p_1d(testCase)
    p = 30;
    Tc2p = precompute_child2parent(p);    
    rvec = chebpts(p, 1);
    x1 = (rvec-1)/2;
    x2 = (rvec+1)/2;
    f = @(x) exp(-(x+1).^2);
    finterp1 = Tc2p.UT{1}' * f(rvec);
    finterp2 = Tc2p.UT{2}' * f(rvec);
    testCase.verifyEqual(finterp1, f(x1), 'reltol', eps(100));
    testCase.verifyEqual(finterp2, f(x2), 'reltol', eps(100));    
end

function test_interp_c2p(testCase)
    % Test that anterp matrices actually interpolate from parent to child
    p = 25;
    rvec = chebpts(p, 1);
    Tc2p = precompute_child2parent(p);    
    tree = octree([0,0,0], 1);
    box_idx = 1;
    box_proxy_points = tree.box_grid(box_idx, rvec);
    for k=1:2
        for j=1:2
            for i=1:2
                child_idx = tree.boxChildren{box_idx}(i,j,k);
                child_proxy_points = tree.box_grid(child_idx, rvec);
                f = @(r) exp(-(r(:,1) + r(:,2)/2 + r(:,3)/3).^2);
                fb = f(box_proxy_points);
                fc = f(child_proxy_points);
                UxT = Tc2p.UT{i};
                UyT = Tc2p.UT{j};
                UzT = Tc2p.UT{k};
                fi = approx.kronmat3_apply(UzT', UyT', UxT', fb);                
                testCase.verifyEqual(fi,fc, 'abstol', eps(100));
            end
        end
    end           
end

function test_upward_pass(testCase)
    N = 64;
    p = 32;
    max_level = 1;
    points = rand(N, 3)-1/2;
    charges = rand(N, 1)-1/2;
    rvec = chebpts(p, 1);
    tree = octree(points, max_level);
    proxy_charges = init_proxy_charges(tree, charges, p);
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
        u = laplace_kernel(target, box_points, box_charges);
        u_proxy = laplace_kernel(target, box_proxy_points, box_proxy_charges);
        testCase.verifyEqual(u, u_proxy, 'reltol', eps(1000));
    end
end

