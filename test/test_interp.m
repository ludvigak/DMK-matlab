%% Main function to generate tests
function tests = test_interp
tests = functiontests(localfunctions);
end

%% Test Functions
function test_interp_p2c(testCase)
    % Test that anterp matrices actually interpolate from parent to child
    p = 25;
    f = @(r) exp(-(r(:,1) + r(:,2)/2 + r(:,3)/3).^2);
    [~, V] = approx.chebvander(p);
    rvec = chebpts(p, 1);
    Tp2c = operator_parent2child(p);    
    tree = octree([0,0,0], 1);
    box_idx = 1;
    box_proxy_points = tree.box_grid(box_idx, rvec);
    fb = f(box_proxy_points);
    expa = approx.kronmat_apply(inv(V), fb, 3);
    for k=1:2
        for j=1:2
            for i=1:2
                child_idx = tree.boxChildren{box_idx}(i,j,k);
                child_proxy_points = tree.box_grid(child_idx, rvec);
                fc = f(child_proxy_points);               
                child_expa = Tp2c(i, j, k, expa);
                fi = approx.kronmat_apply(V, child_expa, 3);
                testCase.verifyEqual(fi,fc, 'abstol', eps(100));
            end
        end
    end           
end
