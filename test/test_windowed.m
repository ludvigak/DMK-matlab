%% Main function to generate tests
function tests = test_windowed
tests = functiontests(localfunctions);
end

%% Test Functions


function test_windowed_kernel(testCase)
    p = 40;
    N = 20;
    tol = 1e-13;
    sigma_0 = 1/sqrt(log(1/tol));
    % Setup test
    points = rand(N, 3) - 1/2;
    charges = rand(N, 1) - 1/2;
    % Run upward pass and collect root proxy charges and points
    tree = octree(points, 1);
    proxy_charges = init_proxy_charges(tree, charges, p);
    [rvec, V] = approx.chebvander(p);
    box_proxy_points = tree.box_grid(1, rvec);
    box_proxy_charges = proxy_charges{1};
    % Fourier setup
    Kmax = ceil( 2*log(1/tol) );
    nf = Kmax;
    hf = Kmax/nf;
    Ctrunc = sqrt(3) + 6*sigma_0;    
    Twin = operator_windowed(p, hf, nf, Ctrunc, sigma_0);
    % Evaluate windowed kernel at proxy points
    far_expa = Twin(box_proxy_charges);
    field_fourier = approx.kronmat_apply(V, far_expa, 3);
    % Evaluate mollified potential directly from source points
    targets = box_proxy_points.';
    r = sqrt( (points(:, 1)-targets(1, :)).^2 + ...
              (points(:, 2)-targets(2, :)).^2 + ...
              (points(:, 3)-targets(3, :)).^2 );
    field_direct = (erf(r/sigma_0)./r).' * charges;
    % Compare
    Emax = norm(field_direct - field_fourier, inf);
    testCase.verifyLessThan(Emax, tol);
end
