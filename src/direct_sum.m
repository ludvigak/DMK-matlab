function u=direct_sum(points, charges)
    N = numel(charges);
    assert(size(points, 1)==N)
    assert(size(points, 2)==3)
    u = zeros(N, 1);
    for i=1:N
        target = points(i, :);
        u(i) = laplace_kernel(target, points, charges);
    end
end

