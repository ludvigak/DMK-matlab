function u = laplace_kernel(target, points, charges)
    N = numel(charges);
    assert(size(points, 1)==N)
    assert(size(points, 2)==3)
    assert(size(target, 1)==1)
    assert(size(target, 2)==3)
    r = sqrt(sum((target - points).^2, 2));
    mask = (r ~= 0);
    u = sum(charges(mask)./r(mask));
end

