function u = laplace_kernel(targets, points, charges)
    N = numel(charges);
    assert(size(points, 1)==N)
    assert(size(points, 2)==3)
    if size(targets, 1) ~= 3
        targets = targets.';
    end
    assert(size(targets, 1)==3);
    u = zeros(size(targets, 2), 1);
    batch_size = 64;
    % Vectorize over batch size
    for idx_begin=1:batch_size:N
        idx_end = min(N, idx_begin+batch_size-1);
        range =  idx_begin:idx_end;
        r = sqrt( (points(range, 1)-targets(1, :)).^2 + ...
                  (points(range, 2)-targets(2, :)).^2 + ...
                  (points(range, 3)-targets(3, :)).^2 );
        K = 1./r;
        K(r==0) = 0;
        u = u + (charges(range).' * K).';
    end
end

