function ures = laplace_reskernel(targets, points, charges, sigma_l)
% Laplace residual kernel R_l(r)
    if size(targets, 1) ~= 3
        targets = targets.';
    end
    assert(size(targets, 1)==3);
    r = sqrt( (points(:, 1)-targets(1, :)).^2 + ...
              (points(:, 2)-targets(2, :)).^2 + ...
              (points(:, 3)-targets(3, :)).^2 );
    Rl = erfc(r/sigma_l)./r;
    Rl(r==0) = 0;
    ures = (charges.' * Rl).';
end

