function udiff = laplace_diffkernel(targets, points, charges, sigma_l)
% Laplace difference kernel D_l(r)
    if size(targets, 1) ~= 3
        targets = targets.';
    end
    assert(size(targets, 1)==3);
    sigma_lp1 = sigma_l/2; % sigma_{l+1} = sigma_l / 2
    r = sqrt( (points(:, 1)-targets(1, :)).^2 + ...
              (points(:, 2)-targets(2, :)).^2 + ...
              (points(:, 3)-targets(3, :)).^2 );
    Dl = (erf(r/sigma_lp1) - erf(r/sigma_l))./r;
    udiff = Dl.' * charges;
end

