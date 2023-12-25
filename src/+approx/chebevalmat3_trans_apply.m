function g = chebevalmat3_trans_apply(x, y, z, p, f, M)
% Fast evaluation of
% g = transpose(chebevalmat3(x, y, z, p)*M)*f
    Ex = approx.chebevalmat(x, p)*M;
    Ey = approx.chebevalmat(y, p)*M;
    Ez = approx.chebevalmat(z, p)*M;
    Exf = (Ex .* f).'; % Precompute Ex columns times f
    g = zeros(p, p, p);
    for i3=1:p
        M = Ey .* Ez(:, i3);
        g(:, :, i3) = Exf * M;
    end
    g = reshape(g, [], 1);
end
