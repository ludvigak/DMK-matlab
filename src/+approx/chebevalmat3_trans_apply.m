function g = chebevalmat3_trans_apply(x, y, z, p, f)
% Fast evaluation of
% g = transpose(chebevalmat3(x, y, z, p))*f
    Ex = approx.chebevalmat(x, p);
    Ey = approx.chebevalmat(y, p);
    Ez = approx.chebevalmat(z, p);
    g = zeros(p^3, 1);
    Exf = (Ex .* f)'; % Precompute Ex columns times f
    i1=1:p;
    for i3=1:p
        for i2=1:p
            k = i1 + p*(i2-1 + p*(i3-1));
            g(k) = Exf * (Ey(:, i2) .* Ez(:, i3));
        end
    end
end
