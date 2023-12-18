function E = chebevalmat3(x, y, z, p)
    Ex = approx.chebevalmat(x, p);
    Ey = approx.chebevalmat(y, p);
    Ez = approx.chebevalmat(z, p);
    [i1, i2, i3] = ndgrid(1:p);
    E = Ex(:, i1(:)) .* Ey(:, i2(:)) .* Ez(:, i3(:));
end
