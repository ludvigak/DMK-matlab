function E = chebevalmat2(x, y, p)
    Ex = approx.chebevalmat(x, p);
    Ey = approx.chebevalmat(y, p);
    [i1, i2] = ndgrid(1:p);
    E = Ex(:, i1(:)) .* Ey(:, i2(:));
end
