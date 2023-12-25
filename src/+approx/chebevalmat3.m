function E = chebevalmat3(x, y, z, p, varargin)
    if isempty(varargin)
        M = 1;
    else
        M = varargin{1};
    end
    Ex = approx.chebevalmat(x, p)*M;
    Ey = approx.chebevalmat(y, p)*M;
    Ez = approx.chebevalmat(z, p)*M;
    [i1, i2, i3] = ndgrid(1:p);
    E = Ex(:, i1(:)) .* Ey(:, i2(:)) .* Ez(:, i3(:));
end
