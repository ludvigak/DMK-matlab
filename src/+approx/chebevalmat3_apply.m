function g = chebevalmat3_apply(x, y, z, p, f)
% Fast evaluation of
% g = chebevalmat3(x, y, z, p)*f
    Ex = approx.chebevalmat(x, p);
    Ey = approx.chebevalmat(y, p);
    Ez = approx.chebevalmat(z, p);
    F = reshape(f, p, p, p);
    % Contract Ex and F
    ExF = pagemtimes(Ex, F); % (N, p, p)
    % Contract with Ey
    ExyF = squeeze( sum(Ey .* ExF, 2) ); % (N, 1, p) squeezed to (N, p)
    % Contract with Ez
    g = sum(ExyF .* Ez, 2);
end
