function g = chebevalmat3_trans_apply(x, y, z, p, f, M)
% Fast evaluation of
% g = transpose(chebevalmat3(x, y, z, p)*M)*f
    N = numel(x);
    Ex = approx.chebevalmat(x, p)*M;
    Ey = approx.chebevalmat(y, p)*M;
    Ez = approx.chebevalmat(z, p)*M;
    % Precompute Ex columns times f
    Exf = Ex .* f; % (N, p)
    % Precompute outer product of Ey and Ez (memory waste...)
    Eyz = reshape(Ey, N, p, 1) .* reshape(Ez, N, 1, p); % (N, p, p)
    % Contract
    g = pagemtimes(Exf, 'transpose', Eyz, 'none');  % (p, p, p)
    g = reshape(g, [], 1);
end
