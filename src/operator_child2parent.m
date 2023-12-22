function Tc2p = operator_child2parent(p)
% Precompute matrices for anterpolation
    [rvec, V] = approx.chebvander(p);
    Vi = inv(V);
    UT = cell(2, 1);
    x1 = (rvec-1)/2; % [-1, 0]
    x2 = (rvec+1)/2; % [ 0, 1]
    UT{1} = transpose(approx.chebevalmat(x1, p)*Vi); % Why more accurate than /V ?
    UT{2} = transpose(approx.chebevalmat(x2, p)*Vi);
    function parent_proxies = apply(i, j, k, child_proxies)
    % Apply operator for child at index (i,j,k) in 2^3 cube
        UxT = UT{i};
        UyT = UT{j};
        UzT = UT{k};
        parent_proxies = approx.kronmat3_apply(UzT, UyT, UxT, child_proxies);
    end
    Tc2p = @apply;
end

