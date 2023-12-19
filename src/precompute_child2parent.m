function Tc2p = precompute_child2parent(p)
% Precompute matrices for anterpolation
    [rvec, V] = approx.chebvander(p);
    Vi = inv(V);
    UT = cell(2, 1);
    x1 = (rvec-1)/2; % [-1, 0]
    x2 = (rvec+1)/2; % [ 0, 1]
    UT{1} = transpose(approx.chebevalmat(x1, p)*Vi); % Why more accurate than /V ?
    UT{2} = transpose(approx.chebevalmat(x2, p)*Vi);
    Tc2p = struct();
    Tc2p.UT = UT;
end

