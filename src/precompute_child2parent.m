function Tc2p = precompute_child2parent(p)
% Precompute matrices for anterpolation
    [rvec, V] = approx.chebvander(p)
    ViT = transpose(inv(V));
    E = cell(2, 1);
    x1 = (rvec-1)/2
    x2 = (rvec+1)/2
    ET{1} = transpose(approx.chebevalmat(x1, p));
    ET{2} = transpose(approx.chebevalmat(x2, p));
    Tc2p = struct();
    Tc2p.ViT = ViT;
    Tc2p.ET = ET;
end

