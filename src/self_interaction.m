function u = self_interaction(tree, charges, sigma_0)
    u = zeros(tree.N, 1);
    for leaf_idx=tree.nonempty_leafs()
        idxs = tree.box_point_idxs(leaf_idx);
        l = tree.boxLevels(leaf_idx);
        sigma_l = sigma_0 / 2^l;
        u(idxs) = -charges(idxs)*2/(sqrt(pi)*sigma_l);
    end
end

