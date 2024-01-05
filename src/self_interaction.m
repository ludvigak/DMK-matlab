function u = self_interaction(tree, charges, sigma_0, kernel)
    u = zeros(tree.N, 1);
    for leaf_idx=tree.nonempty_leafs()
        idxs = tree.box_point_idxs(leaf_idx);
        l = tree.boxLevels(leaf_idx);
        sigma_l = sigma_0 / 2^l;
        u(idxs) = kernel.self_interaction(charges(idxs), sigma_l);
    end
end

