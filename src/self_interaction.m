function u = self_interaction(tree, charges, kernel)
    u = zeros(tree.N, kernel.dim_out);
    for leaf_idx=tree.nonempty_leafs()
        idxs = tree.box_point_idxs(leaf_idx);
        level = tree.boxLevels(leaf_idx);
        u(idxs, :) = kernel.self_interaction(charges(idxs, :), level);
    end
end

