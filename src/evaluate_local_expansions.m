function u = evaluate_local_expansions(local_expansions, tree, p)
    u = zeros(tree.N, 1);
    for leaf_idx=tree.nonempty_leafs()
        idxs = tree.box_point_idxs(leaf_idx);
        box_points = tree.points(idxs, :);
        c = tree.box_center(leaf_idx);
        s = tree.box_size(leaf_idx);
        scaled_points = (box_points-c)*2./s;
        u(idxs) = approx.chebevalmat3_apply(scaled_points(:, 1), ...
                                            scaled_points(:, 2), ...
                                            scaled_points(:, 3), ...
                                            p, ...
                                            local_expansions{leaf_idx});
    end    
end

