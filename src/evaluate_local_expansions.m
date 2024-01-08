function u = evaluate_local_expansions(local_expansions, tree, p)
    leaf_list = tree.nonempty_leafs();
    % Get output dimension from a leaf expansions
    dim = size(local_expansions{leaf_list(1)}, 2);
    u = zeros(tree.N, dim);
    for leaf_idx=leaf_list
        idxs = tree.box_point_idxs(leaf_idx);
        box_points = tree.points(idxs, :);
        c = tree.box_center(leaf_idx);
        s = tree.box_size(leaf_idx);
        scaled_points = (box_points-c)*2./s;
        for d=1:dim
            expa = local_expansions{leaf_idx};
            u(idxs,d) = approx.chebevalmat3_apply(scaled_points(:, 1), ...
                                                  scaled_points(:, 2), ...
                                                  scaled_points(:, 3), ...
                                                  p, ...
                                                  expa(:,d));
        end
    end    
end

