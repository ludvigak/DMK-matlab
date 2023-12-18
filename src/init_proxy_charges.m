function proxy_charges = init_proxy_charges(tree, charges, p)
    proxy_charges = cell(1, tree.numBoxes);
    % Init with zeros
    for i=1:tree.numBoxes
        proxy_charges{i} = zeros(p^3, 1);
    end
    [~, V] = approx.chebvander(p);
    ViT = transpose(inv(V));
    % Fill in at leaf boxes
    for leaf_idx=tree.nonempty_leafs()
        idxs = tree.box_point_idxs(leaf_idx);
        box_points = tree.points(idxs, :);
        box_charges = charges(idxs);
        c = tree.box_center(leaf_idx);
        s = tree.box_size(leaf_idx);
        scaled_points = (box_points-c)*2./s;
        % Evaluation matrix for nonuniform points
        proj = approx.chebevalmat3_trans_apply(...
            scaled_points(:, 1), ...
            scaled_points(:, 2), ...
            scaled_points(:, 3), ...
            p, box_charges);
        proxy_charges{leaf_idx}(:) = approx.kronmat_apply(ViT, proj, 3);
    end
end
