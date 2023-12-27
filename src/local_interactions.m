function u=local_interactions(tree, charges, sigma_0)
    u = zeros(tree.N, 1);
    for leaf_idx=tree.nonempty_leafs()
        l = tree.boxLevels(leaf_idx);
        sigma_l = sigma_0 / 2^l;
        idxs = tree.box_point_idxs(leaf_idx);
        box_points = tree.points(idxs, :);
        clist = tree.boxColleagues{leaf_idx};
        ubox = zeros(numel(idxs), 1);
        for coll=clist
            coll_idxs = tree.box_point_idxs(coll);
            coll_points = tree.points(coll_idxs, :);
            coll_charges = charges(coll_idxs);
            ubox = ubox + laplace_reskernel(box_points, coll_points, coll_charges, sigma_l);
        end
        u(idxs) = ubox;
    end
end

