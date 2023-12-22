function proxy_charges = init_proxy_charges(tree, charges, p)
% DMK upward pass
    Tc2p = operator_child2parent(p);
    proxy_charges = cell(1, tree.numBoxes);
    % Init with zeros
    for i=1:tree.numBoxes
        proxy_charges{i} = zeros(p^3, 1);
    end
    [~, V] = approx.chebvander(p);
    % Fill in at leaf boxes
    for leaf_idx=tree.nonempty_leafs()
        idxs = tree.box_point_idxs(leaf_idx);
        box_points = tree.points(idxs, :);
        box_charges = charges(idxs);
        c = tree.box_center(leaf_idx);
        s = tree.box_size(leaf_idx);
        scaled_points = (box_points-c)*2./s;
        proxy_charges{leaf_idx}(:) = points2proxy(scaled_points, box_charges, p);
    end  
    % Recursively anterpolate up the tree
    collect_proxies(1);
    function collect_proxies(boxNo)
        children = tree.boxChildren{boxNo};
        if isempty(children)
            % Leaf box, do nothing
            return
        end
        % Non-leaf box: let children collect first
        for i=1:numel(children)
            collect_proxies(children(i));
        end
        % Then anterpolate to my proxy points and add
        for k=1:2
            for j=1:2
                for i=1:2
                    child = children(i,j,k);
                    child_proxies = proxy_charges{child};
                    % Anterpolate
                    anterp_proxies = Tc2p(i, j, k, child_proxies);
                    proxy_charges{boxNo} = proxy_charges{boxNo} + anterp_proxies;
                end
            end
        end
    end
end

