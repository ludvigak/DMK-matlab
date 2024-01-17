function u=local_interactions(tree, charges, kernel)
    u = zeros(tree.N, kernel.dim_out);
    leaf_list = tree.nonempty_leafs();
    Nleafs = numel(leaf_list);
    udata = cell(1, Nleafs);
    % Run in parallell if there is a parpool started
    p = gcp('nocreate');
    if isempty(p)
        % There is no parallel pool
        poolsize = 0;
    else
        % There is a parallel pool of <p.NumWorkers> workers
        poolsize = p.NumWorkers;
    end
    %parfor (i=1:Nleafs, poolsize)
    for i=1:Nleafs
        leaf_idx=leaf_list(i);
        l = tree.boxLevels(leaf_idx);
        idxs = tree.box_point_idxs(leaf_idx);
        box_points = tree.points(idxs, :);
        clist = tree.boxColleagues{leaf_idx};
        ubox = zeros(numel(idxs), kernel.dim_out);
        for coll=clist
            coll_idxs = tree.box_point_idxs(coll);
            coll_points = tree.points(coll_idxs, :);
            coll_charges = charges(coll_idxs, :);
            ubox = ubox + kernel.reskernel(box_points, coll_points, coll_charges, l);
        end
        % Store output data from this iteration
        udata{i} = ubox;
    end
    % Copy output data (from threads) to output structure
    for i=1:Nleafs
        leaf_idx=leaf_list(i);
        idxs = tree.box_point_idxs(leaf_idx);    
        u(idxs, :) = udata{i};
    end
end

