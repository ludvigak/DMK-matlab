function local_expansions = collect_local_expansions(tree, proxy_charges, ...
                                                     Troot, Tprox2pw, Tpw2poly, Tpwshift, Tp2c ...
                                                    )
% DMK downward pass
    outgoing_expansions = cell(1, tree.numBoxes);
    local_expansions    = cell(1, tree.numBoxes);
    % Handle root level outside loop (Troot = periodic or windowed)
    root = 1;
    local_expansions{root} = Troot(proxy_charges{root});
    % Downward pass
    for l=0:tree.maxLevel-1
        rl = 1/2^l;
        box_list = find(tree.boxLevels==l);
        Nlevel = numel(box_list); % Boxes on level
        % Form outgoing
        for idx=1:Nlevel
            box = box_list(idx);
            box_proxy_charges = proxy_charges{box};
            outgoing_expansions{box} = Tprox2pw(box_proxy_charges, l);
        end
        % Collect incoming (including self)
        for idx=1:Nlevel
            box = box_list(idx);
            clist = tree.boxColleagues{box};
            slist = tree.boxColleagueShifts{box};
            incoming = 0;
            for cidx=1:numel(clist)
                coll = clist(cidx);
                coll_shift = slist(cidx, :);
                coll_center = tree.box_center(coll);
                pw = outgoing_expansions{coll};
                shift = (coll_center+coll_shift-tree.box_center(box)) / rl;
                incoming = incoming + Tpwshift(pw, shift(1), shift(2), shift(3));
            end
            % Convert to local
            local = Tpw2poly(incoming, l);
            local_expansions{box} = local_expansions{box} + local;
            % Shift to children
            for k=1:2
                for j=1:2
                    for i=1:2
                        child = tree.boxChildren{box}(i,j,k);
                        expa_child = Tp2c(i, j, k, local_expansions{box});
                        local_expansions{child} = expa_child;
                    end
                end
            end
            % Free eexpansion
            local_expansions{box} = [];
        end
        % Free outgoing expansions
        for idx=1:Nlevel
            box = box_list(idx);
            outgoing_expansions{box} = [];
        end
    end % end levels
end
