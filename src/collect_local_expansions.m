function local_expansions = collect_local_expansions(tree, proxy_charges, ...
                                                     Twin, Tprox2pw, Tpw2poly, Tpwshift, Tp2c ...
                                                    )
% DMK downward pass
    outgoing_expansions = cell(1, tree.numBoxes);
    local_expansions    = cell(1, tree.numBoxes);
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
        % Collect incoming
        for idx=1:Nlevel
            box = box_list(idx);
            clist = tree.boxColleagues{box};
            incoming = 0;
            for coll=clist
                pw = outgoing_expansions{coll};
                shift = (tree.box_center(coll)-tree.box_center(box)) / rl;
                incoming = incoming + Tpwshift(pw, shift(1), shift(2), shift(3));
            end
            % Convert to local
            local = Tpw2poly(incoming, l);
            if box==1
                % At root level, also add expansion of windowed kernel
                local = local + Twin(proxy_charges{1});
            end
            if isempty(local_expansions{box})
                local_expansions{box} = local;
            else
                local_expansions{box} = local_expansions{box} + local;
            end
            % Shift to child
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
            outgoing_expansions{box} = [];
        end        
    end % end levels
end

