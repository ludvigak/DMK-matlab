classdef octree < handle
    % Simple linear octree
    
    properties (Access=public)
        N
        maxLevel
        periodic
        points
        boxLevels
        boxParents
        boxChildren
        boxColleagues
        boxColleagueShifts
        boxCorners
        boxPoints
        pointBoxes
        numBoxes
        verbose
    end

    methods
        function tree = octree(points, maxLevel, options)
            arguments
                points (:,3)
                maxLevel
                options.periodic (1,1) logical = false
            end
        % Points are Nx3 in [-1/2, 1/2]^3
            tree.N = size(points, 1);
            tree.maxLevel = maxLevel;
            tree.periodic = options.periodic;
            tree.points = points;
            tree.boxLevels = 0;
            tree.boxParents = 0;
            tree.boxChildren = {[]};
            tree.boxCorners = [-1 -1 -1 1 1 1]/2;
            tree.boxPoints = {(1:tree.N)'};
            tree.pointBoxes = ones(tree.N, 1);
            tree.numBoxes = 1;
            tree.verbose = false;
            % Recursive divide
            divide_box(tree, 1);
            % Finish up
            tree.create_colleague_lists();
        end

        function divide_box(tree, boxNo)
            if tree.boxLevels(boxNo)==tree.maxLevel
                if tree.verbose
                    fprintf("%5d: leaf with %3d points\n", boxNo, numel(tree.boxPoints{boxNo}));
                end
                return
            end
            boxOrigin = tree.boxCorners(boxNo, 1:3);
            boxSize = tree.boxCorners(boxNo, 4:6) - boxOrigin;
            newboxSize = boxSize/2;
            myChildren = zeros(2, 2, 2);
            myPoints = tree.boxPoints{boxNo};
            % Create children
            for k=1:2
                for j=1:2
                    for i=1:2
                        newBoxNo = tree.numBoxes + 1;
                        tree.numBoxes = newBoxNo;
                        myChildren(i, j, k) = newBoxNo;
                        newBoxOrigin = boxOrigin + newboxSize.*([i j k]-1);
                        tree.boxLevels(newBoxNo) = tree.boxLevels(boxNo)+1;
                        tree.boxParents(newBoxNo) = boxNo;
                        tree.boxChildren{newBoxNo} = [];
                        % Set new boundaries
                        newBoxCorners = [newBoxOrigin, newBoxOrigin+newboxSize];
                        tree.boxCorners(newBoxNo, :) = newBoxCorners;
                        % Divide points
                        mask = tree.points(myPoints, 1) >= newBoxCorners(1) & ...
                               tree.points(myPoints, 2) >= newBoxCorners(2) & ...
                               tree.points(myPoints, 3) >= newBoxCorners(3) & ...
                               tree.points(myPoints, 1) < newBoxCorners(4) & ...
                               tree.points(myPoints, 2) < newBoxCorners(5) & ...
                               tree.points(myPoints, 3) < newBoxCorners(6);
                        newBoxPoints = myPoints(mask);
                        tree.boxPoints{newBoxNo} = newBoxPoints;
                        assert(all(tree.pointBoxes(newBoxPoints)==boxNo)); % Assert parent owns these points
                        tree.pointBoxes(newBoxPoints) = newBoxNo;
                    end
                end
            end
            assert(all((tree.pointBoxes(myPoints)==boxNo)==false))
            tree.boxChildren{boxNo} = myChildren;
            tree.boxPoints{boxNo} = [];
            % Recursively subdivide children
            for i=1:8
                child = tree.boxChildren{boxNo}(i);
                divide_box(tree, child);
            end
        end

        function create_colleague_lists(tree)
            tree.boxColleagues = cell(1, tree.numBoxes);
            tree.boxColleagueShifts = cell(1, tree.numBoxes);
            if tree.periodic
                tree.boxColleagues{1}      = [1];
                % 3x3x3 periodic shifts of root box
                [sx,sy,sz] = ndgrid(-1:1);
                S = [sx(:) sy(:) sz(:)];
                tree.boxColleagues{1} = ones(1, 27);
                tree.boxColleagueShifts{1} = S;
            else
                tree.boxColleagues{1} = 1;
                tree.boxColleagueShifts{1} = [0 0 0];
            end
            for l=1:tree.maxLevel
                rl = 1/2^l;
                level_boxes = find(tree.boxLevels == l);
                for i=1:numel(level_boxes)
                    box = level_boxes(i);
                    parent = tree.boxParents(box);
                    parentColleagues = tree.boxColleagues{parent};
                    parentColleagueShifts = tree.boxColleagueShifts{parent};
                    possible_colleagues = zeros(1, 27*8);
                    possible_shifts = zeros(27*8, 3);
                    n = 0;
                    % Parent's colleagues' children are possible colleagues
                    for pc_idx = 1:numel(parentColleagues)
                        pc = parentColleagues(pc_idx);
                        pc_shift = parentColleagueShifts(pc_idx, :);
                        pc_children = reshape(tree.boxChildren{pc}, 1, []);
                        m = numel(pc_children);
                        possible_colleagues((n+1):(n+m)) = pc_children;
                        possible_shifts((n+1):(n+m), :) = repmat(pc_shift, m, 1);
                        n = n+m;
                    end
                    possible_colleagues = possible_colleagues(1:n);
                    possible_shifts = possible_shifts(1:n, :);
                    shifted_centers = tree.box_center(possible_colleagues) + possible_shifts;
                    dist = max(abs( tree.box_center(box) - shifted_centers),...
                               [], 2);
                    % Colleague center-to-center is within r_l + safety factor
                    mask = dist < rl*1.000001;
                    tree.boxColleagues{box} = possible_colleagues(mask);
                    tree.boxColleagueShifts{box} = possible_shifts(mask, :);
                end
            end
        end

        function idx_list = nonempty_leafs(tree)
            mask = false(1, tree.numBoxes);
            for i=1:tree.numBoxes
                mask(i) = ~isempty(tree.boxPoints{i});
            end
            idx_list = find(mask);
        end

        function c = box_center(tree, boxNo)
            c = (tree.boxCorners(boxNo(:), 4:6) + tree.boxCorners(boxNo(:), 1:3)) / 2;
        end
        
        function s = box_size(tree, boxNo)
            s = tree.boxCorners(boxNo, 4:6) - tree.boxCorners(boxNo, 1:3);            
        end

        function idxs = box_point_idxs(tree, boxNo)
            if isempty(tree.boxChildren{boxNo})
                idxs = tree.boxPoints{boxNo};
                return
            end
            child_idxs = cell(1, 8);
            for c=1:8
                child = tree.boxChildren{boxNo}(c);
                child_idxs{c} = tree.box_point_idxs(child);
            end
            idxs = vertcat(child_idxs{:});
        end
        
        function points = box_points(tree, boxNo)
            idxs = tree.box_point_idxs(boxNo);
            points = tree.points(idxs, :);
        end

        function grid = box_grid(tree, boxNo, rvec)
        % Return box tensor product grid with sides rvec
            [xp, yp, zp] = ndgrid(rvec, rvec, rvec);
            c = tree.box_center(boxNo);
            s = tree.box_size(boxNo);
            grid = [xp(:), yp(:), zp(:)].*s/2 + c;
        end
    end
end
