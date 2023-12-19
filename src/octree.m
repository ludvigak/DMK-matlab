classdef octree < handle
    % Simple linear octree
    
    properties (Access=public)
        N
        maxLevel
        points
        boxLevels
        boxParents
        boxChildren
        boxCorners
        boxPoints
        pointBoxes
        numBoxes
        verbose
    end

    methods
        function tree = octree(points, maxLevel)
        % Points are Nx3 in [-1, 1]^3
            tree.N = size(points, 1);
            tree.maxLevel = maxLevel;
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

        function idx_list = nonempty_leafs(tree)
            mask = false(1, tree.numBoxes);
            for i=1:tree.numBoxes
                mask(i) = ~isempty(tree.boxPoints{i});
            end
            idx_list = find(mask);
        end

        function c = box_center(tree, boxNo)
            c = (tree.boxCorners(boxNo, 4:6) + tree.boxCorners(boxNo, 1:3)) / 2;            
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
