%% Main function to generate tests
function tests = test_octree
tests = functiontests(localfunctions);
end

%% Test Functions


function test_colleagues_direct(testCase)
% Compare colleague list (from knnsearch) to direct list
    tree = octree([0,0,0], 3);
    centers = tree.box_center(1:tree.numBoxes);
    for l=0:tree.maxLevel
        rl = 1/2^l;
        level_mask = (tree.boxLevels == l);
        level_centers = centers(level_mask, :);
        level_boxes = find(level_mask);
        for i=1:numel(level_boxes)
            box = level_boxes(i);
            dist = max(abs(level_centers - level_centers(i, :)), [], 2);
            A = sort(level_boxes(dist <= rl+eps()));
            B = sort(tree.boxColleagues{box});
            assert(all(A==B));
        end
    end
    
end

function test_colleagues_are_near(testCase)
    tree = octree([0,0,0], 3);
    for box=1:tree.numBoxes
        Clist = tree.boxColleagues{box};
        cb = tree.box_center(box);
        s = tree.box_size(box);
        rl = 1/2^tree.boxLevels(box);
        assert(all(abs(s-rl)<eps()));
        for c=Clist
            cc = tree.box_center(c);
            d = max(abs(cb-cc));
            assert(d <= rl+eps());
        end
    end
end
