function Tp2c = operator_parent2child(p)
% Precompute matrices for interpolation
    [rvec, V] = approx.chebvander(p);
    Vi = inv(V);
    U = cell(2, 1);
    x1 = (rvec-1)/2; % [-1, 0]
    x2 = (rvec+1)/2; % [ 0, 1]
    U{1} = Vi*approx.chebevalmat(x1, p);
    U{2} = Vi*approx.chebevalmat(x2, p);
    
    function child_expansion = apply(i, j, k, parent_expansion)
    % Apply operator for child at index (i,j,k) in 2^3 cube
        Ux = U{i};
        Uy = U{j};
        Uz = U{k};
        dim = size(parent_expansion, 2);
        child_expansion = zeros(p^3, dim);
        for d=1:dim
            child_expansion(:,d) = approx.kronmat3_apply(Uz, Uy, Ux, parent_expansion(:,d));
        end
    end
    Tp2c = @apply;
end

