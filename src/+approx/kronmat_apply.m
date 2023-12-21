function y = kronmat_apply(A, x, d)
% Fast evaluation of d-dimensional Kronecker tensor product matrix to vector x
% y = (A \otimes A ... \otimes A) * x
    [m, n] = size(A);
    if d==2
        X = reshape(x, n, n);
        Y = A * X * A.';
    else
        X = reshape(x, n^(d-1), n);
        K = zeros(m^(d-1), n);
        for i=1:n
            K(:, i) = approx.kronmat_apply(A, X(:, i), d-1);
        end
        Y = K * A.';
    end
    y = reshape(Y, [], 1);  
end
