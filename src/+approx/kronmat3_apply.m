function y = kronmat3_apply(A, B, C, x)
% Fast evaluation of 3-dimensional Kronecker tensor product matrix to vector x
% y = (A \otimes B \otimes C) * x
    [m1, n1] = size(A);
    [m2, n2] = size(B);
    [m3, n3] = size(C);
    X = reshape(x, n2*n3, n1);
    K = zeros(m2*m3, n1);
    for i=1:n1
        K(:, i) = approx.kronmat2_apply(B, C, X(:, i));
    end
    Y = K * A';
    y = reshape(Y, m1*m2*m3, 1);  
end
