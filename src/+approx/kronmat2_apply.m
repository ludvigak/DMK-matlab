function y = kronmat2_apply(A, B, x)
% Fast evaluation of 2-dimensional Kronecker tensor product matrix to vector x
% y = (A \otimes B) * x
    [~, n1] = size(A);
    [~, n2] = size(B);
    X = reshape(x, n2, n1);
    Y = B * X * A.';
    y = reshape(Y, [], 1);  
end
