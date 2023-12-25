function y = kronmat3_apply(A, B, C, x)
% Fast evaluation of 3-dimensional Kronecker tensor product matrix to vector x
%     y = (A \otimes B \otimes C) * x
% using level 3 BLAS operations as much as possible
    [m1, n1] = size(A);
    [m2, n2] = size(B);
    [m3, n3] = size(C);
    X = reshape(x, n3, n2, n1);
    K = pagemtimes(C, X);
    K = pagemtimes(K, 'none', B, 'transpose');
    Y = reshape(K, m2*m3, n1) * A.';
    y = reshape(Y, m1*m2*m3, 1);  
end
