function dpsi = pswf0_diff(c, n)
% Form chebfun of pswf derivative using integral formula
    psi = pswf(0, c);
    lambda0 = sum(psi, -1, 1) / psi(0);
    t = chebfun(@(t) t);
    function y = fdpsi(x)
        y = zeros(size(x));
        for i=1:numel(x)
            y(i) = real( sum((1i*c*t)^n*psi*exp(1i*c*x(i)*t)) );
        end
        y = y / lambda0;
    end
    dpsi = chebfun(@fdpsi);
end

