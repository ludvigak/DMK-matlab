clear all


% r0 = 1;

figure(1)
clf
tol_list = logspace(-1, -13, 100);
tol_plot = 10.^[-3 -6 -9 -12];

tol_list = fliplr(unique(sort([tol_list, tol_plot])));
sources = rand(1, 100)*2-1;
p_first = [];
err_first = [];
tol_achieved = [];
for tol = tol_list

    sigma0 = 1/sqrt(log(1/tol));
    sigma1 = sigma0 / 2;
    M = @(x, y) erf(abs(x-y)/sigma0)./abs(x-y);
    D = @(x, y) erf(abs(x-y)/sigma1)./abs(x-y) - erf(abs(x-y)/sigma0)./abs(x-y);
    K = @(x) sum(D(x, sources), 2);

    
    xeval = rand(100, 1) - 0.5;
    assert(~any(xeval==0))
    Kref = K(xeval);

    plist = 1:50;
    interp_err = 0*plist;
    for i=1:numel(plist)
        p = plist(i);
        [x, V] = approx.chebvander(p);
        x = x/2;
        E = approx.chebevalmat(2*xeval, p);
        Ki = E*(V\K(x));
        interp_err(i) = max(abs(Ki-Kref)) + eps;
    end
    first = find(interp_err < tol, 1);
    if ~isempty(first)
        p_first(end+1) = plist(first);
        tol_achieved(end+1) = tol;
        err_first(end+1) = interp_err(first);
    end
    if any(tol==tol_plot)
        semilogy(plist', interp_err', '.-', 'DisplayName', sprintf('tol=%.0e', tol))
        hold on
        grid on
        legend
    end
end

figure(2);clf
semilogx(tol_achieved, p_first, '.-')
grid on
