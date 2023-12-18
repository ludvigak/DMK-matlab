clear all

d = 3;
plist1 = [8 9 10 11 12 14 16 18 24 28 32];
plist2 = 2.^[3:8];
time_direct = [];
time_kron   = [];
for p=plist1
    p
    A = rand(p, p);
    x = rand(p^d, 1);
    M = 1;
    for i=1:d
        M = kron(M, A);
    end
    % Kron
    Niter = 1;
    atic = tic;
    for iter=1:Niter
        y1 = M*x;
    end
    time_direct(end+1) = toc(atic) / Niter;
end
M = []
for p=plist2
    p
    A = rand(p, p);
    x = rand(p^d, 1);
    % Fast
    Niter = ceil((128/p)^d)
    atic = tic;
    for iter=1:Niter
        y2 = approx.kronmat_apply(A, x, d);
    end
    time_kron(end+1) = toc(atic) / Niter;
end


clf
loglog(plist1, time_direct, '.-', 'DisplayName', 'Naive')
hold on
loglog(plist2, time_kron, '*-', 'DisplayName', 'KronProd');
loglog(plist1, time_direct(end) * (plist1 ./ plist1(end)).^(2*d), ...
       '--', 'DisplayName', 'O(p^{2d})')
loglog(plist2, time_kron(end) * (plist2 ./ plist2(end)).^(d+1), ...
       '--', 'DisplayName', 'O(p^{d+1})')

grid on
legend('Location', 'best')
xlabel('p')
ylabel('time [s]')
title(sprintf("d=%d", d));
