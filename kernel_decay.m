clear all

tol = 1e-5;
sigma_0 = 1/sqrt(log(1/tol));
sigma_1 = sigma_0 / 2;
sigma_2 = sigma_1 / 2;
sigma_3 = sigma_2 / 2;
sigma_4 = sigma_3 / 2;

r = linspace(0, 1, 1000);

clf
semilogy(r, 1./r, 'DisplayName', '1/r');
hold on


R0 = erfc(r/sigma_0)./r;
R1 = erfc(r/sigma_1)./r;
R2 = erfc(r/sigma_2)./r;


M0 = erf(r/sigma_0)./r;
D0 = (erf(r/sigma_1) - erf(r/sigma_0))./r;
D1 = (erf(r/sigma_2) - erf(r/sigma_1))./r;
D2 = (erf(r/sigma_3) - erf(r/sigma_2))./r;

M1 = M0+D0;
%semilogy(r, M0, '--', 'DisplayName', 'M_0');
%semilogy(r, M1, '--', 'DisplayName', 'M_1');
semilogy(r, D0, 'DisplayName', 'D_0');
semilogy(r, D1, 'DisplayName', 'D_1');
semilogy(r, D2, 'DisplayName', 'D_2');

semilogy(r, R0, '-.', 'DisplayName', 'R_0');
semilogy(r, R1, '-.', 'DisplayName', 'R_1');
semilogy(r, R2, '-.', 'DisplayName', 'R_2');


if false
    set(gca, 'yscale', 'linear')
    ylim([tol 20])
else
    ylim([tol 1e5])
end





legend()
grid on
