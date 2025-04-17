clear

c = 15;
K = c;

psi = pswf(0, c);
lambda = integral(psi)/psi(0);

psi = psi/lambda/psi(0);

dpsi = diff(psi);

psi_hat = chebfun(@(k) psi(k/c)*lambda, [-K K]);
dpsi_hat = diff(psi_hat);

psi_hat2 = ftrans(psi, K);
psi2 = invftrans(psi_hat, K);

kcheb = chebfun(@(k) k, [-K K]);
gamma_hat = psi_hat - kcheb/2*dpsi_hat;
gamma2 = invftrans(gamma_hat, K);

rcheb  = chebfun(@(r) r);
gamma = 3/2*psi + 1/2*rcheb*dpsi;

kernel = kernels.stresslet_pswf(c_pswf=c);

sfigure(1);
clf
subplot(1,2,1);
plot(psi)
hold on
plot(psi2,'.')
plot(gamma2)
plot(gamma,'*')
xlabel('x')

subplot(1,2,2)
plot(psi_hat)
hold on
plot(psi_hat2, '.')
xlabel('k')
k = linspace(0, c, 15);
plot(k, kernel.fourier_scaling(k.^2, 0),'*')
plot(gamma_hat)


return


sfigure(2);
clf

Bres_cheb = compute_Bres(gamma_hat, c);
plot(Bres_cheb)
hold on
rcheb_01 = chebfun(@(r) r, [0 1]);
plot(rcheb_01)
plot(rcheb_01-Bres_cheb)

psi_01 = chebfun(psi, [0 1]);
gamma_01 = chebfun(gamma, [0 1]);

erf_psi = cumsum(psi_01);
erf_gamma = cumsum(gamma_01);
% = erf_psi + 1/2*rcheb*psi_01

f1 = cumsum(erf_gamma);
f1_2 = rcheb_01*erf_psi - 1/2*cumsum(rcheb_01*psi_01);

f2 = cumsum(f1);
f2_2 = 1/2*(rcheb_01^2*erf_psi - rcheb_01*cumsum(rcheb_01*psi_01));

f2 = f2 / rcheb_01 * 4;

Bm = 2*(rcheb_01*cumsum(psi_01) - cumsum(rcheb_01*psi_01));
foo = 2*integral(rcheb_01*psi_01)
const = 1-Bm(1);
plot(Bm+const, '*')
plot(f2-f2(1)+1, '.')
plot(rcheb_01-(f2-f2(1)+1), '.')
grid on

function F=ftrans(f,K)
    xcheb = chebfun(@(x) x);
    ker = @(k) integral(f*exp(-1i*k*xcheb));
    F = real(chebfun(ker, [-K K]));
end

function f=invftrans(F, K)
    kcheb = chebfun(@(k) k, [-K K]);
    ker = @(x) integral(F*exp(1i*kcheb*x)/(2*pi));
    f = real(chebfun(ker, [-1 1]));
end

function Bres_cheb = compute_Bres(gamma_hat, K)
    k = linspace(0, K)';
    hk = k(2)-k(1);
    R = sqrt(3)+1;
    [Bwin, alpha, beta] = windowed_biharmonic(k, R, 1);
    Bmollhat = Bwin .* gamma_hat(k);
    % Fourier integrals for evaluating mollified biharmonic
    f = @(k,r)  k./r   .*(sin(k.*r));
    Bmoll = @(r) 4*pi*hk*sum(  f(k,r) .* Bmollhat, 1) / (2*pi)^3;
    Bres = @(r) r + alpha + beta*r.^2 - Bmoll(r);
    Bres_cheb = chebfun(Bres, [0 1], 'turbo');
end
