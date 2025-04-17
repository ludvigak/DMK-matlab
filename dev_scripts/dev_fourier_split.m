clear all
tol = 1e-14;

kernel = kernels.laplace_ewald(tolerance=tol);

r = linspace(-1.0, 1.0, 200).';

% Setup root level windowed kernel
Kmax_win = ceil( 2.5*log(1/tol));
nf_win = Kmax_win;
hf_win = Kmax_win/nf_win;
Ctrunc = sqrt(3) + 6*kernel.sigma_0;    

k = hf_win*(-nf_win:nf_win);


% LAPLACE radial integral
GR = kernel.reskernel([r r*0 r*0], [0 0 0], 1, 0);
k = linspace(0, Kmax_win, nf_win);
hk = k(2)-k(1);
W0hat_rad_op = kernel.winkernel_fourier(k, k*0, k*0, Ctrunc);
W0radhat = W0hat_rad_op(ones(size(k)));
MRrad = 4*pi*hk*sum(sin(k.*r) .* k .* W0radhat, 2) ./ r / (2*pi)^3;
GRrad = 1./abs(r) - MRrad;
laplace_err = norm(GRrad-GR, inf)

% LAPLACE 3D integral
p = 100;
kbase = hf_win*(-nf_win:nf_win);
[k1, k2, k3] = ndgrid(kbase);
k1 = k1(:); k2 = k2(:); k3 = k3(:);
ksq = k1.^2 + k2.^2 + k3.^2;
w = hf_win^3 / (2*pi)^3;
Nf = numel(k1);
W0hat = kernel.winkernel_fourier(k1, k2, k3,  Ctrunc);
Psi = ones(size(k1));
GRhat = w*W0hat(Psi(:));
[rvec, V] = approx.chebvander(p);
kdotr = exp(-1i*kbase(:).*rvec');
M = V\(kdotr');
GRexpa = real( approx.kronmat_apply(M, GRhat, 3) );
GRint = 1./abs(r) - approx.chebevalmat3_apply(r, r*0, r*0, p, GRexpa);
laplace_3dint_err = norm(GRint-GR, inf)
clf
plot(r, GRint)
hold on
plot(r, GR)




% STOKES
r = linspace(-1.5, 1.5, 200).';

rabs = abs(r);
skernel = kernels.stokeslet_hasimoto(tolerance=1e-14);
sigma_0 = kernel.sigma_0;
xi = 1/sigma_0;
[Sdiag, Soffd] = skernel.real_decay(rabs, 0);

ksq = k.^2;
Bwin = biharmonic_trunc(k, Ctrunc);

Bmollhat = Bwin .* (1 + ksq/(4*xi^2)) .*  exp(-ksq/4/xi^2);

%Bmollhat = Bwin .* exp(-ksq.^2/4/xi^2/Kmax_win^2);

skernel_pswf = kernels.stokeslet_pswf(c_pswf=40);

c_pswf=40;
psi = pswf(0, c_pswf);
dpsi = diff(psi);
d2psi = diff(dpsi);
alpha = -d2psi(0)/psi(0)/c_pswf^2/2;


psi_arg = k/c_pswf;
psi_arg(psi_arg > 1) = 1;

%Bmollhat = Bwin .* psi(psi_arg)/psi(0) .* (1 + alpha*k.^2);

%Bmollhat = Bwin .* psi(ksq/max(ksq))/psi(0);



f   = @(k,r)  k./r   .*(sin(k.*r));
df  = @(k,r) -k./r.^2.*(sin(k.*r)-k.*r.*cos(k.*r));
d2f = @(k,r) -k./r.^3.*(sin(k.*r).*-2.0+k.^2.*r.^2.*sin(k.*r)+k.*r.*cos(k.*r).*2.0);


Bmoll   = 4*pi*hk*sum(  f(k,rabs) .* Bmollhat, 2) / (2*pi)^3;
dBmoll  = 4*pi*hk*sum( df(k,rabs) .* Bmollhat, 2) / (2*pi)^3;
d2Bmoll = 4*pi*hk*sum(d2f(k,rabs) .* Bmollhat, 2) / (2*pi)^3;

Bmoll = -(Bmoll - r.^2/2/Ctrunc - Ctrunc/2);
dBmoll = -(dBmoll - rabs/Ctrunc);
d2Bmoll = -(d2Bmoll - 1/Ctrunc);

Bres = rabs - Bmoll;
dBres = 1 - dBmoll;
d2Bres = -d2Bmoll;

Sdiagrad =  rabs .* d2Bres + dBres;
Soffdrad = -rabs .* d2Bres + dBres;


selfi = skernel.self_interaction([1 1 1], 0)

hold on
plot(r, Sdiagrad)
plot(r, Soffdrad)

sfigure(1); clf
subplot(1,2,1)
semilogy(k, abs(f(k,1) .* Bmollhat))
hold on
semilogy(k, abs(df(k,1) .* Bmollhat))
semilogy(k, abs(d2f(k,1) .* Bmollhat))
axis([0 Kmax_win 1e-15 1e2])

subplot(1,2,2)
semilogy(r, abs(Sdiag), r, abs(Sdiagrad), '.')
hold on
semilogy(r, abs(Soffd), r, abs(Soffdrad), 'o')
diagerr = norm(Sdiag-Sdiagrad, inf)
ylim([1e-16 1])

return

idx_in = 1;
idx_out = 2;
S0hat = skernel.winkernel_fourier(k1, k2, k3, Ctrunc);
Psi = zeros(Nf, 3);
Psi(:, idx_in) = 1;
SMhat = w*S0hat(Psi);

% n = 2*nf_win+1;
% SMhat1 = reshape(SMhat(:, idx_out), n, n, n);
% clf
% subplot(1, 3, 1)
% pcolor(log10(abs(SMhat1(:, :, nf_win+1))))
% axis equal
% shading interp
% subplot(1, 3, 2)
% pcolor(log10(abs(squeeze(SMhat1(:, nf_win+1, :)))))
% axis equal
% shading interp
% subplot(1, 3, 3)
% pcolor(log10(abs(squeeze(SMhat1(nf_win+1, :, :)))))
% axis equal
% shading interp
% return

SMexpa = real( approx.kronmat_apply(M, SMhat(:, idx_out), 3) );
SM = approx.chebevalmat3_apply(r, r, r, p, SMexpa);
fsrc = zeros(1, 3); fsrc(idx_in) = 1;
umoll = skernel.mollkernel([r r r], [0 0 0], fsrc, 0);
umoll = umoll(:, idx_out);
udir = skernel.direct([r r r], [0 0 0], fsrc);

clf
subplot(1,2,1)
semilogy(r, abs(udir(:, idx_out)-umoll))
hold on
semilogy(r, abs(udir(:, idx_out)-SM), '.') 
title('stokes R')

subplot(1,2,2)
ksq = k1.^2 + k2.^2 + k3.^2;
semilogy(sqrt(ksq), abs(SMhat), '.')

function Bwin = biharmonic_trunc(K, R)
    K2 = K.^2;
    Bwin = 8*pi*( ...
        (2-R^2*K.^2).*cos(R*K) + 2*R*K.*sin(R*K) - 2 ...
                ) ./ (2*K.^4);
    Bwin(K==0) = pi*R^4;

    Bwin = 8*pi*((1 + 1/2*cos(R*K) - 3/2*sin(R*K)./(R*K)))./K2.^2;
    Bwin(K==0) = pi*R^4/15;
end
