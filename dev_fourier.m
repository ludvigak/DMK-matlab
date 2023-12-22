clear all
rng(1)
format compact

% Target tolerance
tol = 1e-9;

% Run experiments on tree top level
l = 0;
rl = 1/2^l;
sigma_0 = 1/sqrt(log(1/tol));
sigma_1 = sigma_0 / 2;

% Set up random sources in box [-1/2, 1/2]^3
N = 100;
points = rand(N, 3)-1/2;
charges = rand(N, 1)-1/2;

% Put charges on x-axis at box edges, to measure decay
points(1:2,:) = [ 0.5 0 0
                  -0.5 0 0];
charges(1:2) = 1;    

% Can also study a single source at origin
    
% Measurement points on x-axis
x = linspace(-1.5,1.5, 200);
y = 0.0;
z = 0.0;

% Compute difference kernel D0
r = sqrt( (points(:, 1)-x).^2 + ...
          (points(:, 2)-y).^2 + ...
          (points(:, 3)-z).^2 );
D0 = (erf(r/sigma_1) - erf(r/sigma_0))./r;
% u_diff computed directly
udiff_direct = charges'*D0;

% % Fourier params

% FUDGE = 1.5; % Why??
%     ;
% nf = ceil( 3/pi*log(1/tol) *FUDGE)
% hl = 4*pi/3 /FUDGE


D = 3*rl;
hl = 2*pi/D

Kl = 4/rl * log(1/tol)

nf = ceil(Kl/hl)

Kmax = hl*nf
return
% In Table 3.1 N_1^G is stated as 66 for tol=1e-9
N1 = 2*nf+1 

% Setup Fourier vectors
[m1, m2, m3] = ndgrid(-nf:nf, -nf:nf, -nf:nf);
k1 = hl*m1(:);
k2 = hl*m2(:);
k3 = hl*m3(:);
ksq = k1.^2 + k2.^2 + k3.^2;

% Difference kernel Fourier transform
D0hat = 4*pi*(exp(-ksq*sigma_1^2/4) - exp(-ksq*sigma_0^2/4))./ksq;
% k=0 term
D0hat(k1==0 & k2==0 & k3==0) = pi*(sigma_0^2-sigma_1^2);
% w_l factor 
wl = 1/(2*pi)^3 * D0hat;

% I think this was missing in the paper, but is needed somewhere
wl = wl * hl^3;

% Source expansion
kdoty = k1.*points(:, 1)' + k2.*points(:, 2)' + k3.*points(:, 3)';
ghat = exp(-1i*kdoty) * charges;

% Evaluate at sample points
kdotx = k1.*x + k2.*y + k3.*z;
udiff_fourier = (wl .* ghat).' * exp(1i*kdotx);

fourier_eval_max_err = norm(udiff_fourier - udiff_direct, inf)


% Anterpolate to proxy points
p = 45;
[rvec, V] = approx.chebvander(p);
ViT = transpose(inv(V));
scaled_points = points*2;
proj = approx.chebevalmat3_trans_apply(...
    scaled_points(:, 1), ...
    scaled_points(:, 2), ...
    scaled_points(:, 3), ...
    p, charges);
proxy_charges = approx.kronmat_apply(ViT, proj, 3);
% [xp, yp, zp] = ndgrid(rvec, rvec, rvec);
% xp = xp(:); yp = yp(:); zp = zp(:);

% kdoty = k1.*xp' + k2.*yp' + k3.*zp';
% ghat2 = exp(-1i*kdoty) * proxy_charges;
tic
k = hl*(-nf:nf)';
M = exp(-1i*k.*rvec'/2);
ghat3 = approx.kronmat3_apply(M, M, M, proxy_charges);
toc

udiff_fourier2 = (wl .* ghat3).' * exp(1i*kdotx);
fourier2_eval_max_err = norm(udiff_fourier2 - udiff_direct, inf)


% kron_err = norm(ghat2-ghat3, inf)
g_err = norm(wl.*(ghat3-ghat), inf)


figure(1)
clf
subplot(2, 1, 1)
plot(x, real(udiff_fourier))
hold on
plot(x, udiff_direct, '--')
xlabel('x')
ylabel('u_diff')
legend('fourier', 'direct')

subplot(2, 1, 2)
plot(x, abs(real(udiff_fourier)))
hold on
plot(x, abs(udiff_direct), '--')
plot(x, abs(udiff_direct - udiff_fourier), ':', 'linewidth', 1.5)
plot(x, abs(udiff_direct - udiff_fourier2), ':', 'linewidth', 1.5)
xlabel('x')
ylabel('u_diff')
legend('fourier', 'direct', 'error', 'error via proxy')
set(gca, 'yscale', 'log')

figure(2)
clf

Nf = 2*nf+1;
Uhat = reshape( (wl .* ghat), Nf, Nf, Nf);
U = fftshift(ifftn((Uhat))) * Nf^3;


%[X, Y] = 

pcolor(log10(abs((U(:, :, 41)))))
shading interp
colorbar
