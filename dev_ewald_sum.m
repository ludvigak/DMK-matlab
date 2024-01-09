% Development of Ewald decompositions for Laplace and Stokes
% using both classical methods and PSWF

clear
rng(1);

showplots = true; % Enable/disable to show plots

global c_pswf % WARNING: GLOBAL VARIABLE

% Test system of N points in L^3 cube, measuring everything at center
N =  50;
L = 2*pi;
points = L*rand(N, 3);

target = points(1, :); % With self-interaction
%target = [L L L/2]/2; % Without self-interaction

%% Laplace
disp('* Laplace Classical')
% Setup charge neutral system
charges = rand(N, 1);
charges = charges - sum(charges)/N;
% Check xi independence
xilist = linspace(0.9, 1.5, 10);
ulist = zeros(size(xilist));
for idx=1:numel(xilist)
    xi = xilist(idx);
    u = ewald_sum(target, points, charges, @laplace_real, @laplace_fourier, @laplace_self, L, xi);
    ulist(idx) = u;
end
udiff = ulist - ulist(end);
laplace_udiff_xi = norm(udiff, inf) / norm(u, inf)
assert(laplace_udiff_xi < 1e-13)
u_ewald = ulist(end);

disp('* Laplace PSWF')
c_pswf = 30; % This seems to give us all the accuracy we can get
rl = L;
xi_pswf = 1/rl * c_pswf/2; % Somewhat proportional. This is used throughout
% Check that we get the same value as using regular Ewald
u_pswf = ewald_sum(target, points, charges, @laplace_real_pswf, @laplace_fourier_pswf, @laplace_self_pswf, L, xi_pswf);
err_laplace_pswf = abs(u_pswf - u_ewald)
assert(err_laplace_pswf < 1e-12); 

if showplots
    laplace_plots();
end


%% Stokes
disp('* Stokes: Hasimoto')
strengths = rand(N, 3);
strengths = strengths - sum(strengths, 1)/N;
% Check xi independence
xilist = linspace(0.9, 2, 10);
ulist = zeros(numel(xilist), 3);
for idx=1:numel(xilist)
    xi = xilist(idx);
    u = ewald_sum(target, points, strengths, ...
                  @stokes_real_hasimoto, @stokes_fourier_hasimoto, @stokes_self_hasimoto, L, xi);
    ulist(idx, :) = u;
end
u_hasimoto = ulist(end, :);
udiff = ulist - ulist(end, :);
stokes_udiff_xi = norm(udiff(:), inf) / norm(u, inf)
assert(stokes_udiff_xi < 1e-13)

disp('* Stokes: Beenakker')
% Compare Hasimoto and Beenakker values
u_beenakker = ewald_sum(target, points, strengths, ...
                        @stokes_real_beenakker, @stokes_fourier_beenakker, @stokes_self_beenakker, ...
                        L, xi);

err_beenakker = norm(u_beenakker - u_hasimoto)
assert(err_beenakker < 1e-13)

disp('* Stokes: exp_erf')
% Stokes using the "general window function" framework,
% with the window function being the exponential,
% such that result should be == Bennakker
u_experf = ewald_sum(target, points, strengths, ...
                     stokes_real_window(@exp_erfc), ...
                     stokes_fourier_window(@exp_hat), ...
                     stokes_self_window(@exp_erfc), ...
                     L, xi);
err_experf = norm(u_beenakker - u_experf)
assert(err_experf < 1e-13)

disp('* Stokes: PSWF')
% Finally run Stokes using PSWF.
c_pswf = 40; % We need a bigger c
rl = L;
xi_pswf = 1/rl * c_pswf/2;
u_pswf = ewald_sum(target, points, strengths, ...
                   stokes_real_window(@pswf_erfc), ...
                   stokes_fourier_window(@pswf_hat), ...
                   stokes_self_window(@pswf_erfc), ...
                   L, xi_pswf);
% Compare to previous results
err_pswf = norm(u_hasimoto - u_pswf)
assert(err_pswf < 1e-12)

if showplots
    stokes_plots();
end

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kernel-independent Ewald calculation
function u = ewald_sum(target, points, charges, r_kernel, f_kernel, self_kernel, L, xi)
    % === Real space
    real_check = norm(r_kernel([0 0 0], [L 0 0], charges(1,:).^0, xi));
    assert(real_check < 1e-13) % Sanity check of kernel decay
    uR = 0;
    uF = 0;
    N = size(points, 1);
    fself = [];
    % Collect 3^3 nearest neighbor boxes
    for p1=-1:1
        for p2 = -1:1
            for p3 = -1:1
                shift = L*[p1 p2 p3];
                for i=1:N
                    src = points(i, :)+shift;
                    f = charges(i,:);
                    if norm(src-target)==0
                        fself = f; % We have self interaction;
                        continue
                    end
                    uR = uR + r_kernel(target, src, f, xi);
                end
            end
        end
    end
    % === Fourier space
    % This parameter choice is for exponential, too much for PSWF
    nf = ceil(L/(2*pi) * 2*xi*sqrt(log(1/eps)));
    kvec = (-nf:nf)*2*pi/L;
    kmax = kvec(end);
    fourier_check = norm(f_kernel(kmax, 0, 0, xi, charges(1,:).^0));
    assert(fourier_check < 1e-13); % Sanity check Fourier decay
    [k1, k2, k3] = ndgrid(kvec);
    k1 = k1(:); k2 = k2(:); k3 = k3(:);
    for i=1:N
        rvec = target-points(i,:);
        kr = k1*rvec(1) + k2*rvec(2) + k3*rvec(3);
        GF = f_kernel(k1, k2, k3, xi, charges(i,:));
        uhat = exp(-1i*kr) .* GF;
        uF = uF + 1/L^3 * sum(real(uhat), 1);
    end
    if ~isempty(fself)
        uS = self_kernel(fself, xi);
    end
    u = uR + uF + uS;
end

%% Laplace kernel

% Classical Ewald split: Real
function uR = laplace_real(trg, src, rho, xi)
    r = sqrt(sum((trg-src).^2));
    uR = rho*erfc(xi*r)/r;
end
% Classical Ewald split: Fourier
function uF = laplace_fourier(k1, k2, k3, xi, charge)
    ksq = k1.^2 + k2.^2 + k3.^2;
    uF = 4*pi./ksq .*exp(-ksq/4/xi^2) * charge;
    uF(ksq==0) = 0;
end

function uS = laplace_self(fself, xi)
    uS = -fself*2*xi/sqrt(pi);
end

% PSWF real
function uR = laplace_real_pswf(trg, src, rho, xi)
    r = sqrt(sum((trg-src).^2));
    uR = rho*pswf_erfc(r, xi)/r;
end
% PSWF Fourier
function uF = laplace_fourier_pswf(k1, k2, k3, xi, charge)
    ksq = k1.^2 + k2.^2 + k3.^2;
    kabs = sqrt(ksq);
    psi_hat = pswf_hat(kabs, xi);
    uF = 4*pi./ksq .* psi_hat * charge;
    uF(ksq==0) = 0;
end

function uS = laplace_self_pswf(fself, xi)
    [~, dPhi] = pswf_erf(0, xi);
    uS = -fself*dPhi;
end

%% General window functions
% These functions define an erf-like decay in real space,
% based on a window function of choice

% === Exponential Window
function [Phi, dPhi, ddPhi] = exp_erf(r, xi)
    Phi = erf(xi*r);
    dPhi = 2 * xi * exp(-r.^2*xi^2) / sqrt(pi);
    ddPhi = -4 * xi^3 * r .* exp(-r.^2*xi^2) / sqrt(pi);
end

function [psi_hat, dpsi_hat, ddpsi_hat] = exp_hat(k, xi)
    psi_hat   = exp(-k.^2/(4*xi^2));
    dpsi_hat  = -(k.*exp(-k.^2/(4*xi^2)))/(2*xi^2);
    ddpsi_hat = (k.^2.*exp(-k.^2/(4*xi^2)))/(4*xi^4) - exp(-k.^2/(4*xi^2))/(2*xi^2);
end

function [Phic, dPhic, ddPhic] = exp_erfc(r, xi)
    [Phi, dPhi, ddPhi] = exp_erf(r, xi);
    Phic = 1 - Phi;
    dPhic = -dPhi;
    ddPhic = -ddPhi;
end

% === PSWF Window
function [Phi, dPhi, ddPhi] = pswf_erf(r, xi)
    global c_pswf % WARNING: GLOBAL VARIABLE
    psi  = pswf(0, c_pswf);
    dpsi = diff(psi);
    c0 = integral(psi, 0, 1);
    rl = 1/xi * c_pswf/2;
    if r/rl > 1
        Phi = 1;
        dPhi = 0;
        ddPhi = 0;
    else
        Phi = 1/c0*integral(psi, 0, r/rl);
        dPhi = psi(r/rl)/(c0*rl);
        ddPhi = dpsi(r/rl)/(c0*rl^2);
    end
end

function [psi_hat, dpsi_hat, ddpsi_hat] = pswf_hat(k, xi)
    % Return scaled function \hat\psi(k*r_l/c) / \psi(0), and derivatives
    global c_pswf % WARNING: GLOBAL VARIABLE
    psi   = pswf(0, c_pswf);
    dpsi  = diff(psi);
    ddpsi = diff(dpsi);
    rl = 1/xi * c_pswf/2;
    psi_arg = k*rl/c_pswf;
    [psi_hat, dpsi_hat, ddpsi_hat] = deal(zeros(size(k)));
    mask = (psi_arg <= 1);
    psi_hat(mask)   =   psi(k(mask)*rl/c_pswf)/psi(0);
    dpsi_hat(mask)  =  dpsi(k(mask)*rl/c_pswf)/psi(0) * rl/c_pswf;
    ddpsi_hat(mask) = ddpsi(k(mask)*rl/c_pswf)/psi(0) * (rl/c_pswf)^2;
end

function [Phic, dPhic, ddPhic] = pswf_erfc(r, xi)
    [Phi, dPhi, ddPhi] = pswf_erf(r, xi);
    Phic = 1 - Phi;
    dPhic = -dPhi;
    ddPhic = -ddPhi;
end

%% Stokes kernel

% Decompositions
function [Sdiag, Soffd] = stokes_real_decay_hasimoto(r, xi)
    Sdiag = -2*xi/sqrt(pi)*exp(-xi^2*r^2) + erfc(xi*r)/r;
    Soffd =  2*xi/sqrt(pi)*exp(-xi^2*r^2) + erfc(xi*r)/r;
end

function B = stokes_fourier_scaling_hasimoto(ksq, xi)
    B = 8*pi./ksq.^2 .* (1 + ksq/(4*xi^2)) .*  exp(-ksq/4/xi^2);
    B(ksq==0) = 0;
end

function uS = stokes_self_hasimoto(fself, xi)
    C = -4*xi/sqrt(pi);
    uS =  C * fself;
end

function [Sdiag, Soffd] = stokes_real_decay_beenakker(r, xi)
    Sdiag = erfc(xi*r)/r - 6*xi/sqrt(pi)*exp(-xi^2*r^2) + 4*xi^3*r^2/sqrt(pi)*exp(-xi^2*r^2);
    Soffd = erfc(xi*r)/r + 2*xi/sqrt(pi)*exp(-xi^2*r^2) - 4*xi^3*r^2/sqrt(pi)*exp(-xi^2*r^2);
end

function B = stokes_fourier_scaling_beenakker(ksq, xi)
    B = 8*pi./ksq.^2 .* (1 + ksq/(4*xi^2) + ksq.^2/(8*xi^4)) .*  exp(-ksq/4/xi^2);
    B(ksq==0) = 0;
end

function uS = stokes_self_beenakker(fself, xi)
    C = -8*xi/sqrt(pi);
    uS =  C * fself;
end

% Using a general window function
function [Sdiag, Soffd] = stokes_real_decay_window(r, xi, win_erfc)
    [Phic, dPhic, ddPhic] = win_erfc(r, xi);
    %Theta = r*Phic; % Not used
    dTheta = Phic + r*dPhic;
    ddTheta = dPhic + dPhic + r*ddPhic;
    Sdiag =  ddTheta + dTheta/r;
    Soffd = -ddTheta + dTheta/r;
end
function B = stokes_fourier_scaling_window(ksq, xi, win_hat)
    k = sqrt(ksq);
    [psi_hat, dpsi_hat, ddpsi_hat] = win_hat(k, xi);
    B = 8*pi./ksq.^2 .* ( ksq.*ddpsi_hat/2 - k.*dpsi_hat + psi_hat );
    B(ksq==0) = 0;
end
function C = stokes_self_scaling_window(xi, win_erfc)
    [~, dPhic, ~] = win_erfc(0, xi);
    C = 4 * dPhic;
end

% === Define functions for use with Ewald summation
function uR = stokes_real_hasimoto(trg, src, f, xi)
    uR = stokes_real(trg, src, f, xi, @stokes_real_decay_hasimoto);
end
function uF = stokes_fourier_hasimoto(k1, k2, k3, xi, f)
    uF = stokes_fourier(k1, k2, k3, xi, f, @stokes_fourier_scaling_hasimoto);
end
function uR = stokes_real_beenakker(trg, src, f, xi)
    uR = stokes_real(trg, src, f, xi, @stokes_real_decay_beenakker);
end
function uF = stokes_fourier_beenakker(k1, k2, k3, xi, f)
    uF = stokes_fourier(k1, k2, k3, xi, f, @stokes_fourier_scaling_beenakker);
end
% General window: real
function real_fun = stokes_real_window(win_erfc)
    real_decay = @(r, xi) stokes_real_decay_window(r, xi, win_erfc);
    real_fun = @(trg, src, f, xi) stokes_real(trg, src, f, xi, real_decay);
end
% General window: Fourier
function fourier_fun = stokes_fourier_window(win_hat)
    fourier_scaling = @(ksq, xi) stokes_fourier_scaling_window(ksq, xi, win_hat);
    fourier_fun = @(k1, k2, k3, xi, f) stokes_fourier(k1, k2, k3, xi, f, fourier_scaling);
end
% General window: real
function self_fun = stokes_self_window(win_erfc)
    self_scaling = @(xi) stokes_self_scaling_window(xi, win_erfc);
    self_fun = @(fself, xi) fself * self_scaling(xi);
end


% === Helper routines for vector components
function uR = stokes_real(trg, src, f, xi, decay_fun)
    assert(numel(f)==3);
    rvec = trg-src;
    r = sqrt(sum(rvec.^2));
    rhat = rvec/r;
    [Sdiag, Soffd] = decay_fun(r, xi);
    S = zeros(3,3);
    for i=1:3
        for j=1:3
            S(i,j) = Sdiag*(i==j) + Soffd * rhat(i) * rhat(j);
        end
    end
    uR = (S*f(:)).';
end
function uF = stokes_fourier(k1, k2, k3, xi, f, scaling_fun)
    assert(numel(f)==3);
    ksq = k1.^2 + k2.^2 + k3.^2;
    B = scaling_fun(ksq, xi);
    Nf = numel(k1);
    % Layout is suboptmal in memory,
    % just trying to be readable
    S = zeros(3, 3, Nf);
    k = {k1, k2, k3};
    for i=1:3
        for j=1:3
            S(i, j, :) = ksq*(i==j) - k{i}.*k{j};
        end
    end
    uF = pagemtimes(S, f(:)); % size (3, 1, Nf)
    uF = squeeze(uF).'; % size (3, 1, Nf) -> (3, Nf) -> (Nf, 3)
    uF = uF .* B;
end


%% MISC PLOTTING
function laplace_plots()
    % Plot things Laplace
    L = 1;
    xi_ewald = 2*pi/L;
    xi_pswf = 4*xi_ewald;
    kmax = 2*xi_ewald*sqrt(log(1/eps));
    kmax_pswf = 2*xi_pswf;
    % Make sure last point of PSWF bandwidth is included
    kvec = [linspace(0.1, kmax_pswf, 50) kmax_pswf:kmax_pswf/50:kmax].';
    figure(1); clf
    subplot(1, 2, 1)
    semilogy(kvec, abs(laplace_fourier(kvec, kvec*0, kvec*0, xi_ewald, 1)), displayname="Laplace Fourier Ewald")
    hold on
    semilogy(kvec, abs(laplace_fourier_pswf(kvec, kvec*0, kvec*0, xi_pswf, 1)), '.-', displayname="Laplace Fourier PSWF")
    axis([0 kvec(end) eps 1e5])
    xlabel('|k|')
    grid on
    legend
    subplot(1, 2, 2)
    r = linspace(0, L, 50).';
    real_ewald = 0*r;
    real_pswf = 0*r;
    for i=1:numel(r)
        real_ewald(i) = norm(laplace_real([0 0 0], [r(i) 0 0], 1, xi_ewald));
        real_pswf(i)= norm(laplace_real_pswf([0 0 0], [r(i) 0 0], 1, xi_pswf));
    end
    semilogy(r, real_ewald, displayname="Laplace Real Ewald")
    hold on
    semilogy(r, real_pswf, '.-', displayname="Laplace Real PSWF")
    axis([0 L eps inf])
    legend
    grid on
end    

function stokes_plots()
    % Plot things Stokes
    global c_pswf % WARNING: GLOBAL VARIABLE
    c_pswf = 40;
    L = 1;
    xi_ewald = 2*pi/L;
    xi_pswf = 4*xi_ewald;
    real_pswf_fun = stokes_real_window(@pswf_erfc);
    fourier_pswf_fun = stokes_fourier_window(@pswf_hat);    
    kmax = 2*xi_ewald*sqrt(log(1/eps));
    kmax_pswf = 2*xi_pswf;
    % Make sure last point of PSWF bandwidth is included
    kvec = [linspace(0.1, kmax_pswf, 50) kmax_pswf:kmax_pswf/50:kmax].';    
    fourier_hasimoto = max(abs(...
        stokes_fourier_hasimoto(kvec, kvec*0, kvec*0, xi_ewald, [1,1,1])), [], 2);
    fourier_beenakker = max(abs(...
        stokes_fourier_beenakker(kvec, kvec*0, kvec*0, xi_ewald, [1,1,1])), [], 2);
    fourier_pswf = max(abs(...
        fourier_pswf_fun(kvec, kvec*0, kvec*0, xi_pswf, [1,1,1])), [], 2);
    figure(2); clf
    subplot(1, 2, 1)
    semilogy(kvec, fourier_beenakker, displayname="Stokes Fourier Beenakker")
    hold on
    semilogy(kvec, fourier_hasimoto, '-', displayname="Stokes Fourier Hasimoto")
    semilogy(kvec, fourier_pswf, '.-', displayname="Stokes Fourier PSWF")
    stem(kmax_pswf, fourier_pswf(find(kvec==kmax_pswf, 1)), displayname="PSWF bandwidth")
    axis([0 kvec(end) eps 1e5])
    xlabel('|k|')
    grid on
    legend
    subplot(1, 2, 2)
    rl = 1/xi_pswf * c_pswf/2;
    r = [linspace(0, rl, 50) rl:rl/50:L].';
    [real_hasimoto, real_beenakker, real_pswf] = deal(0*r);
    for i=1:numel(r)
        real_hasimoto(i) = norm(stokes_real_hasimoto([0 0 0], [r(i) 0 0], [1 1 1], xi_ewald));
        real_beenakker(i)= norm(stokes_real_beenakker([0 0 0], [r(i) 0 0], [1 1 1], xi_ewald));
        real_pswf(i)= norm(real_pswf_fun([0 0 0], [r(i) 0 0], [1 1 1], xi_pswf));
    end
    semilogy(r, real_beenakker, displayname="Stokes Real Beenakker")
    hold on
    semilogy(r, real_hasimoto, '-', displayname="Stokes Real Hasimoto")
    semilogy(r, real_pswf, '.-', displayname="Stokes Real PSWF")
    stem(rl, real_pswf(find(r==rl, 1)), displayname="PSWF cutoff")
    axis([0 L eps inf])
    legend
    grid on
end    
