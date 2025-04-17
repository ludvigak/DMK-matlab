clear all classes
rng(1);

kernel = @kernels.laplace_ewald;
%kernel = @kernels.laplace_pswf;
%kernel = @kernels.stokeslet_hasimoto;
%kernel = @kernels.stokeslet_pswf;
%kernel = @kernels.stokeslet_pswf2;
%kernel = @kernels.stresslet_hasimoto;
%kernel = @kernels.stresslet_pswf;
%kernel = @kernels.rotlet_ewald;
%kernel = @kernels.rotlet_pswf;

tol = 1e-10;
N = 1000;
max_level = 1;

dmk_opt = dmk_default_opts(tolerance=tol, verbose=true, kernel=kernel, periodic=true);

points = rand(N, 3)-1/2;
charges = rand(N, dmk_opt.kernel.dim_in);
charges = charges - sum(charges, 1)/N;
charges(end, :) = charges(end, :)-sum(charges, 1);
assert(all(abs(sum(charges)) < 1e-14))

disp("* DMK")
atic = tic();clf
dmk_state = dmk_init(points, max_level, dmk_opt);
[u_dmk, ufar, ures, uself] = dmk_apply(charges, dmk_state);
toc(atic)





if false
    p = dmk_opt.p;
    dim_in = dmk_opt.kernel.dim_in;

    u_dmk1 = u_dmk(1)
    u_res1 = ures(1)
    u_far1 = ufar(1)
    u_self1 = uself(1)
    L=1;


    [rvec, V] = approx.chebvander(dmk_opt.p);
    proxy_points = dmk_state.tree.box_grid(1, rvec);
    proxy_charges = init_proxy_charges(dmk_state.tree, charges, p, dmk_state.opt.kernel);

    u_four = zeros(size(u_dmk));
    nf = ceil(L/(2*pi) * 2*sqrt(log(1/eps))/dmk_opt.kernel.sigma_0);
    kvec = (-nf:nf)*2*pi/L;
    kmax = kvec(end);
    [k1, k2, k3] = ndgrid(kvec);
    k1 = k1(:); k2 = k2(:); k3 = k3(:);
    origin = find(k1==0 & k2==0 & k3==0);
    src_points = points;
    src_charges = charges;

    %src_points = proxy_points;
    %src_charges = proxy_charges{1};

    ksq = k1.^2 + k2.^2 + k3.^2;
    ksrc = k1*(src_points(:,1)).' + k2*(src_points(:,2)).' + k3*(src_points(:,3)).';
    ghat1 = exp(1i*ksrc) .* (src_charges).';
    GF = 4*pi./ksq .*exp(-ksq/4*dmk_opt.kernel.sigma_0^2) .* ghat1;
    GF(ksq==0,:) = 0;

    for targ_idx=1:N
        rvec = (points(targ_idx,:)).';
        kr = k1*rvec(1,:) + k2*rvec(2,:) + k3*rvec(3,:);
        uhat = exp(-1i*kr) .* GF;
        u_four(targ_idx) = u_four(targ_idx) + 1/L^3 * sum(sum(real(uhat)));
    end


    src_points = proxy_points;
    src_charges = proxy_charges{1};
    origin = find(k1==0 & k2==0 & k3==0);
    [rvec, V] = approx.chebvander(p);
    kdotr = exp(1i*kvec(:).*rvec'/2);
    ghat2 = approx.kronmat_apply(kdotr, src_charges, 3);
    Mhat = 4*pi./ksq .*exp(-ksq/4*dmk_opt.kernel.sigma_0^2) .* ghat2;
    Mhat(origin) = 0;
    u_four2 = u_four*0;
    for targ_idx=1:N
        t = points(targ_idx,:);
        kr = k1*t(1) + k2*t(2) + k3*t(3);
        u_four2(targ_idx) = sum(real(exp(-1i*kr).*Mhat));
    end


    fourdiff = norm(u_four - u_four2, inf)

    ufar


    u_dmk

    u_per = ures + uself + u_four;
end


disp('* Direct Ewald')
atic = tic();
[u_ewald, ue_moll, ue_res, ue_self] = ewald_sum(points, charges, dmk_opt.kernel);
toc(atic)

error = norm(u_ewald - u_dmk, inf) / norm(u_ewald, inf)

if max_level == 0
    res_err = norm(ue_res - ures, inf)
    self_err = norm(ue_self - uself, inf)
    moll_err = norm(ue_moll - ufar, inf)
end


return

% [sx,sy,sz] = ndgrid([0 -1 1]);
% S = [sx(:) sy(:) sz(:)];
% repeated_points = zeros(27*N, 3);
% repeated_charges = zeros(27*N, dim_in);
% for i=1:27
%     m = N*(i-1) + (1:N);
%     repeated_points(m, :) = points + S(i, :);
%     repeated_charges(m, :) = charges;
% end



% disp('* Direct eval')
% atic = tic();
% u_ref = dmk_opt.kernel.direct(repeated_points, repeated_points, repeated_charges);
% toc(atic)
% u_ref = u_ref(1:N, :);
% max_rel_err = norm(u_ref(:) - u_dmk(:), inf) / norm(u_ref(:), inf)





%%%%%%%%%%%%
L = 1;
target = points(1, :); % With self-interaction
                       %target = [0 0 0]; % Without self-interaction

xi = 1/dmk_opt.kernel.sigma_0
u_ewald = my_ewald_sum(target, points, charges, @laplace_real, @laplace_fourier, @laplace_self, L, xi)


err_dmk = u_ewald - u_dmk(1)
err_per = u_ewald - u_per(1)


return

%% Laplace
disp('* Laplace Ewald')
% Check xi independence
xilist = linspace(0.9, 1.5, 10)*2*pi/L;
ulist = zeros(size(xilist));
for idx=1:numel(xilist)
    xi = xilist(idx);
    u = ewald_sum(target, points, charges, @laplace_real, @laplace_fourier, @laplace_self, L, xi);
    ulist(idx) = u;
end
udiff = ulist - ulist(end);
laplace_udiff_xi = norm(udiff, inf) / norm(u, inf);
assert(laplace_udiff_xi < 1e-13)
u_ewald = ulist(end);
u1 = ulist(1)








%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kernel-independent Ewald calculation
function u = my_ewald_sum(target, points, charges, r_kernel, f_kernel, self_kernel, L, xi)
    % === Real space
    real_check = norm(r_kernel([0 0 0], [L 0 0], charges(1,:).^0, xi));
    assert(real_check < 1e-10) % Sanity check of kernel decay
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
    assert(fourier_check < 1e-10); % Sanity check Fourier decay
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
    uF
    uS
    uR
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
    global c_pswf
    % Plot things Laplace
    L = 1;
    xi_ewald = 0.8*2*pi/L;
    xi_pswf = L/2 * c_pswf;
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
    rl = L/xi_pswf * c_pswf/2;
    r = [linspace(0, rl, 50)].';
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
    drawnow
end    

function stokes_plots()
    % Plot things Stokes
    global c_pswf % WARNING: GLOBAL VARIABLE
    L = 1;
    xi_ewald = 0.95*2*pi/L;
    xi_pswf = L/2 * c_pswf;
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
    semilogy(kvec, fourier_pswf, '.-', displayname="Stokes Fourier PSWF/erf")
    stem(kmax_pswf, fourier_pswf(find(kvec==kmax_pswf, 1)), displayname="PSWF bandwidth")
    axis([0 kvec(end) eps 1e5])
    xlabel('|k|')
    grid on
    legend
    subplot(1, 2, 2)
    rl = L/xi_pswf * c_pswf/2;
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
    semilogy(r, real_pswf, '.-', displayname="Stokes Real PSWF/erf")
    stem(rl, real_pswf(find(r==rl, 1)), displayname="PSWF cutoff")
    axis([0 L eps inf])
    legend
    grid on
    drawnow
end    
