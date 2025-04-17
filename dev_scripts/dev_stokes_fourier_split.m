clear all
rng(1);
tol = 1e-14;

disp("STOKESLET")
kernel = kernels.stokeslet_hasimoto(tolerance=tol);

% Setup root level windowed kernel
Kmax_win = ceil( kernel.Kmax );
Ctrunc = sqrt(3) + 1;

nf_win = Kmax_win;
hf_win = Kmax_win/nf_win;


k = hf_win*(-nf_win:nf_win);
[k1, k2, k3] = ndgrid(k);
k1 = k1(:); k2 = k2(:); k3 = k3(:);
Nf = numel(k1);
w = hf_win^3 / (2*pi)^3;

src = [0 0 0];


% Stokeslet
trg = [0.25 0.2 0.3];
f = [1 1 1];
fhat = ones(Nf, 3);
W0hat_op = kernel.winkernel_fourier(k1, k2, k3, Ctrunc);
[W0hat const] = W0hat_op(fhat);
SR = kernel.reskernel(trg, src, f, 0);
SW = real(sum(w .* W0hat .* exp(-1i*(k1*trg(1) + k2*trg(2) + k3*trg(3)))));

S = kernel.direct(trg, src, f);
err=SR+SW+const- S

disp("STRESSLET")


p = 40;
N = 1;
tol = 1e-12;


kernel = kernels.stresslet_hasimoto(tolerance=tol);

% Setup root level windowed kernel
Kmax_win = ceil( kernel.Kmax );
Ctrunc = sqrt(3) + 1;

nf_win = Kmax_win;
hf_win = Kmax_win/nf_win;


k = hf_win*(-nf_win:nf_win);
[k1, k2, k3] = ndgrid(k);
k1 = k1(:); k2 = k2(:); k3 = k3(:);
Nf = numel(k1);
w = hf_win^3 / (2*pi)^3;


src = [pi exp(1) sqrt(2)]/10;
trg = [1 3 5]/10;
f = [1 0 0 1 0 0];
%f = rand(1,6);
q = f(1:3);
n = f(4:6);


tree = octree([0 0 0], 1);
proxy_charges = init_proxy_charges(tree, [0 0 0 0 0 0], p);
[rvec, V] = approx.chebvander(p);
box_proxy_points = tree.box_grid(1, rvec);
box_proxy_charges = proxy_charges{1};

box_proxy_charges(:) = 0;
box_proxy_charges(end, :) = f;

trg = box_proxy_points(1,:)
src = box_proxy_points(end, :)


Twin = operator_windowed(p, hf_win, nf_win, Ctrunc, kernel);
[far_expa, uhat] = Twin(box_proxy_charges);

% 
fhat = ones(Nf, 3, 3);
for i1=1:3
    for i2=1:3
        fhat(:,i1,i2) = exp(-1i*(k1*src(1) + k2*src(2) + k3*src(3))) .* q(:, i1) .* n(:, i2);
    end
end
%fhat(:, 1, 1) = 1;

W0hat_op = kernel.winkernel_fourier(k1, k2, k3, Ctrunc);
[W0hat const] = W0hat_op(fhat);
%uRref = stresslet_real_space([trg;src], [f;f], 1/kernel.sigma_0, 100)
TR = kernel.reskernel(trg, src, f, 0)


my_uhat = w .* W0hat;
TW = real(sum(my_uhat .* exp(1i*(k1*trg(1) + k2*trg(2) + k3*trg(3))) ))
const
T = kernel.direct(trg, src, f)
Ttot = TR+TW+const

err = Ttot - T


function [u varargout] = stresslet_real_space(x, f, xi, rc)
    q = f(:,[1 2 3]);
    n = f(:,[4 5 6]);

    N = size(x, 1);
    timestamp = tic();
    [idx, d] = rangesearch(x, x, rc);
    walltime.nblist = toc(timestamp);

    timestamp = tic();
    MATLAB = true;
    if MATLAB
        u = zeros(N, 3);
        for target = 1:N
            nnb_idx = idx{target}(2:end);
            rvec = bsxfun(@minus, x(target,:), x(nnb_idx,:));
            r = d{target}(2:end)';
            qsrc = q(nnb_idx,:);
            nsrc = n(nnb_idx,:);
            
            r2 = r.^2;
            c = xi^2*r2;    
            % Beenakker
            %C=-2./r2.^2.*( 3.0./r.*erfc(xi*r) + 2.0*xi/sqrt(pi)*(3.0+2.0*c-4.0*c.^2).*exp(-c) );
            %D=8/sqrt(pi)*xi^3*(2.0-c).*exp(-c);
            % Hasimoto
            C=-2./r2.^2.*( 3.0./r.*erfc(xi*r) + 2.0*xi/sqrt(pi)*(3.0+2.0*c).*exp(-c) );
            D=4/sqrt(pi)*xi^3.*exp(-c);
            
            rdotn = sum(rvec .* nsrc, 2);
            rdotq = sum(rvec .* qsrc, 2);
            ndotq = sum(nsrc .* qsrc, 2);
            
            u(target,:) = sum(...
                bsxfun(@times, C.*rdotn.*rdotq + D.*ndotq, rvec) + ...
                bsxfun(@times, D.*rdotq, nsrc) +...
                bsxfun(@times, D.*rdotn, qsrc) ...
                , 1);
        end
    else
        u = stresslet_rsrc_mex(x,f,idx,d,xi);
    end
    walltime.eval = toc(timestamp);

    if nargout==2
        walltime.total = walltime.nblist + walltime.eval;
        varargout{1} = walltime;
    end

end

