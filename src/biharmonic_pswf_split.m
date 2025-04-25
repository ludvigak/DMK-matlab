function bsplit = biharmonic_pswf_split(c_pswf)
% Return chebfuns with biharmonic split for given PSWF shape
bsplit = struct();
% Init PSWF
    psi = pswf(0, c_pswf);
    dpsi = pswf0_diff(c_pswf, 1);
    % Fourier space
    f_scaling = psi(0);
    psi_hat = psi / f_scaling;
    dpsi_hat = dpsi / f_scaling;
    k_cheb = chebfun(@(k) k);
    bsplit.psi_hat = psi_hat;
    bsplit.dpsi_hat = dpsi_hat;
    bsplit.gamma_hat = psi_hat  - 1/2*k_cheb*dpsi_hat;
    % Real space
    r_scaling = 1 / sum(psi, -1, 1);
    %lambda = sum(psi, -1, 1) / psi(0);
    %r_scaling = r_scaling * lambda^2*c_pswf/(2*pi); % "correct" transform
    psi = restrict(psi * r_scaling, [0 1]);
    dpsi = restrict(dpsi * r_scaling, [0 1]);
    r = chebfun(@(r) r, [0 1]);
    bsplit.gamma = 3/2*psi + 1/2*r*dpsi;
    bsplit.psi = psi;
    bsplit.dpsi = dpsi;
    bsplit.Bmoll = 2*(r*cumsum(psi) - cumsum(r*psi)) + 2*integral(r*psi);
    bsplit.dBmoll = 2*cumsum(psi);
    bsplit.d2Bmoll = 2*psi;
    bsplit.d3Bmoll = 2*dpsi;
    bsplit.Bres = r  - bsplit.Bmoll;
    bsplit.dBres = 1 - bsplit.dBmoll;
    bsplit.d2Bres = -bsplit.d2Bmoll;
    bsplit.d3Bres = -bsplit.d3Bmoll;
    bsplit.r = r;
end

