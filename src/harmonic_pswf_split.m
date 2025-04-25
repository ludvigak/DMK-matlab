function hsplit = harmonic_pswf_split(c_pswf)
% Return chebfuns with harmonic split for given PSWF shape
    hsplit = struct();
    % Init PSWF
    psi = pswf(0, c_pswf);
    dpsi = pswf0_diff(c_pswf, 1);
    d2psi = pswf0_diff(c_pswf, 2);
    % Fourier space
    f_scaling = psi(0);
    psi_hat = psi / f_scaling;
    d2psi_hat = d2psi / f_scaling;
    hsplit.gamma_hat = psi_hat;
    hsplit.d2gamma_hat = d2psi_hat;
    % Real space
    r_scaling = sum(psi, -1, 1);
    psi = restrict(psi / r_scaling, [0 1]);
    dpsi = restrict(dpsi / r_scaling, [0 1]);
    r = chebfun(@(r) r, [0 1]);
    hsplit.Phi = 1-2*cumsum(psi);
    hsplit.dPhi = -2*psi;
    hsplit.d2Phi = -2*dpsi;
    hsplit.r = r;
    hsplit.psi = psi;
end

