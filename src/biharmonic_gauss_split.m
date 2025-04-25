function bsplit = biharmonic_gauss_split(sigma, Kmax)
% Return chebfuns with biharmonic split for given Gauss shape

    bsplit = struct();

    % Fourier space
    psi_hat = chebfun(@(k) exp(-k.^2*sigma^2/4), [0 Kmax]);
    dpsi_hat = chebfun(@(k) -2*k*sigma^2/4.*exp(-k.^2*sigma^2/4), [0 Kmax]);
    k_cheb = chebfun(@(k) k, [0 Kmax]);
    bsplit.psi_hat = psi_hat;
    bsplit.dpsi_hat = psi_hat;
    bsplit.gamma_hat = psi_hat  - 1/2*k_cheb*dpsi_hat;
    % Real space
    psi = chebfun(@(r) exp(-r.^2/sigma^2) ./ (sigma*sqrt(pi)), [0 2]);
    dpsi = chebfun(@(r) -2*r/sigma^2.*exp(-r.^2/sigma^2) ./ (sigma*sqrt(pi)), [0 2]);
    r = chebfun(@(r) r, [0 2]);
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

