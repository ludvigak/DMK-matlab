clear all
global c_pswf
c_pswf = 30


xi = 2;



[psi_hat, dpsi_hat, d2psi_hat, d3psi_hat, d4psi_hat] = pswf_hat(0, 2*xi);
alpha2 = -d2psi_hat / (2*psi_hat)



[psi_hat, dpsi_hat, d2psi_hat, d3psi_hat, d4psi_hat] = pswf_hat(0, xi);
alpha = -d2psi_hat / (2*psi_hat)

limit = (6*d2psi_hat^2-psi_hat*d4psi_hat)/(24*psi_hat^2)

k = linspace(0, 10);
[psi_hat, dpsi_hat, d2psi_hat] = pswf_hat(k, xi);
[psi_hat2, ~, ~] = pswf_hat(k, 2*xi);

[psi_hat_sq] = pswf_hat(k.^2, xi);


semilogy(k, (1 - psi_hat.*(1 + alpha*k.^2)) ./ k.^4, ...
         k, k.^0*limit, '--', ...
         k, (psi_hat.*(1 + alpha*k.^2)) ./ k.^4, ...
         k, exp(-(k/xi).^4) ./ k.^4, ...
         k, psi_hat_sq ./ k.^4)

ylim([1e-15 10])

function [psi_hat, dpsi_hat, d2psi_hat, d3psi_hat, d4psi_hat] = pswf_hat(k, xi)
    % Return scaled function \hat\psi(k*r_l/c) / \psi(0), and derivatives
    global c_pswf % WARNING: GLOBAL VARIABLE
    psi   = pswf(0, c_pswf);
    dpsi  = diff(psi);
    d2psi = diff(dpsi);
    d3psi = diff(d2psi);
    d4psi = diff(d3psi);
    rl = 1/xi * c_pswf/2;
    psi_arg = k*rl/c_pswf;
    [psi_hat, dpsi_hat, d2psi_hat, d3psi_hat, d4psi_hat] = deal(zeros(size(k)));
    mask = (psi_arg <= 1);
    psi_hat(mask)   =   psi(k(mask)*rl/c_pswf)/psi(0);
    dpsi_hat(mask)  =  dpsi(k(mask)*rl/c_pswf)/psi(0) * rl/c_pswf;
    d2psi_hat(mask) = d2psi(k(mask)*rl/c_pswf)/psi(0) * (rl/c_pswf)^2;
    d3psi_hat(mask) = d3psi(k(mask)*rl/c_pswf)/psi(0) * (rl/c_pswf)^3;
    d4psi_hat(mask) = d4psi(k(mask)*rl/c_pswf)/psi(0) * (rl/c_pswf)^4;

    if k==0
        alpha = -d2psi(0)/(2*psi(0)) * (rl/c_pswf)^2
        limit = (6*d2psi(0)^2-psi(0)*d4psi(0))/(24*psi(0)^2) * (rl/c_pswf)^4
    end
end
