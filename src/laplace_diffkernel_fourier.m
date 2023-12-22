function Dlhat = laplace_diffkernel_fourier(k1, k2, k3, sigma_l)
% Fourier transform of Laplace difference kernel
    sigma_lp1 = sigma_l/2; % sigma_{l+1} = sigma_l / 2
    ksq = k1.^2 + k2.^2 + k3.^2;
    Dlhat = 4*pi*(exp(-ksq*sigma_lp1^2/4) - exp(-ksq*sigma_l^2/4))./ksq;
    Dlhat(ksq==0) = pi*(sigma_l^2-sigma_lp1^2);
end

