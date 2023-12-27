function W0hat = laplace_winkernel_fourier(k1, k2, k3, sigma_0, Ctrunc)
% Fourier transform of Laplace windowedo mollified kernel
    ksq = k1.^2 + k2.^2 + k3.^2;
    k = sqrt(ksq);
    W0hat = 8*pi*(sin(Ctrunc*k/2)./k).^2 .* exp(-ksq * sigma_0^2/4);
    W0hat(ksq==0) = 8*pi * Ctrunc^2/4;
end

