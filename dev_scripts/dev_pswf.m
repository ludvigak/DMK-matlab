clear

tol = 1e-6;
sigma_0 = 1/sqrt(log(1/tol))
c = log(1/tol)*1.2

sigma_0 = 1/sqrt(log(1e6))
c = 16.893999099731445

nf = 100;
Nf = 2*nf+1;
x = linspace(-3, 3, Nf);

psi = pswf(0, c);

c_0 = integral(psi, 0, 1, abstol=1e-13)/2;

psi_x = psi(x)/psi(0);

gauss_x = exp(-x.^2/sigma_0.^2);
psi_x(abs(x)>1)=0;
psi_trunc = psi_x(end)
gauss_trunc = gauss_x(end)

phi_diff = -psiintdiff(x, c) ./ abs(x);

clf
subplot(1, 2, 1)
plot(x, psi_x, displayname='\psi')
hold on
plot(x, gauss_x, displayname='erfc')

plot(x, phi_diff, displayname='\Phi')

legend
set(gca, 'yscale', 'linear')
xlim([0, 1])
grid on

subplot(1, 2, 2)

oversample = 60;
nfft = oversample*nf;
Nfft = 2*nfft+1;
psi_hat   = fftshift(fft(psi_x, Nfft));
gauss_hat = fftshift(fft(gauss_x, Nfft));

L = x(end)-x(1)
hf = 2*pi/(L*oversample);
kvec = hf*(-nfft:nfft);


semilogy(kvec, abs(psi_hat))
hold on
semilogy(x*c, psi_x*max(abs(psi_hat)), '.')

semilogy(kvec, abs(gauss_hat))
stem(c, 1)
xlim([0 2*c])
grid on

function Psi = psiint(r, c)
    psi = pswf(0, c);
    c_0 = integral(psi, 0, 1, abstol=1e-13)/2;
    Psi = zeros(size(r));
    for i=1:numel(r)
        arg = abs(r(i));
        if arg > 1
            Psi(i) = 1;
        elseif arg == 0
            Psi(i) = 0;
        else
            Psi(i) = integral(psi, 0, abs(arg), abstol=1e-13)/c_0;
        end
    end
end

function Psi = psiintdiff(r, c)
    psi = pswf(0, c);
    c_0 = integral(psi, 0, 1)/2;
    Psi = zeros(size(r));
    for i=1:numel(r)
        arg = abs(r(i));
        if arg > 1
            Psi(i) = 1;
        elseif arg == 0
            Psi(i) = 0;
        else
            Psi(i) = integral(psi, 0, abs(arg)/2)/c_0 - ...
                     integral(psi, 0, abs(arg))/c_0;
        end
    end
end
