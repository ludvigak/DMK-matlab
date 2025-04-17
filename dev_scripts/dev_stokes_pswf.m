clear
clf

c = 20.736000061035156; % 9 digits in paper
xi = 4;
rl = 1.5;
n = 100;
N = 2*n+1;

r = linspace(-rl, rl, N);




r = linspace(0.1, rl, 50);
[f, df, ddf] = splitfun_pswf(r, rl, c);



h = r(2)-r(1);
plot(r,f, ...
     r(2:end-1), (f(3:end)-f(1:end-2))./(r(3:end)-r(1:end-2)), '.', ...
     r, df, ...
     r(2:end-1), (f(1:end-2)-2*f(2:end-1)+f(3:end))/h^2, '.', ...
     r, ddf)
xlim([0 inf])

return

[A1_pswf, A2_pswf] = stokeslet_funcs(r, f, df, ddf);
plot(r, abs(A1_pswf), displayname='A1_{pswf}'); hold on
plot(r, abs(A2_pswf), displayname='A2_{pswf}')
[f_erfc, df_erfc, ddf_erfc] = splitfun_erfc(r, rl, xi);
[A1_erfc, A2_erfc] = stokeslet_funcs(r, f_erfc, df_erfc, ddf_erfc);
plot(r, abs(A1_erfc), displayname='A1_{erfc}')
plot(r, abs(A2_erfc), displayname='A2_{erfc}')
legend
set(gca, yscale='linear')
xlim([0 rl])
return


oversample = 60;
nfft = oversample*n;
Nfft = 2*nfft+1;
f_hat   = fftshift(fft(f, Nfft));
df_hat   = fftshift(fft(df, Nfft));
ddf_hat   = fftshift(fft(ddf, Nfft));
L = r(end)-r(1)
hf = 2*pi/(L*oversample);
kvec = hf*(-nfft:nfft);


clf
subplot(1, 2, 1)
semilogy(r, f, displayname='\Phi')
hold on
semilogy(r, abs(df), displayname="\Phi'")
semilogy(r, abs(ddf), displayname="\Phi''")
legend
grid on
xlim([0, rl])

subplot(1,2,2)
semilogy(kvec, abs(f_hat), displayname='\Phi')
hold on
semilogy(kvec, abs(df_hat), displayname="\Phi'")
semilogy(kvec, abs(ddf_hat), displayname="\Phi''")
legend
grid on


function [f, df, ddf] = splitfun_erfc(r, rl, xi)
    r = abs(r);
    f = erfc(xi*r/rl);
    df = -(2*xi.*exp(-(r.^2*xi^2)/rl^2))/(rl*pi^(1/2));
    ddf = (4*r*xi^3.*exp(-(r.^2*xi^2)/rl^2))/(rl^3*pi^(1/2));
end

function [f, df, ddf] = splitfun_pswf(r, rl, c)
    psi = pswf(0, c);
    dpsi = diff(psi);
    c0 = integral(psi, 0, 1);
    f = zeros(size(r));
    df = zeros(size(r));
    ddf = zeros(size(r));
    r = abs(r);
    for i=1:numel(f)
        f(i) = integral(psi, r(i)/rl, 1)/c0;
        df(i) = -psi(r(i)/rl)/(c0*rl);
        ddf(i) = -dpsi(r(i)/rl)/(c0*rl^2);
    end
end

function [A1, A2] = stokeslet_funcs(r, f, df, ddf)
    A1 = abs(r).*ddf + df;
    A2 = abs(r).*ddf - df;
end
