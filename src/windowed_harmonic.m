function Hwin = windowed_harmonic(k, Ctrunc)
    Hwin = 8*pi*(sin(Ctrunc*k/2)./k).^2;
    Hwin(k==0) = 8*pi * Ctrunc^2/4;
end

