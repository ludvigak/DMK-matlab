function [Bwin, alpha, beta] = windowed_biharmonic(k, Ctrunc, type)
    arguments
        k
        Ctrunc
        type = 1
    end
    % type:
    %  1 = Tornberg & Bagge. Faster decay, ~1 digit better.
    %  2 = Vico et al.
    ksq = k.^2;
    switch type
      case 1
        Bwin = -8*pi./ksq.^2 .* (1 + 1/2*cos(Ctrunc*k) - 3/2*sin(Ctrunc*k)./(Ctrunc*k));
        Bwin(k==0) = -pi*Ctrunc^4/15;
        alpha = -Ctrunc/2;
        beta = -1/(2*Ctrunc);
      case 2
        Bwin = 4*pi./ksq.^2 .* ((2-Ctrunc^2*k.^2).*cos(Ctrunc*k) + 2*Ctrunc*k.*sin(Ctrunc*k) - 2);
        Bwin(k==0) = -pi*Ctrunc^4;
        alpha = 0;
        beta = 0;
    end
end

