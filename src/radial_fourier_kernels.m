function [f, df, d2f, d3f] = radial_fourier_kernels()
% Return function handles for the kernel of the radial Fourier transform
% f(k, r) = k/r*sin(k*r)
% and derivatives w.r.t. radius r.
% === Generating code:
% syms k r; f = k/r*sin(k*r);
% for d=0:3
%     fprintf('d%df     = %s\n', d, func2str(matlabFunction(simplify(diff(f,d)))))
%     fprintf('d%df_lim = %s\n', d, func2str(matlabFunction(limit(diff(f,d)), r, 0)))
% end
% ==================
    f = @(k,r)  k./r   .*(sin(k.*r));
    function y=df_(k,r)
        assert(size(r, 1)==1)
        y = -k./r.^2.*(sin(k.*r)-k.*r.*cos(k.*r));
        if any(r==0)
            y(:,r==0) = 0;
        end
    end
    function y=d2f_(k,r)
        assert(size(r, 1)==1)
        y = -k./r.^3.*(sin(k.*r).*-2.0+k.^2.*r.^2.*sin(k.*r)+k.*r.*cos(k.*r).*2.0);
        if any(r==0)
            y(:, r==0) = -k.^4/3;
        end
    end
    function y=d3f_(k,r)
        assert(size(r, 1)==1)
        y = -k.*1.0./r.^4.*(sin(k.*r).*6.0+k.^3.*r.^3.*cos(k.*r)-k.^2.*r.^2.*sin(k.*r).*3.0-k.*r.*cos(k.*r).*6.0);
        if any(r==0)
            y(:, r==0) = 0;
        end
    end
    df = @df_;
    d2f = @d2f_;
    d3f = @d3f_;
end
