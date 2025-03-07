function [f, df, d2f, d3f] = radial_fourier_kernels_gen()
    f   = construct_kernel(0);
    df  = construct_kernel(1);
    d2f = construct_kernel(2);
    d3f = construct_kernel(3);
end

function mfunc = construct_kernel(D)
    [srcpath, ~, ~] = fileparts(mfilename('fullpath'));
    base_name = sprintf('%s/+gen/kernel_d%df', srcpath, D);
    % Re-generate kernel files if they do not exist
    if exist(base_name, 'file') ~= 2
        warning('Generating kernel for D=%d', D)
        taylor_order = 32; % Experimentally verified
        k = sym('k');
        r = sym('r');
        % Base function
        f = k/r*sin(k*r);
        % Differentiate
        dDf = diff(f, D, r);
        series = taylor(dDf, r, 0, Order=taylor_order);
        series2 = taylor(dDf, r, 0, Order=taylor_order+4);
        errest = abs(simplify(series-series2));
        % Construct functions
        matlabFunction(simplify(dDf, 1), Vars=[k,r], file=base_name);
        matlabFunction(limit(dDf, r, 0), Vars=k, file=[base_name '_zerolim']);
        matlabFunction(series, Vars=[k,r], file=[base_name '_series']);
        matlabFunction(errest, Vars=[k,r], file=[base_name '_series_err']);
    end
    func              = str2func(sprintf('gen.kernel_d%df', D));
    lim               = str2func(sprintf('gen.kernel_d%df_zerolim', D));
    dDf_series        = str2func(sprintf('gen.kernel_d%df_series', D));
    dDf_series_errest = str2func(sprintf('gen.kernel_d%df_series_err', D));
    function [y, err] = dDf_func(k,r)
        assert(size(r, 1)==1)
        y = func(k, r);
        err = zeros(size(y));
        mask = r < (10^(D-1) ./ k.^D).^(1/D);
        if any(mask(:))
            S = dDf_series(k, r);
            y(mask) = S(mask);
            E = dDf_series_errest(k, r);
            err(mask) = E(mask);
        end
        if any(r==0)
            y(:,r==0) = lim(k);
        end
    end
    mfunc = @dDf_func;
end
