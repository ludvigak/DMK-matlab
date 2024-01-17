classdef SplitKernelInterface < handle
% Abstract class for interfacing to split kernels
    properties
        tolerance
        % Tolerance used in kernel
        Kmax
        % Maximum Fourier mode needed to resolve kernel
    end
    methods
        function obj = SplitKernelInterface(args)
        % Constructur, takes keyword argument tolerance
            arguments
                args.tolerance (1,1) double = 0 % tolerance default=0 (unset)
            end
            obj = obj@handle();
            obj.tolerance = args.tolerance;
        end
    end
        
    properties (Abstract, Constant)
        name
        % Name of the split kernel
        dim_in
        % Dimension of kernel input data, eg. 1 for Laplace, 3 for Stokeslet
        dim_out
        % Dimension of kernel output data, eg. 1 for Laplace, 3 for Stokeslet
    end
    methods (Abstract)
        u = direct(targets, sources, data)
        % Direct evaluation of the kernel
        % targets: (Ntrg, 3)
        % sources: (Nsrc, 3)
        % data: (Nsrc, dim_in)
        % returns: (Ntrg, dim_out)

        ures = reskernel(targets, points, charges, level)
        % Residual kernel R_l(r)
        
        udiff = diffkernel(targets, points, charges, level)
        % Difference kernel D_l(r)

        Dlhat_fun = diffkernel_fourier(k1, k2, k3, level)
        % Fourier transform o difference kernel D_l(r)
        % returns: function handle that applies kernel
        
        W0hat_fun = winkernel_fourier(k1, k2, k3, Ctrunc)
        % Fourier transform of windowed mollified kernel at 0 level
        % returns: function handle that applies kernel

        uself = self_interaction(charges, level)
        % Self interaction
    end
end
