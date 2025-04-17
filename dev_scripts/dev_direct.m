clear;
rng(1)


kernel = kernels.stokeslet_hasimoto();
N = 2000;

points = rand(N, 3)-1/2;
charges = rand(N, kernel().dim_in)-1/2;

u_ref = kernel.direct(points, points, charges); % Warmup

atic = tic();
M = 4;
for i=1:4
    u_ref = kernel.direct(points, points, charges);
end
toc(atic)
