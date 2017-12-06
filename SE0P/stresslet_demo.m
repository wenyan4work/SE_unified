clear
rng(1);

L = 1;
box = [L L L];
N = 1000;
[x, f] = vector_system(N, box, 6);

% Maxmimize point-to-point distance
x(1,:) = 0;
x(2,:) = box;

M0 = 30; % Set M0, the rest is auto

opt.M = M0*box;
opt.xi = pi*M0 / 13;
opt.P = 26;
opt.oversampling = 1+sqrt(3);
opt.rc = 6 / opt.xi;
opt.box = box;
opt

% Precomp
tic
pre = stresslet_precomp(opt);
toc
% Ewald
tic
uf = stresslet_fourier_space(x, f, opt, pre);
ur = stresslet_real_space(x, f, opt);
ue = uf+ur;
toc
% > 1.3

% Direct
tic
u = stresslet_direct(x, f, box);
toc

error = norm(u-ue, inf) / norm(u,inf)