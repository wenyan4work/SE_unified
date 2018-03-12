rng(1)
box = [1 1 1];   % domain
N = 100;         % number of charged particles

M0 = 28;

opt.M = M0*box;
opt.xi = pi*M0 / 12;
opt.rc = 6 / opt.xi;
opt.box = box;

opt.layers = (opt.M-1)/2;

opt.beta = 2.3;

% charge-neutral system
[x q] = SE_charged_system(N,box,'scalar');

idx = 1:10;
% compute FD Ewald sum
ref = SE3P_direct_fd_mex(idx,x,q,opt);

opt.P = 16;
[u time]= se3p_fourier_space_kaiser(1:N,x,q,opt);

e_rms   = rms(u(idx)-ref)/rms(ref)
