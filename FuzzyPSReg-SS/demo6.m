clc
clear

P_F = load('data/fixedSet6.txt')';   % fixed  point set P_F (3 x N_PF)
P_M = load('data/movingSet6.txt')';  % moving point set P_M (3 x N_PM)

%%%%%% Initial pose setting for P_M
r0 = -pi + 2*pi*rand(1,3);
theta = norm(r0);
P_M = axang2rotm([r0/theta, theta])*P_M;

%%%%%% Parameter setting
noisyF = 0;          % 1/0 for the condition of P_F with/without noisy outliers
noisyM = 0;          % 1/0 for the condition of P_M with/without noisy outliers
trimRatio = 0;       % trimming ratio for partially overlapping point sets
N_C = 30;            % number of fuzzy clusters of each point set in the coarse registration
gsF = 0;             % grid step of the box grid filter to dowmsample P_F for the fuzzy clustering
gsM = 0;             % grid step of the box grid filter to dowmsample P_M for the fuzzy clustering
gsF_fine = 0;        % grid step of the box grid filter to dowmsample P_F for the fine registration
gsM_fine = 0;        % grid step of the box grid filter to dowmsample P_M for the fine registration
RRange = pi;         % half side-length of the initial rotation    cube [-RRange, RRange]^3 (radians)
TRange = 0.5;        % half side-length of the initial translation cube [-TRange, TRange]^3
sigma_r = 0.01;      % half side-length of the minimum rotation    cube for stopping the global search
sigma_t = TRange/4;  % half side-length of the minimum translation cube for stopping the inner  search
epsilon = 0;         % threshold value for stopping the global search/inner search

%%%%%% Fuzzy cluster-based point set registration for similar-sized scans
[R,t] = fuzzyPSReg_SS(P_F,P_M,...
    noisyF,noisyM,trimRatio,N_C,gsF,gsM,gsF_fine,gsM_fine,RRange,TRange,sigma_r,sigma_t,epsilon);

P_MRnT = R*P_M + t;  % transformed moving set

%%%%%% Plot the results
figure
subplot(121)
hold on
title('Initial poses');
plot3(P_F(1,:),P_F(3,:),P_F(2,:),'.r','MarkerSize',2);
plot3(P_M(1,:),P_M(3,:),P_M(2,:),'.b','MarkerSize',2);

subplot(122)
hold on
title('After registration');
plot3(P_F(1,:),P_F(3,:),P_F(2,:),'.r','MarkerSize',2);
plot3(P_MRnT(1,:),P_MRnT(3,:),P_MRnT(2,:),'.b','MarkerSize',2);




