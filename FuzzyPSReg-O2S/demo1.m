clc
clear

Po = load('data/Po1.txt')';  % object model (3 x N_Po)
Ps = load('data/Ps1.txt')';  % scene  scan  (3 x N_Ps)

%%%%%% Parameter setting
trimRatio1 = 0.4;  % \xi_1, trimming ratio of Po in the coarse registration
N_Co = 50;         % number of fuzzy clusters of Po
N_Cs = 50;         % number of fuzzy clusters of Ps
gs_Po = 0.005;     % grid step of for the box grid filter to downsample Po for improving efficiency
gs_Ps = 0.005;     % grid step of for the box grid filter to downsample Ps for improving efficiency
p1 = 1.05;         % (1 + dJ_1) for the extended registration quality assessment
p2 = 3;            % (1 + dJ_1) for the extended registration quality assessment
r_inits = [0  0  0;0  pi  0]';  % N_r initial rotations (3 x N_r) of Po in the multi-start local optimization

%%%%%% Fuzzy cluster-based point set registration for object pose estimation
[R,t] = fuzzyPSReg_O2S(Po,Ps,trimRatio1,N_Co,N_Cs,gs_Po,gs_Ps,p1,p2,r_inits);

PoRnT = R*Po + t;  % transformed object model

%%%%%% Plot the results
figure
subplot(121)
hold on
title('Initial poses');
plot3(Ps(1,:),Ps(2,:),Ps(3,:),'.r','MarkerSize',2);
plot3(Po(1,:),Po(2,:),Po(3,:),'.b','MarkerSize',4);

subplot(122)
hold on
title('After registration');
plot3(Ps(1,:),Ps(2,:),Ps(3,:),'.r','MarkerSize',2);
plot3(PoRnT(1,:),PoRnT(2,:),PoRnT(3,:),'.b','MarkerSize',4);

