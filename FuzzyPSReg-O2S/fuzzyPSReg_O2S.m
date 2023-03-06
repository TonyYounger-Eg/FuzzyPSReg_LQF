%%%%%% The main function
function [R,t,timeReg] = fuzzyPSReg_O2S(Po,Ps,trimRatio1,N_Co,N_Cs,gs_Po,gs_Ps,p1,p2,r_inits)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FuzzyPSReg-O2S: Fuzzy cluster-based point set registration for registering an object to a scene.
%   This function estimates the rigid transformation parameters, R and t, to register the two given 3D
%   point sets Po and Ps, where Po (3 x N_Po) is the object model and Ps (3 x N_Ps) is the scene scan.
%
%   Inputs:
%     'Po'                   object model (3 x N_Po).
%
%     'Ps'                   scene  scan  (3 x N_Ps).
%
%     'trimRatio1'           trimming ratio \xi_1 of Po in the multi-start optimization
%                            of the coarse registration. trimRatio1 = 0 for no trimming case.
%
%     'N_Co'                 number of fuzzy clusters of Po.
%
%     'N_Cs'                 number of fuzzy clusters of Ps.
%
%     'gs_Po'                grid step of the box grid filter to downsample Po for improving efficiency,
%                            gs_Po = 0 for no downsampling.
%
%     'gs_Ps'                grid step of the box grid filter to downsample Ps for improving efficiency,
%                            gs_Ps = 0 for no downsampling.
%
%     'p1'                   (1 + dJ_1) for the extended registration quality assessment.
%
%     'p2'                   (1 + dJ_2) for the extended registration quality assessment.
%
%     'r_inits'              N_r initial rotations (3 x N_r) of Po in the multi-start local optimization
%                            of the fine registration.
%
%   Outputs:
%     'R'                    3 x 3 rotation    matrix to transform Po for registration.
%
%     't'                    3 x 1 translation vector to transform Po for registration.
%
%     'timeReg'              time costs (sec) of the registration.
%
%   References:
%     [1] Q. Liao, D. Sun, and H. Andreasson, "FuzzyPSReg: Strategies of Fuzzy Cluster-based
%                  Point Set Registration", IEEE Transactions on Robotics, 2021.
%
%     [2] Q. Liao, D. Sun, and H. Andreasson, "Point Set Registration for 3D Range Scans
%                  Using Fuzzy Cluster-based Metric and Efficient Global Optimization",
%                  IEEE Transactions on Pattern Analysis and Machine Intelligence, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[D_Po,N_Po] = size(Po);
[D_Ps,N_Ps] = size(Ps);
[D_r,N_r] = size(r_inits);

if ( N_Po == 0 )||( N_Ps == 0 )||( D_Po ~= 3 )||( D_Ps ~= 3 )
    error('Something wrong with the input point sets. ');
end
if ( N_r == 0 )||( D_r ~= 3 )
    error('Incorrect setting for r_inits. ');
end

objMax = max(Po,[],2);  % users can choose a proper value for objMax
objMin = min(Po,[],2);  % users can choose a proper value for objMin
objMax2Min = objMax - objMin;
objMax_s = objMax - 0.1*objMax2Min;
objMin_s = objMin + 0.1*objMax2Min;
objMax_l = objMax + 0.1*objMax2Min;
objMin_l = objMin - 0.1*objMax2Min;

Po = downsamplePoints(Po,gs_Po);
[Co_fcm,AFPCD_fcm,U] = fuzzyClusteringFCM(Po,N_Co,100);
[Co_gk,AFPCD_gk,pInvCov_gk] = fuzzyClusteringGK(Po,N_Co,30,U);

P_F1 = Ps;
N_PF1 = N_Ps;
gs_PF1 = gs_Ps;
N_CF1 = N_Cs;
C_M1 = Co_fcm;
N_CM1 = N_Co;
N_CM1Trim = round( N_CM1*(1-trimRatio1) );
C_F2 = Co_gk;
N_CF2 = N_Co;
gs_PM2 = gs_Ps;
trimRatio2 = 0;  % trimming ratio \xi_2

timeReg = 0;
N_PsMin = 3000;  % threshold for the number of points in Ps
zeros6xNr = zeros(6,N_r);
largeNum1xNr = 10^6*ones(1,N_r);
options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','off');

fprintf('Registration starts. \n\n')

while (1)
    
    if N_PF1 < N_PsMin
        fprintf('   The number of points in the scene scan is below the threshold. \n')
        fprintf('   Rescan the scene. \n\n')
        P_F1 = Ps;
        N_PF1 = size(P_F1,2);
        continue
    end
    
    fprintf('A coarse registration starts. \n')
    tic
    [R1,t1,P_F1,N_PF1,flag_ms] = coarseReg(P_F1,N_CF1,C_M1,N_CM1,N_CM1Trim,gs_PF1,...
                                           N_PsMin,objMax,objMin,objMax_s,objMin_s);
    timeReg = timeReg + toc;
    fprintf('The coarse registration is complete. \n')
    
    if flag_ms == 0
        continue
    end
    % R1
    % t1
    
    ones1xNPF1 = ones(1,N_PF1);
    inv_P_F1 = R1'*(P_F1 - t1*ones1xNPF1);
    idxs_PM2 = ( (inv_P_F1(1,:)<objMax(1))&(inv_P_F1(1,:)>objMin(1))...
                &(inv_P_F1(2,:)<objMax(2))&(inv_P_F1(2,:)>objMin(2))...
                &(inv_P_F1(3,:)<objMax(3))&(inv_P_F1(3,:)>objMin(3)) );
    N_PM2 = sum(idxs_PM2);
    
    if N_PM2 <= 1000
        fprintf('   The result is not good, prune the scan and rerun the coarse registration. \n\n')
        P_F1(:,idxs_PM2) = [];
        N_PF1 = size(P_F1,2);
        continue
    else
        P_M2 = inv_P_F1(:,idxs_PM2);
        P_M2 = downsamplePoints(P_M2,gs_PM2);
        N_PM2 = size(P_M2,2);
    end
    N_PM2Trim = round( N_PM2*(1-trimRatio2) );
    
    AFPCD2xN_PM2Trim = AFPCD_gk*N_PM2Trim;
    p1AFPCD2xN_PM2Trim = p1*AFPCD2xN_PM2Trim;
    p2AFPCD2xN_PM2Trim = p2*AFPCD2xN_PM2Trim;
    
    fprintf('A fine registration starts. \n')
    rt2s = zeros6xNr;
    rt2Costs = largeNum1xNr;
    
    tic
    for i = 1:N_r
        rt2_init = [r_inits(:,i)'  0  0  0]';
        [rt2s(:,i), rt2Costs(i)]= fminunc(@(rt)fuzzyMetricGK(rt,C_F2,P_M2,...
            N_CF2,N_PM2,N_PM2Trim,pInvCov_gk),rt2_init,options);
        if rt2Costs(i) <= p1AFPCD2xN_PM2Trim
            break
        end
    end
    timeReg = timeReg + toc;
    fprintf('The fine registration is complete. \n\n')
    
    [rt2Cost,idx_rt2] = min(rt2Costs);
    rt2 = rt2s(:,idx_rt2);
    r2 = rt2(1:3);
    theta = norm(r2);
    if theta == 0
        R2 = [1 0 0;0 1 0;0 0 1];
    else
        R2 = axang2rotm([r2'/theta, theta]);
    end
    t2 = rt2(4:6);
    % R2
    % t2
    
    R = R1*R2';
    t = t1 - R*t2;
    
    if rt2Cost > p2AFPCD2xN_PM2Trim
        
        fprintf('q_gk = -1, the result is wrong. Thus, the scan is pruned for a new round. \n\n')
        %%%%%%
        % P_F1(:,idxs_PM2) = [];
        % N_PF1 = size(P_F1,2);
        %%%%%%
        inv_P_F1 = R'*(P_F1 - t*ones1xNPF1);
        idxs_PM2 = ( (inv_P_F1(1,:)<objMax_s(1))&(inv_P_F1(1,:)>objMin_s(1))...
                    &(inv_P_F1(2,:)<objMax_s(2))&(inv_P_F1(2,:)>objMin_s(2))...
                    &(inv_P_F1(3,:)<objMax_s(3))&(inv_P_F1(3,:)>objMin_s(3)));
        P_F1(:,idxs_PM2) = [];
        N_PF1 = size(P_F1,2);
        
    elseif rt2Cost > p1AFPCD2xN_PM2Trim
        
        fprintf('q_gk = 0, please check the result. \n\n')
        plot4checking(P_F1,Po,R,t)
        
        choice = input('Is the result OK (at least an acceptable position)? (1: Yes, other numbers: No): ');
        if isempty(choice)
            choice = 1;
        end
        
        if choice == 1
            
            fprintf('The result needs further correction. FuzzyPSReg-SS is applied. \n')
            close
            inv_P_F1 = R'*(P_F1 - t*ones1xNPF1);
            idxs_PM2 = ( (inv_P_F1(1,:)<objMax_l(1))&(inv_P_F1(1,:)>objMin_l(1))...
                        &(inv_P_F1(2,:)<objMax_l(2))&(inv_P_F1(2,:)>objMin_l(2))...
                        &(inv_P_F1(3,:)<objMax_l(3))&(inv_P_F1(3,:)>objMin_l(3)));
            P_M2 = inv_P_F1(:,idxs_PM2);
            
            plot4checking(Po,P_M2)
            while (1)
                trimRatio3 = input('Please set a proper trimming ratio for FuzzyPSReg-SS: ');
                if ( trimRatio3 < 0 )||( trimRatio3 >= 1 )  % trimming ratio \xi_3
                    fprintf('Incorrect setting for the trimming ratio. \n\n')
                else
                    fprintf('\n')
                    break
                end
            end
            close
            
            fprintf('A registration using FuzzyPSReg-SS starts. \n')
            tic
            [R3,t3] = fuzzyPSReg_SS(Co_fcm,AFPCD_fcm,Co_gk,pInvCov_gk,N_Co,P_M2,trimRatio3,gs_PM2);
            timeReg = timeReg + toc;
            fprintf('The registration using FuzzyPSReg-SS is complete. \n\n')
            % R3
            % t3
            
            R = R*R3';
            t = t - R*t3;
            
            plot4checking(Po,P_M2,R3,t3)
            % R
            % t
            break
            
        else
            
            fprintf('The result is wrong. Thus, the scan is pruned for a new round. \n\n')
            close
            %%%%%%
            % P_F1(:,idxs_PM2) = [];
            % N_PF1 = size(P_F1,2);
            %%%%%%
            inv_P_F1 = R'*(P_F1 - t*ones1xNPF1);
            idxs_PM2 = ( (inv_P_F1(1,:)<objMax_s(1))&(inv_P_F1(1,:)>objMin_s(1))...
                        &(inv_P_F1(2,:)<objMax_s(2))&(inv_P_F1(2,:)>objMin_s(2))...
                        &(inv_P_F1(3,:)<objMax_s(3))&(inv_P_F1(3,:)>objMin_s(3)));
            P_F1(:,idxs_PM2) = [];
            N_PF1 = size(P_F1,2);
            
        end
        
    else
        
        fprintf('q_gk = 1, the result is correct. \n\n')
        break
        
    end
    
end

fprintf('Registration is complete, please see the following results. \n');
fprintf('The rotation matrix is: \n');
R
fprintf('The translation vector is: \n');
t
fprintf('The time cost (sec, excluding the checking time of the user) is: \n');
timeReg

end


%%%%%% The other functions
function Pds = downsamplePoints(P,gs)  % point set downsampling using box grid filter

if gs == 0  % no downsampling is applied when gs is 0
    
    Pds = P;
    
else
    
    ptCloud = pointCloud(P');
    ptCloud = pcdownsample(ptCloud,'gridAverage',gs);
    Pds = double(ptCloud.Location');  % 3 x N_Pds
    
end

end


function [C,AFPCD,U] = fuzzyClusteringFCM(P,N_C,maxIter)  % FCM clustering

[D_P,N_P] = size(P);

if ( N_C < 2 )||( N_C > N_P )
    error('Incorrect setting for the number of fuzzy clusters. ');
end

onesDPx1 = ones(D_P,1);
ones1xNP = ones(1,N_P);
onesNCx1 = ones(N_C,1);
minDiff = 0;  % minimum difference

U = rand(N_C,N_P);
U = U./(onesNCx1*sum(U,1));  % random initialization for fuzzy membership matrix
U2 = U.^2;  % m = 2
rowSumU2 = sum(U2,2);  % N_C x 1
C = P*U2'./(onesDPx1*rowSumU2');  % fuzzy cluster centers

if N_P < 20000
    e2 = @(a,b) sum(bsxfun(@minus,permute(a,[3 2 1]),permute(b,[2 3 1])).^2,3);
    
    dist2 = e2(P,C);  % N_C x N_P, squared distance matrix
    obj = sum(U2.*dist2,'all');  % FCM clustering objective function (distance loss)
    maxIter = maxIter - 1;
    
    for kIter = 1:maxIter
        
        tmp = 1./dist2;
        U = tmp./(onesNCx1*sum(tmp,1));
        U2 = U.^2;
        rowSumU2 = sum(U2,2);
        C = P*U2'./(onesDPx1*rowSumU2');  % updated fuzzy cluster centers
        dist2 = e2(P,C);
        obj_1 = obj;
        obj = sum(U2.*dist2,'all');
        
        if abs(obj - obj_1) <= minDiff
            break;
        end
        
    end
else
    dist2 = 0*U;  % N_C x N_P, squared distance matrix
    for i = 1:N_C
        dist2(i,:) = sum((C(:,i)*ones1xNP - P).^2,1);
    end
    obj = sum(U2.*dist2,'all');  % FCM clustering objective function (distance loss)
    maxIter = maxIter - 1;
    
    for kIter = 1:maxIter
        
        tmp = 1./dist2;
        U = tmp./(onesNCx1*sum(tmp,1));
        U2 = U.^2;
        rowSumU2 = sum(U2,2);
        C = P*U2'./(onesDPx1*rowSumU2');  % updated fuzzy cluster centers
        
        for i = 1:N_C
            dist2(i,:) = sum((C(:,i)*ones1xNP - P).^2,1);
        end
        
        obj_1 = obj;
        obj = sum(U2.*dist2,'all');
        
        if abs(obj - obj_1) <= minDiff
            break;
        end
        
    end
end

AFPCD = obj/N_P;

end


function [C,AFPCD,pInvCov] = fuzzyClusteringGK(P,N_C,maxIter,U)  % GK clustering

[D_P,N_P] = size(P);

onesDPx1 = ones(D_P,1);
ones1xNP = ones(1,N_P);
onesNCx1 = ones(N_C,1);
minDiff = 0;  % minimum difference
obj = 0;  % GK clustering objective function (distance loss)

if nargin < 4  % if no initialization is provided, random initialization is used
    U = rand(N_C,N_P);
    U = U./(onesNCx1*sum(U,1));
end

MD2 = 0*U;  % N_C x N_P, squared Mahalanobis distance matrix, dist*(p*inv(Cov))*dist
maxIter = maxIter - 1;
a = 1/D_P;

for kIter = 1:maxIter
    
    U2 = U.^2;  % fuzzy membership matrix squared (m = 2)
    rowSumU2 = sum(U2,2);  % N_C x 1
    C = P*U2'./(onesDPx1*rowSumU2');  % updated fuzzy cluster centers
    
    for i = 1:N_C
        dist = P - C(:,i)*ones1xNP;  % D_P x N_P
        U2D = (onesDPx1*U2(i,:)).*dist;
        cov_i = dist*U2D'/rowSumU2(i);  % D_P x D_P, fuzzy covariance matrix
        pInvCov_i = (det(cov_i)^a)*cov_i^(-1);  % D_P x D_P, norm-inducing matrix
        MD2(i,:) = sum(dist.*(pInvCov_i*dist),1);
    end
    
    obj_1 = obj;
    obj = sum(U2.*MD2,'all');
    
    if abs(obj - obj_1) <= minDiff
        break;
    end
    
    tmp = 1./MD2;
    U = tmp./(onesNCx1*sum(tmp,1));
    
end

U2 = U.^2;
rowSumU2 = sum(U2,2);  % N_C x 1
C = P*U2'./(onesDPx1*rowSumU2');  % updated fuzzy cluster centers
pInvCov = zeros(D_P,D_P,N_C);

for i = 1:N_C
    
    dist = P - C(:,i)*ones1xNP;
    U2D = (onesDPx1*U2(i,:)).*dist;
    cov_i = dist*U2D'/rowSumU2(i);
    pInvCov_i = (det(cov_i)^a)*cov_i^(-1);
    MD2(i,:) = sum(dist.*(pInvCov_i*dist),1);
    pInvCov(:,:,i) = pInvCov_i;
    
end

obj = sum(U2.*MD2,'all');
AFPCD = obj/N_P;

end


function [R,t,P_F,N_PF,flag_ms] = coarseReg(P_F,N_CF,C_M,N_CM,N_CMTrim,gs_PF,...
    N_PsMin,objMax,objMin,objMax0,objMin0)  % coarse registration using multi-start optimization

N = 100;
Nf = 10;
flag_ms = 1;
rtCosts = zeros(1,N);
rtfs = zeros(6,Nf);
rtfCosts = rtfs(1,:);
R = [1 0 0;0 1 0;0 0 1];
t = [0 0 0]';
options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','off');

while (1)
    
    N_PF = size(P_F,2);
    if N_PF < N_PsMin
        flag_ms = 0;
        break
    end
    
    Pds_F = downsamplePoints(P_F,gs_PF);
    [C_F,AFPCD,~] = fuzzyClusteringFCM(Pds_F,N_CF,100);
    AFPCDxN_CMTrim = AFPCD*N_CMTrim;
    
    PFMax = max(P_F,[],2) + 0.1;
    PFMin = min(P_F,[],2) - 0.1;
    rtMax = [ pi  pi  pi  (PFMax-objMax)']';
    rtMin = [-pi -pi -pi  (PFMin-objMin)']';
    rts = rtMin + rand(6,N).*(rtMax - rtMin);
    for i = 1:N
        rtCosts(i) = fuzzyMetricFCM(rts(:,i),C_F,C_M,N_CF,N_CM,N_CMTrim);
    end
    for i = 1:Nf
        [~,idx_minCost] = min(rtCosts);
        rtCosts(idx_minCost) = inf;
        [rtfs(:,i), rtfCosts(i)]= fminunc(@(rt)fuzzyMetricFCM(rt,C_F,C_M,...
            N_CF,N_CM,N_CMTrim),rts(:,idx_minCost),options);
    end
    [rtCost,idx_minCost] = min(rtfCosts);
    rt = rtfs(:,idx_minCost);
    
    r = rt(1:3)';
    theta = norm(r);
    if theta == 0
        R = [1 0 0;0 1 0;0 0 1];
    else
        R = axang2rotm([r/theta  theta]);
    end
    t = rt(4:6);
    
    if rtCost < AFPCDxN_CMTrim
        break
    else
        fprintf('   Prune the scan to rerun the multi-start optimization. \n')
        inv_P_F = R'*(P_F - t*ones(1,N_PF));
        idxs = ( (inv_P_F(1,:)<objMax0(1))&(inv_P_F(1,:)>objMin0(1))...
                &(inv_P_F(2,:)<objMax0(2))&(inv_P_F(2,:)>objMin0(2))...
                &(inv_P_F(3,:)<objMax0(3))&(inv_P_F(3,:)>objMin0(3)));
        P_F(:,idxs) = [];
    end
    
end

end


function [J,g] = fuzzyMetricFCM(rt,C_F,C_M,N_CF,N_CM,N_CMTrim)  % the FCM-based metric for 3D scan registration

r = rt(1:3);
theta = norm(r);
if theta == 0
    theta = 10^(-6);
end
invth = 1/theta;
sinth = sin(theta);
costh = cos(theta);
M = invth*[0  -r(3)  r(2);r(3)  0  -r(1);-r(2)  r(1)  0];
R = [1 0 0;0 1 0;0 0 1] + M*sinth + M*M*(1 - costh);
t = rt(4:6);

C_MRnT = R*C_M + t*ones(1,N_CM);
e2 = @(a,b) sum(bsxfun(@minus,permute(a,[3 2 1]),permute(b,[2 3 1])).^2,3);
dist2 = e2(C_MRnT,C_F);
obj_CM = 1./sum(1./dist2);
if N_CM == N_CMTrim
    J = sum(obj_CM);
else
    [obj_CM_sort,CM_idxs] = sort(obj_CM);
    J = sum(obj_CM_sort(1:N_CMTrim));
end

if nargout > 1
    
    g = [0 0 0 0 0 0]';
    dC_MRnT = [0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];
    ones1xNCF = ones(1,N_CF);
    invth2 = invth*invth;
    term01 = invth*sinth;
    term02 = invth2*(1 - costh);
    term03 = invth2*(costh - term01);
    term04 = 2*invth2*invth2*(costh - 1) + invth2*term01;
    rxterm02 = r*term02;
    
    if N_CM == N_CMTrim
        
        for j = 1:N_CM
            
            C_Mj = C_M(:,j);
            rDotC_Mj = r(1)*C_Mj(1) + r(2)*C_Mj(2) + r(3)*C_Mj(3);
            term05 = rDotC_Mj*term02;
            term06 = rDotC_Mj*term04;
            C_Mjxterm01 = C_Mj*term01;
            
            term1_1 = r(1)*term06 - C_Mjxterm01(1) + (C_Mj(3)*r(2) - C_Mj(2)*r(3))*term03;
            term1_2 = r(2)*term06 - C_Mjxterm01(2) + (C_Mj(1)*r(3) - C_Mj(3)*r(1))*term03;
            term1_3 = r(3)*term06 - C_Mjxterm01(3) + (C_Mj(2)*r(1) - C_Mj(1)*r(2))*term03;
            
            dC_MRnT(1,1) = r(1)*term1_1 + C_Mj(1)*rxterm02(1) + term05;
            dC_MRnT(2,1) = r(2)*term1_1 + C_Mj(2)*rxterm02(1) + C_Mjxterm01(3);
            dC_MRnT(3,1) = r(3)*term1_1 + C_Mj(3)*rxterm02(1) - C_Mjxterm01(2);
            
            dC_MRnT(1,2) = r(1)*term1_2 + C_Mj(1)*rxterm02(2) - C_Mjxterm01(3);
            dC_MRnT(2,2) = r(2)*term1_2 + C_Mj(2)*rxterm02(2) + term05;
            dC_MRnT(3,2) = r(3)*term1_2 + C_Mj(3)*rxterm02(2) + C_Mjxterm01(1);
            
            dC_MRnT(1,3) = r(1)*term1_3 + C_Mj(1)*rxterm02(3) + C_Mjxterm01(2);
            dC_MRnT(2,3) = r(2)*term1_3 + C_Mj(2)*rxterm02(3) - C_Mjxterm01(1);
            dC_MRnT(3,3) = r(3)*term1_3 + C_Mj(3)*rxterm02(3) + term05;
            
            dd2 = 2*dC_MRnT*(C_MRnT(:,j)*ones1xNCF - C_F);
            g = g + (obj_CM(j)^2)*(dd2*(1./dist2(:,j).^2));
            
        end
        
    else
        
        for j = 1:N_CMTrim
            
            j1 = CM_idxs(j);
            C_Mj = C_M(:,j1);
            rDotC_Mj = r(1)*C_Mj(1) + r(2)*C_Mj(2) + r(3)*C_Mj(3);
            term05 = rDotC_Mj*term02;
            term06 = rDotC_Mj*term04;
            C_Mjxterm01 = C_Mj*term01;
            
            term1_1 = r(1)*term06 - C_Mjxterm01(1) + (C_Mj(3)*r(2) - C_Mj(2)*r(3))*term03;
            term1_2 = r(2)*term06 - C_Mjxterm01(2) + (C_Mj(1)*r(3) - C_Mj(3)*r(1))*term03;
            term1_3 = r(3)*term06 - C_Mjxterm01(3) + (C_Mj(2)*r(1) - C_Mj(1)*r(2))*term03;
            
            dC_MRnT(1,1) = r(1)*term1_1 + C_Mj(1)*rxterm02(1) + term05;
            dC_MRnT(2,1) = r(2)*term1_1 + C_Mj(2)*rxterm02(1) + C_Mjxterm01(3);
            dC_MRnT(3,1) = r(3)*term1_1 + C_Mj(3)*rxterm02(1) - C_Mjxterm01(2);
            
            dC_MRnT(1,2) = r(1)*term1_2 + C_Mj(1)*rxterm02(2) - C_Mjxterm01(3);
            dC_MRnT(2,2) = r(2)*term1_2 + C_Mj(2)*rxterm02(2) + term05;
            dC_MRnT(3,2) = r(3)*term1_2 + C_Mj(3)*rxterm02(2) + C_Mjxterm01(1);
            
            dC_MRnT(1,3) = r(1)*term1_3 + C_Mj(1)*rxterm02(3) + C_Mjxterm01(2);
            dC_MRnT(2,3) = r(2)*term1_3 + C_Mj(2)*rxterm02(3) - C_Mjxterm01(1);
            dC_MRnT(3,3) = r(3)*term1_3 + C_Mj(3)*rxterm02(3) + term05;
            
            dd2 = 2*dC_MRnT*(C_MRnT(:,j1)*ones1xNCF - C_F);
            g = g + (obj_CM_sort(j)^2)*(dd2*(1./dist2(:,j1).^2));
            
        end
        
    end
    
end

end


function [J,g] = fuzzyMetricGK(rt,C,P,N_C,N_P,N_PTrim,pInvCov)  % the  GK-based metric for 3D scan registration

r = rt(1:3);
theta = norm(r);
if theta == 0
    theta = 10^(-6);
end
invth = 1/theta;
sinth = sin(theta);
costh = cos(theta);

M = invth*[0  -r(3)  r(2);r(3)  0  -r(1);-r(2)  r(1)  0];
R = [1 0 0;0 1 0;0 0 1] + M*sinth + M*M*(1 - costh);
t = rt(4:6);

ones1xNP = ones(1,N_P);
PRnT = R*P + t*ones1xNP;
MD2 = zeros(N_C,N_P);
pD = zeros(3,N_C,N_P);

for i = 1:N_C
    dist = PRnT - C(:,i)*ones1xNP;  % 3 x N_CM
    pD_i = pInvCov(:,:,i)*dist;  % 3 x N_CM
    MD2(i,:) = sum(dist.*pD_i,1);
    pD(:,i,:) = pD_i;
end

obj_P = 1./sum(1./MD2,1);

if N_P == N_PTrim
    J = sum(obj_P);
else
    [obj_P_sort,P_idxs] = sort(obj_P);
    J = sum(obj_P_sort(1:N_PTrim));
end

if nargout > 1
    
    g = [0 0 0 0 0 0]';
    dPRnT = [0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];
    invth2 = invth*invth;
    term01 = invth*sinth;
    term02 = invth2*(1 - costh);
    term03 = invth2*(costh - term01);
    term04 = 2*invth2*invth2*(costh - 1) + invth2*term01;
    rxterm02 = r*term02;
    
    if N_P == N_PTrim
        
        for j = 1:N_P
            
            Pj = P(:,j);
            rDotPj = r(1)*Pj(1) + r(2)*Pj(2) + r(3)*Pj(3);
            term05 = rDotPj*term02;
            term06 = rDotPj*term04;
            Pjxterm01 = Pj*term01;
            
            term1_1 = r(1)*term06 - Pjxterm01(1) + (Pj(3)*r(2) - Pj(2)*r(3))*term03;
            term1_2 = r(2)*term06 - Pjxterm01(2) + (Pj(1)*r(3) - Pj(3)*r(1))*term03;
            term1_3 = r(3)*term06 - Pjxterm01(3) + (Pj(2)*r(1) - Pj(1)*r(2))*term03;
            
            dPRnT(1,1) = r(1)*term1_1 + Pj(1)*rxterm02(1) + term05;
            dPRnT(2,1) = r(2)*term1_1 + Pj(2)*rxterm02(1) + Pjxterm01(3);
            dPRnT(3,1) = r(3)*term1_1 + Pj(3)*rxterm02(1) - Pjxterm01(2);
            
            dPRnT(1,2) = r(1)*term1_2 + Pj(1)*rxterm02(2) - Pjxterm01(3);
            dPRnT(2,2) = r(2)*term1_2 + Pj(2)*rxterm02(2) + term05;
            dPRnT(3,2) = r(3)*term1_2 + Pj(3)*rxterm02(2) + Pjxterm01(1);
            
            dPRnT(1,3) = r(1)*term1_3 + Pj(1)*rxterm02(3) + Pjxterm01(2);
            dPRnT(2,3) = r(2)*term1_3 + Pj(2)*rxterm02(3) - Pjxterm01(1);
            dPRnT(3,3) = r(3)*term1_3 + Pj(3)*rxterm02(3) + term05;
            
            dd2 = 2*dPRnT*pD(:,:,j);
            g = g + (obj_P(j)^2)*(dd2*(1./MD2(:,j).^2));
            
        end
        
    else
        
        for j = 1:N_PTrim
            
            j1 = P_idxs(j);
            Pj = P(:,j1);
            rDotPj = r(1)*Pj(1) + r(2)*Pj(2) + r(3)*Pj(3);
            term05 = rDotPj*term02;
            term06 = rDotPj*term04;
            Pjxterm01 = Pj*term01;
            
            term1_1 = r(1)*term06 - Pjxterm01(1) + (Pj(3)*r(2) - Pj(2)*r(3))*term03;
            term1_2 = r(2)*term06 - Pjxterm01(2) + (Pj(1)*r(3) - Pj(3)*r(1))*term03;
            term1_3 = r(3)*term06 - Pjxterm01(3) + (Pj(2)*r(1) - Pj(1)*r(2))*term03;
            
            dPRnT(1,1) = r(1)*term1_1 + Pj(1)*rxterm02(1) + term05;
            dPRnT(2,1) = r(2)*term1_1 + Pj(2)*rxterm02(1) + Pjxterm01(3);
            dPRnT(3,1) = r(3)*term1_1 + Pj(3)*rxterm02(1) - Pjxterm01(2);
            
            dPRnT(1,2) = r(1)*term1_2 + Pj(1)*rxterm02(2) - Pjxterm01(3);
            dPRnT(2,2) = r(2)*term1_2 + Pj(2)*rxterm02(2) + term05;
            dPRnT(3,2) = r(3)*term1_2 + Pj(3)*rxterm02(2) + Pjxterm01(1);
            
            dPRnT(1,3) = r(1)*term1_3 + Pj(1)*rxterm02(3) + Pjxterm01(2);
            dPRnT(2,3) = r(2)*term1_3 + Pj(2)*rxterm02(3) - Pjxterm01(1);
            dPRnT(3,3) = r(3)*term1_3 + Pj(3)*rxterm02(3) + term05;
            
            dd2 = 2*dPRnT*pD(:,:,j1);
            g = g + (obj_P_sort(j)^2)*(dd2*(1./MD2(:,j1).^2));
            
        end
        
    end
    
end

end


function [R,t] = fuzzyPSReg_SS(C_F_fcm,AFPCD_fcm,C_F_gk,pInvCov_gk,N_C,P_M,trimRatio,gs_PM)

[C_F_fcm,P_M,tMidF,tMidM,scaleFactor] = normPoints(C_F_fcm,P_M);

C_F_gk = scaleFactor*(C_F_gk - tMidF);

sf2 = scaleFactor^2;
AFPCD_fcm = AFPCD_fcm*sf2;
pInvCov_gk = pInvCov_gk*sf2;

P_M = downsamplePoints(P_M,scaleFactor*gs_PM);
[C_M,~,~] = fuzzyClusteringFCM(P_M,N_C,100);

RRange = pi;
TRange = 0.5;
sigma_r = 0.01;
sigma_t = TRange/4;
epsilon = 0;
rt = globalSearch(C_F_fcm,C_M,AFPCD_fcm,trimRatio,RRange,TRange,sigma_r,sigma_t,epsilon);

N_PM = size(P_M,2);
N_PMTrim = round( N_PM*(1-trimRatio) );
options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','off');
rt = fminunc(@(rt)fuzzyMetricGK(rt,C_F_gk,P_M,N_C,N_PM,N_PMTrim,pInvCov_gk),rt,options);

r = rt(1:3);
theta = norm(r);
if theta == 0
    R = [1 0 0;0 1 0;0 0 1];
else
    R = axang2rotm([r'/theta  theta]);
end
t = tMidF - R*tMidM + rt(4:6)/scaleFactor;

end


function [P_F,P_M,tMidF,tMidM,scaleFactor] = normPoints(P_F,P_M)  % normalization of P_F and P_M

choice = 1;  % choice = {0, 1, 2, 3}

if choice == 0  % no change is applied to the two point sets
    
    tMidF = [0 0 0]';
    tMidM = [0 0 0]';
    scaleFactor = 1;
    fprintf('The two point sets are not normalized for this registration. \n');
    
else
    
    if choice == 1  % method 1
        tMidF = 0.5*(max(P_F,[],2) + min(P_F,[],2));
        tMidM = 0.5*(max(P_M,[],2) + min(P_M,[],2));
    elseif choice == 2  % method 2
        tMidF = sum(P_F,2)/size(P_F,2);
        tMidM = sum(P_M,2)/size(P_M,2);
    elseif choice == 3  % method 3
        tMidF = [0 0 0]';
        tMidM = [0 0 0]';
    else
        error('Incorrect choice for the normalization method. ')
    end
    
    P_F = P_F - tMidF;
    P_M = P_M - tMidM;
    scaleFactor = 1/max(max(abs(P_F),[],'all'), max(abs(P_M),[],'all'));
    P_F = scaleFactor*P_F;  % normalized fixed  set
    P_M = scaleFactor*P_M;  % normalized moving set
    
end

end


function rt = globalSearch(C_F,C_M,AFPCD,trimRatio,RRange,TRange,sigma_r,sigma_t,epsilon)

N_CF = size(C_F,2);
N_CM = size(C_M,2);
N_CMTrim = round( N_CM*(1-trimRatio) );
normCM = sqrt(sum(C_M.^2));
epsilonxN_CMTrim = epsilon*N_CMTrim;
AFPCDxN_CMTrim = AFPCD*N_CMTrim;
fprintf('The value of AFPCD x N_CMTrim is %f. \n', AFPCDxN_CMTrim);
sqrt3 = 1.73205080756888;
octree = [-1 -1 -1 -1  1  1  1  1;-1 -1  1  1 -1 -1  1  1;-1  1 -1  1 -1  1 -1  1];
ones1x8 = [1 1 1 1 1 1 1 1];

rt = [0 0 0 0 0 0]';
N_Qr = 6000;
Qr = Inf(N_Qr,6);
Infs1x6 = Qr(1,:);
Qr(1,1:3) = rt(1:3)';
Qr(1,4) = RRange;
Qr(1,6) = 0;
n_Qr = 1;

options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','off');
[rt,bestDistLoss] = fminunc(@(rt)fuzzyMetricFCM(rt,C_F,C_M,N_CF,N_CM,N_CMTrim),rt,options);
fprintf('The value of AFCCD x N_CMTrim (current best value of the FCM-based metric) \n');
fprintf('                                   after the first local convergence is %f. \n', bestDistLoss);

while (1)
    
    if bestDistLoss < AFPCDxN_CMTrim
        fprintf('Great! Now it has AFCCD < AFPCD (rho_fcm < 1). \n');
        break
    end
    if n_Qr == 0
        fprintf('The number of cubes in the queue is 0. \n');
        break
    end
    [minLB,minLB_index] = min(Qr(:,6));
    
    if bestDistLoss - minLB <= epsilonxN_CMTrim
        fprintf('The difference between the current best cost and the lowest lower bound is less than or equal to the threshold. \n');
        break
    end
    if Qr(minLB_index,4) < sigma_r
        fprintf('The divided cubes are sufficiently small. \n');
        break
    end
    
    rParentInfo = Qr(minLB_index,:);
    Qr(minLB_index:N_Qr-1,:) = Qr(minLB_index+1:N_Qr,:);
    Qr(N_Qr,:) = Infs1x6;
    n_Qr = n_Qr - 1;
    
    rRange = 0.5*rParentInfo(4);
    rKids = rParentInfo(1:3)'*ones1x8 + rRange*octree;
    sqrt3xrRange = sqrt3*rRange;
    
    for s = 1:8
        theta = norm(rKids(:,s));
        if theta - sqrt3xrRange > RRange
            continue
        end
        
        if theta == 0
            C_MR = C_M;
        else
            C_MR = axang2rotm([rKids(:,s)'/theta  theta])*C_M;  % rotated C_M
        end
        
        %%%%%% for lower bound
        [LB,~] = innerSearch(C_F,C_MR,N_CF,N_CM,N_CMTrim,AFPCDxN_CMTrim,normCM,...
                             bestDistLoss,rRange,TRange,sigma_t,epsilonxN_CMTrim);
        
        if ( LB >= bestDistLoss )||( LB >= AFPCDxN_CMTrim )
            continue;
        end
        
        %%%%%% for upper bound
        [UB,tBest] = innerSearch(C_F,C_MR,N_CF,N_CM,N_CMTrim,AFPCDxN_CMTrim,normCM,...
                                 bestDistLoss,0,TRange,sigma_t,epsilonxN_CMTrim);
        
        if UB < bestDistLoss
            rtNew = [rKids(:,s);tBest];
            fprintf('A better solution is found by BnB, now the value of AFCCD x N_CMTrim is %f. \n', UB);
            
            %%%%%% local convergence for refinement
            [rt,bestDistLoss] = fminunc(@(rt)fuzzyMetricFCM(rt,C_F,C_M,N_CF,N_CM,N_CMTrim),rtNew,options);
            fprintf('After refinement by local convergence, the value of AFCCD x N_CMTrim is %f. \n', bestDistLoss);
            
            if bestDistLoss < AFPCDxN_CMTrim
                break;
            end
            
            if n_Qr ~= 0
                Qr_New = Qr(Qr(:,6) < bestDistLoss,:);
                n_Qr = size(Qr_New,1);
                Qr(1:n_Qr,:) = Qr_New;
                Qr(n_Qr+1:N_Qr,:) = 10^6;
            end
        end
        
        n_Qr = n_Qr + 1;
        if n_Qr >= N_Qr
            Qr(N_Qr+1:N_Qr+100,:) = 10^6;
            N_Qr = N_Qr + 100;
        end
        Qr(n_Qr,1:3) = rKids(:,s)';
        Qr(n_Qr,4) = rRange;
        Qr(n_Qr,5) = UB;
        Qr(n_Qr,6) = LB;
    end
    
end

end


function [bound,tBest] = innerSearch(C_F,C_MR,N_CF,N_CM,N_CMTrim,AFPCDxN_CMTrim,normCM,...
    bestDistLoss,rRange,TRange,sigma_t,epsilonxN_CMTrim)

sqrt3 = 1.73205080756888;
octree = [-1 -1 -1 -1  1  1  1  1;-1 -1  1  1 -1 -1  1  1;-1  1 -1  1 -1  1 -1  1];
ones1x8 = [1 1 1 1 1 1 1 1];
ones1xNCM = ones(1,N_CM);
zeros1xNCM = 0*ones1xNCM;
onesNCFxNCM = ones(N_CF,N_CM);
e2 = @(a,b) sum(bsxfun(@minus,permute(a,[3 2 1]),permute(b,[2 3 1])).^2,3);

if rRange ~= 0
    gamma_r = 2*ones(N_CF,1)*sin(0.5*min(sqrt3*rRange,pi))*normCM;  % N_CF x N_CM
end

obj_CMr = zeros1xNCM;
if rRange == 0
    dist2_r = e2(C_MR,C_F);
    minD2r = min(dist2_r);
    minD2r_idxs = find(minD2r);
    obj_CMr(minD2r_idxs) = 1./sum(1./dist2_r(:,minD2r_idxs));
else
    dist_r = sqrt(e2(C_MR,C_F)) - gamma_r;
    minDr = min(dist_r);
    minDr(minDr < 0) = 0;
    minDr_idxs = find(minDr);
    obj_CMr(minDr_idxs) = 1./sum(1./dist_r(:,minDr_idxs).^2);
end

if N_CM == N_CMTrim
    UB = sum(obj_CMr);
else
    obj_CMr_sort = sort(obj_CMr);
    UB = sum(obj_CMr_sort(1:N_CMTrim));
end
if UB < bestDistLoss
    bound = UB;
else
    bound = bestDistLoss;
end

tBest = [0 0 0]';
N_Qt = 1000;
Qt = 10^6*ones(N_Qt,6);
largeNum1x6 = Qt(1,:);
Qt(1,1:3) = tBest';
Qt(1,4) = TRange;
Qt(1,5) = UB;
Qt(1,6) = 0;
n_Qt = 1;

while (1)
    
    if ( n_Qt == 0 )||( bound < AFPCDxN_CMTrim )
        break;
    end
    
    [minLB,minLB_idx] = min(Qt(:,6));  % read out the cube with the lowest lower bound
    
    if ( bound - minLB <= epsilonxN_CMTrim )||( Qt(minLB_idx,4) < sigma_t )
        break;
    end
    
    tParentInfo = Qt(minLB_idx,:);
    Qt(minLB_idx:N_Qt-1,:) = Qt(minLB_idx+1:N_Qt,:);
    Qt(N_Qt,:) = largeNum1x6;
    n_Qt = n_Qt - 1;
    
    tRange = 0.5*tParentInfo(4);
    tKids = tParentInfo(1:3)'*ones1x8 + tRange*octree;
    gamma_t = sqrt3*tRange*onesNCFxNCM;
    
    for s = 1:8
        C_MRnT = C_MR + tKids(:,s)*ones1xNCM;
        
        obj_CMrt = zeros1xNCM;
        obj_CMr  = zeros1xNCM;
        if N_CM == N_CMTrim
            if rRange == 0
                dist2_r = e2(C_MRnT,C_F);
                dist_rt = sqrt(dist2_r) - gamma_t;
                minDrt = min(dist_rt);
                minDrt(minDrt < 0) = 0;
                minDrt_idxs = find(minDrt);
                obj_CMrt(minDrt_idxs) = 1./sum(1./dist_rt(:,minDrt_idxs).^2);
                LB = sum(obj_CMrt);
                
                if ( LB >= bound )||( LB >= AFPCDxN_CMTrim )
                    continue;
                end
                
                minD2r = min(dist2_r);
                minD2r_idxs = find(minD2r);
                obj_CMr(minD2r_idxs) = 1./sum(1./dist2_r(:,minD2r_idxs));
                UB = sum(obj_CMr);
            else
                dist_r = sqrt(e2(C_MRnT,C_F)) - gamma_r;
                dist_rt = dist_r - gamma_t;
                minDrt = min(dist_rt);
                minDrt(minDrt < 0) = 0;
                minDrt_idxs = find(minDrt);
                obj_CMrt(minDrt_idxs) = 1./sum(1./dist_rt(:,minDrt_idxs).^2);
                LB = sum(obj_CMrt);
                
                if ( LB >= bound )||( LB >= AFPCDxN_CMTrim )
                    continue;
                end
                
                minDr = min(dist_r);
                minDr(minDr < 0) = 0;
                minDr_idxs = find(minDr);
                obj_CMr(minDr_idxs) = 1./sum(1./dist_r(:,minDr_idxs).^2);
                UB = sum(obj_CMr);
            end
        else
            if rRange == 0
                dist2_r = e2(C_MRnT,C_F);
                minD2r = min(dist2_r);
                minD2r_idxs = find(minD2r);
                obj_CMr(minD2r_idxs) = 1./sum(1./dist2_r(:,minD2r_idxs));
                dist_rt = sqrt(dist2_r) - gamma_t;
            else
                dist_r = sqrt(e2(C_MRnT,C_F)) - gamma_r;
                minDr = min(dist_r);
                minDr(minDr < 0) = 0;
                minDr_idxs = find(minDr);
                obj_CMr(minDr_idxs) = 1./sum(1./dist_r(:,minDr_idxs).^2);
                dist_rt = dist_r - gamma_t;
            end
            minDrt = min(dist_rt);
            minDrt(minDrt < 0) = 0;
            minDrt_idxs = find(minDrt);
            obj_CMrt(minDrt_idxs) = 1./sum(1./dist_rt(:,minDrt_idxs).^2);
            [obj_CMr_sort,CMr_idxs] = sort(obj_CMr);
            LB = sum(obj_CMrt( CMr_idxs(1:N_CMTrim) ));
            
            if ( LB >= bound )||( LB >= AFPCDxN_CMTrim )
                continue;
            end
            
            UB = sum(obj_CMr_sort(1:N_CMTrim));
        end
        
        if UB < bound
            bound = UB;
            tBest = tKids(:,s);
            if bound < AFPCDxN_CMTrim
                break;
            end
        end
        
        n_Qt = n_Qt + 1;
        if n_Qt >= N_Qt
            Qt(N_Qt+1:N_Qt+100,:) = 10^6;
            N_Qt = N_Qt + 100;
        end
        Qt(n_Qt,1:3) = tKids(:,s)';
        Qt(n_Qt,4) = tRange;
        Qt(n_Qt,5) = UB;
        Qt(n_Qt,6) = LB;
    end
    
end

end


function plot4checking(P_F,P_M,R,t)

if nargin == 4
    P_MRt = R*P_M + t;
elseif nargin == 2
    P_MRt = P_M;
else
    error('Incorrect setting ')
end
figure
hold on
plot3(P_F(1,:),P_F(2,:),P_F(3,:),'.r','MarkerSize',4);
plot3(P_MRt(1,:),P_MRt(2,:),P_MRt(3,:),'.b','MarkerSize',4);

end

