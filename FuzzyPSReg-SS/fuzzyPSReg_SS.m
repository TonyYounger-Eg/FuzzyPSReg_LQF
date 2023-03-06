%%%%%% The main function
function [R,t,timeReg,timeFuzzyClustering,timeCoarseReg,timeFineReg] = fuzzyPSReg_SS(P_F,P_M,...
    noisyF,noisyM,trimRatio,N_C,gsF,gsM,gsF_fine,gsM_fine,RRange,TRange,sigma_r,sigma_t,epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FuzzyPSReg-SS: Fuzzy cluster-based point set registration for similar-sized range scans.
%   This function estimates the rigid transformation parameters, R and t, to register the two given 3D
%   point sets P_F and P_M, where P_F (3 x N_PF) is the fixed set and P_M (3 x N_PM) is the moving set.
%
%   Inputs:
%     'P_F'                  fixed  (target) point set (3 x N_PF).
%
%     'P_M'                  moving (source) point set (3 x N_PM).
%
%     'noisyF'               1/0 denoting the condition of P_F with/without noisy outliers.
%
%     'noisyM'               1/0 denoting the condition of P_M with/without noisy outliers.
%
%     'trimRatio'            trimming ratio for registering partially overlapping point sets,
%                            trimRatio = 0 for no trimming cases.
%
%     'N_C'                  number of fuzzy clusters of each point set in the coarse registration.
%
%     'gsF'                  grid step of the box grid filter to downsample P_F for improving the efficiency
%                            of the fuzzy clustering, gsF = 0 for no downsampling.
%
%     'gsM'                  grid step of the box grid filter to downsample P_M for improving the efficiency
%                            of the fuzzy clustering, gsM = 0 for no downsampling.
%
%     'gsF_fine'             grid step of the box grid filter to downsample P_F for improving the efficiency
%                            of the fine registration, gsF_fine = 0 for no downsampling.
%
%     'gsM_fine'             grid step of the box grid filter to downsample P_M for improving the efficiency
%                            of the fine registration, gsM_fine = 0 for no downsampling.
%
%     'RRange'               half side-length of the initial rotation    cube [-RRange, RRange]^3 (radians).
%
%     'TRange'               half side-length of the initial translation cube [-TRange, TRange]^3.
%
%     'sigma_r'              half side-length of the minimum rotation    cube for stopping the global search.
%
%     'sigma_t'              half side-length of the minimum translation cube for stopping the inner  search.
%
%     'epsilon'              threshold value for stopping the global search/inner search.
%
%   Outputs:
%     'R'                    3 x 3 rotation    matrix to transform P_M for registration.
%
%     't'                    3 x 1 translation vector to transform P_M for registration.
%
%     'timeReg'              total time cost (sec) of the registration.
%
%     'timeFuzzyClustering'  time cost (sec) of the fuzzy clustering.
%
%     'timeCoarseReg'        time cost (sec) of the coarse registration.
%
%     'timeFineReg'          time cost (sec) of the fine registration.
%
%   References:
%     [1] Q. Liao, D. Sun, and H. Andreasson, "FuzzyPSReg: Strategies of Fuzzy Cluster-based
%                  Point Set Registration", IEEE Transactions on Robotics, 2021.
%
%     [2] Q. Liao, D. Sun, and H. Andreasson, "Point Set Registration for 3D Range Scans
%                  Using Fuzzy Cluster-based Metric and Efficient Global Optimization",
%                  IEEE Transactions on Pattern Analysis and Machine Intelligence, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[D_PF,N_PF] = size(P_F);
[D_PM,N_PM] = size(P_M);

if ( N_PF == 0 )||( N_PM == 0 )||( D_PF ~= 3 )||( D_PM ~= 3 )
    error('Something wrong with the input point sets. ');
end
if ( noisyF ~= 0 )&&( noisyF ~= 1 )
    error('Incorrect setting for noisyF. ');
end
if ( noisyM ~= 0 )&&( noisyM ~= 1 )
    error('Incorrect setting for noisyM. ');
end
if ( trimRatio < 0 )||( trimRatio >= 1 )
    error('Incorrect setting for trimRatio. ');
end

fprintf('Registration starts. \n\n');

[P_F,P_M,tMidF,tMidM,scaleFactor] = normPoints(P_F,P_M);  % normalize the point sets to fit the cube [-1, 1]^3

%%%%%% Point set downsampling to improve efficiency
if noisyF == 0
    Pds_F = downsamplePoints(P_F,scaleFactor*gsF);
else
    Pds_F = P_F;  % no downsampling is applied if P_F has noisy outliers
end
if noisyM == 0
    Pds_M = downsamplePoints(P_M,scaleFactor*gsM);
else
    Pds_M = P_M;  % no downsampling is applied if P_M has noisy outliers
end

fprintf('Fuzzy clustering (FCM and GK clustering) ... \n');
%%%%%%%%%%%%%%%%%%%%%%%%%% Fuzzy clustering (N_CF = N_CM = N_C, see Strategy 1 in reference [2])
N_CF = N_C;  % number of fuzzy clusters of P_F
N_CM = N_C;  % number of fuzzy clusters of P_M

%%%%%% FCM clustering for the two point sets
tic
maxIter = 100;
[C_F,AFPCD_F,Pp_F,U_F] = fuzzyClusteringFCM(Pds_F,N_CF,maxIter,noisyF);
[C_M,AFPCD_M,Pp_M,U_M] = fuzzyClusteringFCM(Pds_M,N_CM,maxIter,noisyM);

%%%%%% GK clustering for one of the point sets
maxIter = 30;
if AFPCD_F >= AFPCD_M
    baseIsF = 1;  % 'baseIsF' relates to Strategy 1 in reference [2]
    AFPCD = AFPCD_F;
    if noisyF == 0
        [C_fine,~,pInvCov] = fuzzyClusteringGK(Pds_F,N_C,maxIter,U_F);
    else
        [C_fine,~,pInvCov] = fuzzyClusteringGK(Pp_F, N_C,maxIter,U_F);
    end
else
    baseIsF = 0;
    AFPCD = AFPCD_M;
    if noisyM == 0
        [C_fine,~,pInvCov] = fuzzyClusteringGK(Pds_M,N_C,maxIter,U_M);
    else
        [C_fine,~,pInvCov] = fuzzyClusteringGK(Pp_M, N_C,maxIter,U_M);
    end
    tmp = C_M;
    C_M = C_F;
    C_F = tmp;
end
timeFuzzyClustering = toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Fuzzy clustering is complete. \n\n');

fprintf('Coarse registration (using the FCM-based metric) ... \n');
%%%%%%%%%%%%%%%%%%%%%%%%%% Coarse registration
[rt,timeCoarseReg] = globalSearch(C_F,C_M,AFPCD,trimRatio,RRange,TRange,sigma_r,sigma_t,epsilon);
%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Coarse registration is complete. \n\n');

fprintf('Fine registration (using the GK-based metric) ... \n');
%%%%%%%%%%%%%%%%%%%%%%%%%% Fine registration
if baseIsF == 1
    if noisyM == 0
        P_fine = downsamplePoints(P_M, scaleFactor*gsM_fine);
    else
        P_fine = downsamplePoints(Pp_M,scaleFactor*gsM_fine);
    end
else
    if noisyF == 0
        P_fine = downsamplePoints(P_F, scaleFactor*gsF_fine);
    else
        P_fine = downsamplePoints(Pp_F,scaleFactor*gsF_fine);
    end
end

%%%%%% Equation (18) in reference [2], which is just a suggestion.
if trimRatio == 0
    trimRatio_fine = 0;
elseif trimRatio < 0.1
    trimRatio_fine = 0.75*trimRatio + 0.075;
elseif trimRatio < 0.2
    trimRatio_fine = 0.5*trimRatio + 0.1;
else
    trimRatio_fine = trimRatio;
end

N_P_fine = size(P_fine,2);
N_PTrim_fine = round( N_P_fine*(1-trimRatio_fine) );
options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','off');

tic
rt = fminunc(@(rt)fuzzyMetricGK(rt,C_fine,P_fine,N_C,N_P_fine,N_PTrim_fine,pInvCov),rt,options);
timeFineReg = toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Fine registration is complete. \n\n');

if baseIsF == 1
    theta = norm(rt(1:3));
    if theta == 0
        R = [1 0 0;0 1 0;0 0 1];
    else
        R = axang2rotm([rt(1:3)'/theta, theta]);
    end
else  % calculate the inverse transformation if baseIsF is 0
    theta = norm(rt(1:3));
    if theta == 0
        R = [1 0 0;0 1 0;0 0 1];
        rt = -rt;
    else
        R = axang2rotm([rt(1:3)'/theta, theta]);
        R = R';
        rt = [-rt(1:3); -R*rt(4:6)];
    end
end
t = tMidF - R*tMidM + rt(4:6)/scaleFactor;

timeReg = timeFuzzyClustering + timeCoarseReg + timeFineReg;

fprintf('Registration is complete, please see the following results. \n');
%%%%%%%%%%%%%%%%%%%%%%%%%% Display of results
fprintf('The rotation matrix is: \n');
R
fprintf('The translation vector is: \n');
t
fprintf('The total time cost (sec) is: \n');
timeReg
%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function [C,AFPCD,Pp,U] = fuzzyClusteringFCM(P,N_C,maxIter,noisy)  % FCM clustering

[D_P,N_P] = size(P);

if ( N_C < 2 )||( N_C > N_P )
    error('Incorrect setting for N_C. ');
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
        U2D2 = U2.*dist2;
        obj = sum(U2D2,'all');
        
        if abs(obj - obj_1) <= minDiff
            break;
        end
        
    end
    
else
    
    dist2 = 0*U;  % N_C x N_P, squared distance matrix
    for i = 1:N_C
        dist2(i,:) = sum((C(:,i)*ones1xNP - P).^2,1);  % squared distance matrix
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
        U2D2 = U2.*dist2;
        obj = sum(U2D2,'all');
        
        if abs(obj - obj_1) <= minDiff
            break;
        end
        
    end
    
end

if noisy == 0
    
    AFPCD = obj/N_P;
    Pp = P;
    
else  % noisy outlier pruning (see reference [2] for details)
    
    %%%%%% the first step
    eta2 = sum(U2D2,2)./rowSumU2;
    normD2 = dist2./(eta2*ones1xNP);
    minND2 = min(normD2);
    minND2_idxs = find(minND2 < 1);
    N_Pp = size(minND2_idxs,2);
    Pp = P(:,minND2_idxs);
    U = U(:,minND2_idxs);
    U2D2p = U2D2(:,minND2_idxs);
    obj_Pp = sum(U2D2p,1);
    %%%%%% the second step
    pruneRatio = 0.15;  % this value can be chosen by users
    [obj_Pp_sort,Pp_idxs] = sort(obj_Pp);
    N_Pp = round( N_Pp*(1-pruneRatio) );
    AFPCD = sum(obj_Pp_sort(1:N_Pp))/N_Pp;
    Pp = Pp(:,Pp_idxs(1:N_Pp));  % pruned point set
    U = U(:,Pp_idxs(1:N_Pp));
    
end

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


function [rt,timeCoarseReg] = globalSearch(C_F,C_M,AFPCD,trimRatio,RRange,TRange,sigma_r,sigma_t,epsilon)

tic
N_CF = size(C_F,2);
N_CM = size(C_M,2);
N_CMTrim = round( N_CM*(1-trimRatio) );
normCM = sqrt(sum(C_M.^2));
epsilonxN_CMTrim = epsilon*N_CMTrim;
AFPCDxN_CMTrim = AFPCD*N_CMTrim;
fprintf('The value of AFPCD x N_CMTrim is %f. \n', AFPCDxN_CMTrim);
sqrt3 = 1.73205080756888;  % sqrt(3)
octree = [-1 -1 -1 -1  1  1  1  1;-1 -1  1  1 -1 -1  1  1;-1  1 -1  1 -1  1 -1  1];
ones1x8 = [1 1 1 1 1 1 1 1];

rt = [0 0 0 0 0 0]';  % 3D rotation parameters (axis-angle representation) and 3D translation parameters
N_Qr = 6000;
Qr = 10^6*ones(N_Qr,6);
largeNum1x6 = Qr(1,:);
Qr(1,1:3) = rt(1:3)';
Qr(1,4) = RRange;
Qr(1,6) = 0;
n_Qr = 1;

options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','off');
[rt,bestDistLoss] = fminunc(@(rt)fuzzyMetricFCM(rt,C_F,C_M,N_CF,N_CM,N_CMTrim),rt,options);
%%%%%% 'bestDistLoss' is the current best value of the FCM-based metric (AFCCD x N_CMTrim).
fprintf('The value of AFCCD x N_CMTrim (current best value of the FCM-based metric) \n');
fprintf('                                   after the first local convergence is %f. \n', bestDistLoss);
% n_div = 0;
% minLB = 0;
% stopCondition = 0;
% cubeInfo = zeros(10000,6);

while (1)
    
    % n_div = n_div + 1;
    % cubeInfo(n_div,6) = toc;
    if bestDistLoss < AFPCDxN_CMTrim
        fprintf('Great! Now it has AFCCD < AFPCD (rho_fcm < 1). \n');
        % stopCondition = 1;
        % cubeInfo(n_div,1:5) = [0  bestDistLoss  minLB  n_Qr  stopCondition];
        break
    end
    if n_Qr == 0
        fprintf('The number of cubes in the queue is 0. \n');
        % stopCondition = 2;
        % cubeInfo(n_div,1:5) = [0  bestDistLoss  minLB  n_Qr  stopCondition];
        break
    end
    [minLB,minLB_idx] = min(Qr(:,6));  % read out the cube with the lowest lower bound
    % cubeInfo(n_div,1:4) = [0  bestDistLoss  minLB  n_Qr];
    if bestDistLoss - minLB <= epsilonxN_CMTrim
        fprintf('The difference between the current best cost and the lowest lower bound is less than or equal to the threshold. \n');
        % stopCondition = 3;
        % cubeInfo(n_div,5) = stopCondition;
        break
    end
    if Qr(minLB_idx,4) < sigma_r
        fprintf('The divided cubes are sufficiently small. \n');
        % stopCondition = 4;
        % cubeInfo(n_div,5) = stopCondition;
        break
    end
    
    rParentInfo = Qr(minLB_idx,:);
    Qr(minLB_idx:N_Qr-1,:) = Qr(minLB_idx+1:N_Qr,:);
    Qr(N_Qr,:) = largeNum1x6;
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
            % cubeInfo(n_div,1) = UB;
            
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

timeCoarseReg = toc;
% cubeInfo = cubeInfo(1:n_div,:);

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

