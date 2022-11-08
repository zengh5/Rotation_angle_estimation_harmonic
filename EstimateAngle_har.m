% -------------------------------------------------------------------------
%   Written by Hui Zeng and Kun Yu
%   Copyright 2022 by Hui Zeng and Kun Yu
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The end-user understands that the program was
% developed for research purposes and is advised not to rely exclusively on
% the program for any reason.
% -------------------------------------------------------------------------
% Contact: zengh5@mail2.sysu.edu.cn
% -------------------------------------------------------------------------
% Input:       P .... input image spectrum, e.g., Txx
%           JPEG .... flag for whether the image undergone JPEG compression
%          N_har ....the highest order of the harmonics to be aggregated
% Output:  angle .... estimated angle
% -------------------------------------------------------------------------
% heavily borrowed from the following code:
% https://github.com/ChenglongChen/image-rotation-angle-estimation
function angle = EstimateAngle_har(P,JPEG, N_har)
% 1. Since the spectrum is centrosymmetry, we can add P with its 90-degree
% rotated version. This way we only need search peaks from one arc instead
% of four.
[~,N]=size(P);
P_rot = imrotate(P, 90);
% P_rot = [zeros(1,N) ; P_rot(1:(N-1),:)];  % note, this is wrong.
P_rot = [P_rot(N,:) ; P_rot(1:(N-1),:)];
P = P_rot + P;

% 2. avoiding compoments near DC with square window of size 9 x 9 (8 x 8 indeed)
% shift it to [-0.5,0.5]^2
P = fftshift(P);
% remove components near DC
w = 4; %w = 4;
x0 = N/2; y0 = N/2;
P(y0-w+1:y0+w,x0-w+1:x0+w)=0;
% shift it back to [0,1]^2
P = fftshift(P);
% figure, imshow(P, [])

% 3. Core routine---------------
% nearest neighbor around the arcs in the sampling grid
Neighbor = 3;    % Recommendation: using smaller neighborhood, e.g., Neighbor = 2, for small images, e.g., 64 by 64 images.  
angle = ArcEstimate_1_zh(P,JPEG,Neighbor, N_har);
end


% ----------------------------------------------------------------------- %
function EstimatedAngle = ArcEstimate_1_zh(P,JPEG,Neighbor, N_har)
    [~,N]=size(P);

    % peak searching around arcs  -----------------------------------------
    theta_step = 0.2;
    AngleRange = -89:theta_step:-1;
    c = cosd(AngleRange)';
    s = sind(AngleRange)';  
    NormalizedFlag = 1;
    
    % image origins from bottom left
    % -- arc 2: 90-degree arc center at upper-right (1,1) 
    x0 = 1-c;
    y0 = 1+s;
    figure, imshow(zeros(N)), hold on,
    plot(round(y0*N)+1, round((1-x0)*N)+1)
    
    [X2,Y2] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag);
    
    Fx = [X2];
    % Matlab origins from the upper-left
    Fy = [N+1-Y2];
    
    % remove duplicate rows  ----------
    Fxy = unique([Fx,Fy], 'rows');
    Fx = Fxy(:,1);
    Fy = Fxy(:,2);
    
    % get the magnitude aroud the arcs  ----------
    linearInd = sub2ind([N,N],Fy,Fx);
    Q = P(linearInd);
%     figure, plot(Q)
    
    % peak searching  ----------
    % since we only search one arc, the code is significantly simplified
    [fx,fy, EstimatedAngle] = FindPeak_harmo(P, Q, N, Fx, Fy, JPEG, N_har);   
end


% ----------------------------------------------------------------------- %
function [X,Y] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag)

if(NormalizedFlag==0.5)
    x0 = x0+0.5;
    y0 = y0+0.5;
end
if(Neighbor==2)
    % 2 nearest neighbor of x0 in the sampling grid
    x1 = ceil(x0*N+1);
    x2 = floor(x0*N+1);
    
    % 2 nearest neighbor of y0 in the sampling grid
    y1 = ceil(y0*N+1);
    y2 = floor(y0*N+1);
    
    % avoiding out of boundry
    x1(x1<1) = 1;
    x2(x2<1) = 1;
    y1(y1<1) = 1;
    y2(y2<1) = 1;
    
    % 2x2 nearest neighbor of arc l^(1,0) in the sampling grid
    X = [x1;x1;x2;x2];
    Y = [y1;y2;y1;y2];
elseif(Neighbor==3)
    % 3 nearest neighbor of x0 in the sampling grid
    x1 = round(x0*N+1);
    x2 = x1-1;
    x3 = x1+1;
    
    % 3 nearest neighbor of y0 in the sampling grid
    y1 = round(y0*N+1);
    y2 = y1-1;
    y3 = y1+1;
    
    % avoiding out of boundry
    x1(x1<1) = 1;
    x2(x2<1) = 1;
    y1(y1<1) = 1;
    y2(y2<1) = 1;
    
    x3(x3>N) = N;
    y3(y3>N) = N;
    
    % 3x3 nearest neighbor of arc l^(1,0) in the sampling grid
    X = [x1;x1;x1;x2;x2;x2;x3;x3;x3];
    Y = [y1;y2;y3;y1;y2;y3;y1;y2;y3];
end
end


function [fx,fy, EstimatedAngle] = FindPeak_harmo(P,Q,N,X,Y,JPEG, N_har)
    N_candi = 5;    % It doesn't matter what you set N_candi = 5, 10 or something else
    EstimatedAngle = 0;
    [magnitude,ind]=sort(Q(:),'descend');
    if(JPEG==0)
        bestmagnitude = 0;
        for i = 1 : N_candi
            % 1. Candidate peak location
            fx = X(ind(i));
            fy = Y(ind(i));
            theta = atand((fy-1)/(N+1-fx));
            % project it to the Arch
            fx_proj = N+1-N*cosd(theta);
            fy_proj = 1+N*sind(theta);
            
            % 2. Expected locations of the 2-5th order harmonics
            fx_Harmo(1) = mod(1 + 2 * (fx_proj-1), N);
            fy_Harmo(1) = mod(1 + 2 * (fy_proj-1), N);
            fx_Harmo(2) = mod(1 + 3 * (fx_proj-1), N);
            fy_Harmo(2) = mod(1 + 3 * (fy_proj-1), N);
            fx_Harmo(3) = mod(1 + 4 * (fx_proj-1), N);
            fy_Harmo(3) = mod(1 + 4 * (fy_proj-1), N);
            fx_Harmo(4) = mod(1 + 5 * (fx_proj-1), N);
            fy_Harmo(4) = mod(1 + 5 * (fy_proj-1), N);
            
            % 3. Search the peaks around the expected locations of the harmonics
            for i_har = 2 : N_har
                x1 = round(fx_Harmo(i_har-1));
                y1 = round(fy_Harmo(i_har-1));
                if x1 == 0
                    x1 = N;
                end
                x2 = x1-1; x3 = x1+1;
                
                if y1 == 0
                    y1 = N;
                end
                y2 = y1-1; y3 = y1+1;
                
                % avoiding out of boundry
                x2(x2<1) = 1; y2(y2<1) = 1;
                x3(x3>N) = N; y3(y3>N) = N;
                
                % 3x3 nearest neighbor of arc l^(1,0) in the sampling grid
                neighbors = [P(y1,x1), P(y1,x2), P(y1,x3), P(y2,x1),P(y2,x2), ...
                             P(y2,x3), P(y3,x1), P(y3,x2), P(y3,x3)];
                V_harmo(i_har-1) = max(neighbors);               
            end  
            
            % 4 Aggregate the magnitudes of harmonics
            if N_har == 1
                % Without harmonics, basically reduce to Chen et al.'s method
                magnitude_combine(i) = magnitude(i);
            else
                magnitude_combine(i) = magnitude(i) + sum(V_harmo);
            end
            
            % 5 Update the best peak according to the aggregated magnitude
            if magnitude_combine(i) > (bestmagnitude+1e-4)
                bestmagnitude = magnitude_combine(i);
                bestfx = fx;
                bestfy = fy;
                EstimatedAngle = theta;
            end
            
        end  
        % normalizing to [0,1)^2
        fx = (bestfx-1)/N;
        fy = (bestfy-1)/N;        
    else        
        A = 1; % To do
    end
end
