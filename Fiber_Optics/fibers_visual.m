function fibers_visual(numFibers,psi,n_vect,d,Lam)

% example script: fibers_visual(2,60,[1.66,1.65,1.51],10000,550); ***Variables for two aspirin fibers.
% numFibers = number of fibers
% psi = total maximum misalignment angle (degrees)
% n_vect = refractive index vector [n_x, n_y, n_z] of lamella corresponding to a twist angle of zero
% ***In this arbitrary system, n_x is the radial direction, n_z is the direction out of paper.***
% d = thickness of one lamellae, in nm
% Lam = wavelength of light, in nm
% The imaginary part should be positive.

% fiber dimensions & refractive index
a = 1; 
b = 5; 
c = 30; 
psi = psi*pi/180; 
psi_0 = psi/(numFibers-1);
epsilon_c = diag(n_vect.^2); 
M = eye(4);

% define overlapping
overlapRegions = cell(numFibers, 1); 
for i = 1:numFibers
    overlapRegions{i} = polyshape(); 
end

% define optical properties for each layer
cbValues = zeros(numFibers,1); 
lbValues = zeros(numFibers,1);

figure;
hold on;
axis equal;
view(3);
colorbar
blueWhiteRed = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
                ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
blackToGreen = [zeros(256,1), linspace(0,1,256)', zeros(256,1)];

% visualizing fibers
        for k = 1:numFibers
            baseZ = (k-1) * a;
            tilt = (k-1) * psi_0;
            corners = [-c/2, -b/2; c/2, -b/2; c/2, b/2; -c/2, b/2];
            Ro = [cos(tilt), -sin(tilt); sin(tilt), cos(tilt)];
            rotatedCorners = corners * Ro;
            currentPoly = polyshape(rotatedCorners(:,1), rotatedCorners(:,2));
            vertices = [rotatedCorners, baseZ * ones(4,1); rotatedCorners, (baseZ + a) * ones(4,1)];
            faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
            patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'none', 'EdgeColor', 'k');
            
            % overlapping
            if k == 1
                overlapRegions{1} = currentPoly;
            else
                for j = k-1:-1:1
                    if ~isempty(overlapRegions{j}.Vertices)
                        newOverlap = intersect(overlapRegions{j}, currentPoly);
                        if ~isempty(newOverlap.Vertices)
                            if j+1 <= numFibers
                                overlapRegions{j+1} = union(overlapRegions{j+1}, newOverlap);
                            end
                        end
                    end
                end
                overlapRegions{1} = union(overlapRegions{1}, currentPoly);
            end
        end
        
% Mueller matrix calculation      
        for j = 1:numFibers
            t = -psi*((j-1)/(numFibers-1)); % average splay angle; why -psi? To be consistent with Xiaoyan and Melissa's convention of psi/handedness correlation
            R = rotz(t)*rotx(0); % rotation matrix
            epsilon = R*epsilon_c*R.';
            L  = (sqrt(epsilon(1,1))-sqrt(epsilon(2,2)))*2*pi*d/Lam; % LB calculated from refractive indices at each pixel
            Lp = (sqrt((epsilon(1,1)+epsilon(1,2)+epsilon(2,1)+epsilon(2,2))./2)-...
                   sqrt((epsilon(1,1)-epsilon(1,2)-epsilon(2,1)+epsilon(2,2))./2))*2*pi*d/Lam; % LB'
            LB = real(L);
            LD = imag(L);
            LBp = real(Lp);
            LDp = imag(Lp);
            Mat = exp(-sqrt(LD^2+LDp^2))*expm([0,-LD,-LDp,0;-LD,0,0,LBp;-LDp,0,0,-LB;0,-LBp,LB,0]);
            %Mat = exp(-sqrt(LD^2+LDp^2))*expm([0,-LD,-LDp,0;-LD,0,0,-LBp;-LDp,0,0,LB;0,LBp,-LB,0]);
            M = Mat * M; % Mueller matrix;
            J = nearestJones(M); % Jones matrix
            O = jonesAnisotropy(J);
            LB_M = real(1i.*O.*( J(1,1,:) - J(2,2,:) ));
            LBp_M = real(1i.*O.*( J(1,2,:) + J(2,1,:) ));
            lb = sqrt(LB_M.^2 + LBp_M.^2); % linear retardance
            lbValues(j) = lb;
            cb = (real(O.*( J(1,2,:) - J(2,1,:) ))); % circular retardance
            cbValues(j) = cb;
        end

% cb visualization
colormap(blueWhiteRed);
caxis([-0.5 0.5]);
    for j = 1:numFibers
        if ~isempty(overlapRegions{j}.Vertices)
            x = overlapRegions{j}.Vertices(:,1);
            y = overlapRegions{j}.Vertices(:,2);
            z = a * numFibers * ones(size(x)); 
            fill3(x, y, z, cbValues(j), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end
disp(cbValues);

% lb visualization, need to draw the fibers in another panel
figure;
hold on;
axis equal;
view(3);
colorbar;
colormap(blackToGreen);
caxis([0 1.57]); % 0 to pi/2
    for k = 1:numFibers
            baseZ = (k-1) * a;
            tilt = (k-1) * psi_0;
            corners = [-c/2, -b/2; c/2, -b/2; c/2, b/2; -c/2, b/2];
            Ro = [cos(tilt), -sin(tilt); sin(tilt), cos(tilt)];
            rotatedCorners = corners * Ro;
            %currentPoly = polyshape(rotatedCorners(:,1), rotatedCorners(:,2));
            vertices = [rotatedCorners, baseZ * ones(4,1); rotatedCorners, (baseZ + a) * ones(4,1)];
            faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
            patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'none', 'EdgeColor', 'k');
    end

    for j = 1:numFibers
        if ~isempty(overlapRegions{j}.Vertices)
            x = overlapRegions{j}.Vertices(:,1);
            y = overlapRegions{j}.Vertices(:,2);
            z = a * numFibers * ones(size(x));  
            fill3(x, y, z, lbValues(j), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end
disp(lbValues);

% save the lb and cb values
fid = fopen('fiber_values.txt', 'wt');
if fid == -1
    error('Cannot open file for writing.');
end
fprintf(fid, 'LB Values:\n');
fprintf(fid, '%f\n', lbValues); 
fprintf(fid, 'CB Values:\n');
fprintf(fid, '%f\n', cbValues); 
fclose(fid); 

end

% \\ LOCAL FUNCTIONS \\

    function R = MMrot(t)
        R = [1,0,0,0;0,cos(2*t),-sin(2*t),0;0,sin(2*t),cos(2*t),0;0,0,0,1]; % rotation for MMI
    end

    function matrix = rotx(angle)
        matrix = [1,0,0;0,cos(angle),-sin(angle);0,sin(angle),cos(angle)]; % rotation matrix on x
    end

    function matrix = rotz(angle)
        matrix = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1]; % rotation matrix on z
    end


function J = MJ2J(M)  % Mueller-Jones to Jones
J(1,1,:) = ((M(1,1,:)+M(1,2,:)+M(2,1,:)+M(2,2,:))/2).^(1/2);
k = 1./(2.*J(1,1,:));
J(1,2,:) = k.*(M(1,3,:)+M(2,3,:)-1i.*(M(1,4,:)+M(2,4,:)));
J(2,1,:) = k.*(M(3,1,:)+M(3,2,:)+1i.*(M(4,1,:)+M(4,2,:)));
J(2,2,:) = k.*(M(3,3,:)+M(4,4,:)+1i.*(M(4,3,:)-M(3,4,:)));
end


function C = M2Cov(M) % Mueller to Cloude covariance
C(1,1,:) = M(1,1,:) + M(1,2,:) + M(2,1,:) + M(2,2,:);
C(1,2,:) = M(1,3,:) + M(1,4,:)*1i + M(2,3,:) + M(2,4,:)*1i;
C(1,3,:) = M(3,1,:) + M(3,2,:) - M(4,1,:)*1i - M(4,2,:)*1i;
C(1,4,:) = M(3,3,:) + M(3,4,:)*1i - M(4,3,:)*1i + M(4,4,:);
C(2,1,:) = M(1,3,:) - M(1,4,:)*1i + M(2,3,:) - M(2,4,:)*1i;
C(2,2,:) = M(1,1,:) - M(1,2,:) + M(2,1,:) - M(2,2,:);
C(2,3,:) = M(3,3,:) - M(3,4,:)*1i - M(4,3,:)*1i - M(4,4,:);
C(2,4,:) = M(3,1,:) - M(3,2,:) - M(4,1,:)*1i + M(4,2,:)*1i;
C(3,1,:) = M(3,1,:) + M(3,2,:) + M(4,1,:)*1i + M(4,2,:)*1i;
C(3,2,:) = M(3,3,:) + M(3,4,:)*1i + M(4,3,:)*1i - M(4,4,:);
C(3,3,:) = M(1,1,:) + M(1,2,:) - M(2,1,:) - M(2,2,:);
C(3,4,:) = M(1,3,:) + M(1,4,:)*1i - M(2,3,:) - M(2,4,:)*1i;
C(4,1,:) = M(3,3,:) - M(3,4,:)*1i + M(4,3,:)*1i + M(4,4,:);
C(4,2,:) = M(3,1,:) - M(3,2,:) + M(4,1,:)*1i - M(4,2,:)*1i;
C(4,3,:) = M(1,3,:) - M(1,4,:)*1i - M(2,3,:) + M(2,4,:)*1i;
C(4,4,:) = M(1,1,:) - M(1,2,:) - M(2,1,:) + M(2,2,:);
C = C./2;
end

function M = Cov2M(C) % Cloude covariance to Mueller
M(1,1,:) = C(1,1,:) + C(2,2,:) + C(3,3,:) + C(4,4,:);
M(1,2,:) = C(1,1,:) - C(2,2,:) + C(3,3,:) - C(4,4,:);
M(1,3,:) = C(1,2,:) + C(2,1,:) + C(3,4,:) + C(4,3,:);
M(1,4,:) = ( -C(1,2,:) + C(2,1,:) - C(3,4,:) + C(4,3,:) )*1i;
M(2,1,:) = C(1,1,:) + C(2,2,:) - C(3,3,:) - C(4,4,:);
M(2,2,:) = C(1,1,:) - C(2,2,:) - C(3,3,:) + C(4,4,:);
M(2,3,:) = C(1,2,:) + C(2,1,:) - C(3,4,:) - C(4,3,:);
M(2,4,:) = ( -C(1,2,:) + C(2,1,:) + C(3,4,:) - C(4,3,:) )*1i;
M(3,1,:) = C(1,3,:) + C(2,4,:) + C(3,1,:) + C(4,2,:);
M(3,2,:) = C(1,3,:) - C(2,4,:) + C(3,1,:) - C(4,2,:);
M(3,3,:) = C(1,4,:) + C(2,3,:) + C(3,2,:) + C(4,1,:);
M(3,4,:) = ( -C(1,4,:) + C(2,3,:) - C(3,2,:) + C(4,1,:) )*1i;
M(4,1,:) = ( C(1,3,:) + C(2,4,:) - C(3,1,:) - C(4,2,:) )*1i;
M(4,2,:) = ( C(1,3,:) - C(2,4,:) - C(3,1,:) + C(4,2,:) )*1i;
M(4,3,:) = ( C(1,4,:) + C(2,3,:) - C(3,2,:) - C(4,1,:) )*1i;
M(4,4,:) = C(1,4,:) - C(2,3,:) - C(3,2,:) + C(4,1,:);
M = real(M)./2;
end

function J = nearestJones(M)
C = M2Cov(M);
J = zeros(2,2,size(C,3));
for n=1:size(C,3)
    [V,D] = eig(C(:,:,n),'vector');
    [~,mx] = max(D);
    J(:,:,n) = sqrt(D(mx))*reshape(V(:,mx),2,2).';
end
end

function Mfiltered = filterM(M)  % M to nearest physical M
C_raw = M2Cov(M);
C = zeros(size(C_raw));
for n=1:size(C_raw,3)
    [V,D] = eig(C_raw(:,:,n),'vector');
    list = find(D > 0.00001).';
    idx = 0;
    temp = zeros(4,4,length(list));
    for j = list
        idx = idx + 1;
        temp(:,:,idx) = D(j)*V(:,j)*V(:,j)';
    end
    C(:,:,n) = sum(temp,3);
end
Mfiltered = Cov2M(C);
end

function O = jonesAnisotropy(J)
K = ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(-1/2);
T = acos( K.*( J(1,1,:) + J(2,2,:) )./2);
O = (T.*K)./(sin(T));
end
