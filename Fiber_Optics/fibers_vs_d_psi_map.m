function fibers_vs_d_psi_map(nFiber, n_vect, Lam)
    psi = linspace(0, 180, 91); % psi in degrees
    d = linspace(0, 3000, 301);
    cbValues = zeros(length(psi), length(d));
    lbValues = zeros(length(psi), length(d));  
    epsilon_c = diag(n_vect.^2);

    for i = 1:length(d)
        for j = 1:length(psi)
            psi_rad = psi(j) * pi / 180; % psi converted to rad
            M = eye(4);

            for k = 1:nFiber
                t = psi_rad * ((k - 1) / (nFiber - 1)); % average splay angle
                R = rotz(t) * rotx(0); % rotation matrix
                epsilon = R * epsilon_c * R';
                L = (sqrt(epsilon(1,1)) - sqrt(epsilon(2,2))) * 2 * pi * d(i) / Lam;
                Lp = (sqrt((epsilon(1,1) + epsilon(1,2) + epsilon(2,1) + epsilon(2,2)) / 2) - ...
                     sqrt((epsilon(1,1) - epsilon(1,2) - epsilon(2,1) + epsilon(2,2)) / 2)) * 2 * pi * d(i) / Lam;
                LB = real(L);
                LD = imag(L);
                LBp = real(Lp);
                LDp = imag(Lp);
                Mat = exp(-sqrt(LD^2 + LDp^2)) * expm([0, -LD, -LDp, 0; -LD, 0, 0, LBp; -LDp, 0, 0, -LB; 0, -LBp, LB, 0]);
                M = Mat * M;
            end            
            J = nearestJones(M);
            O = jonesAnisotropy(J); 
            LB_M = real(1i * O * (J(1,1) - J(2,2)));
            LBp_M = real(1i * O * (J(1,2) + J(2,1)));
            lb = sqrt(LB_M^2 + LBp_M^2); % linear retardance
            cb = real(O * (J(1,2) - J(2,1))); % circular retardance

            lbValues(j, i) = lb;
            cbValues(j, i) = cb;
        end
    end
    
% \\ FIGURE PLOTTING \\
blueWhiteRed = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
                ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
blackToGreen = [zeros(256,1), linspace(0,1,256)', zeros(256,1)];
[PsiGrid, DGrid] = meshgrid(psi, d);
figure;

% circular retardance
subplot(2, 1, 1);
surf(PsiGrid, DGrid, cbValues', 'EdgeColor', 'none');
title('CR', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 22);
xlabel('Psi (°)', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 22);
ylabel('d (\mum)', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 22);
zlabel('Circular Birefringence', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 22);
colormap(blueWhiteRed);
caxis([-0.5 0.5]);
colorbar;
set(gca, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 22);

% linear retardance
subplot(2, 1, 2);
surf(PsiGrid, DGrid, lbValues', 'EdgeColor', 'none');
title('|LR|', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 22);
xlabel('Psi (°)', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 22);
ylabel('d (\mum)', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 22);
zlabel('Linear Birefringence', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 22);
colormap(gca, blackToGreen);
caxis([0 3.14]);
colorbar;
set(gca, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 22);

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