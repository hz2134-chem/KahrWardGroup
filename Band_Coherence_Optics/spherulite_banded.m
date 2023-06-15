function [DI,Zdata,Zdatadiff] = spherulite_banded(nFiber,psi,n_vect,d,Lam)

% example script:[DI,Zdata,Zdatadiff] = spherulite_banded(10,13,[1.55,1.656,1.452],150,550);

% nFiber = number of fibers
% psi = total maximum misalignment angle (degrees)
% n_vect = refractive index vector [n_x, n_y, n_z] of lamella corresponding to a twist angle of zero
% ***In this arbitrary system, n_x is the radial direction, n_z is the direction out of paper.***
% d = thickness of one lamellae, in nm
% Lam = wavelength of light, in nm

% Input material properties
psi = psi*pi/180;
epsilon_c = diag(n_vect.^2); 

% Input image, convert to binary image, read positions of endpoints of edge-on area
img = imread('band_1350-1500_sq_white.png');
[height, ~ , ~] = size(img);
N = height; %reading dimension,number of pixels along x and y, need to be *square*
binary_img = im2bw(img, graythresh(img));
%inverted_binary_img = ~binary_img; %For the case where color scale is reversed
[rows, cols] = find(binary_img);
black_pixels_coords = [cols, rows]; %Reversed x,y for convenience, can be ignored if no need

p = 280; % pitch size, in pixel, dimensionless, bassed on image input

% If incoherence is uniform across the bands
UniInc = 0.187;

% If incoherence is changing along the image: read incoherence information of each band from Excel
%data = xlsread('data.xlsx', 'Sheet1', 'A2:B9'); % read the data from the Excel file
%max_i = max(data(:, 1)); % find the maximum value of i in the data
%Inc = zeros(max_i, 1); % create a vector to store the j values
%for i = 1:max_i % loop through all values of i
  %  row_index = find(data(:, 1) == i, 1, 'first'); % find the index of the row where i equals i
  %  Inc(i) = data(row_index, 2); % extract the j value from the row where i equals i and store it in H
%end

% General idea: to modify the position of LB and misalignment

for X = 1:N
    found = false;
    NBand = 0;  %Band number, related to incoherence by yyf's simulation

    for Y = 1:N
        phi = -pi/4;
        
        if (ismember([X, Y], black_pixels_coords, 'rows'))
            found = true;
            x = 0;
            phi = 2*pi*(x)/p;
            NBand = NBand + 1;
            
        elseif found
            x = x + 1;
            phi = 2*pi*(x)/p;
        end
        %If polar coordinate
        %[theta,r] = cart2pol(y,x); 
        %theta = -theta;
        
        M = eye(4);
        
        for k = 1:nFiber
                
        % Generate normal distribution
%        if max_i - NBand == 0
%            random_INC = 0;
%        else
        sigma = UniInc * p/2;    
        %sigma = Inc(max_i - NBand) * p/2; 
        random_numbers = sigma * randn(1, 50);
        random_index = randi(length(random_numbers));
        random_INC = random_numbers(random_index);
        positive_or_negative = randi([1,2]);
            if positive_or_negative == 2
            random_INC = -random_INC;
            end
%        end
            
            t = sin(2*(phi))*psi*((k-1)/(nFiber-1)-(1/2)); %average splay angle, related to misalignment
            %R = rotz(t)*rotx(phi + pi/4); %If no incoherence
            R = rotz(t)*rotx(phi + pi/4 + 2*pi*(random_INC)/p);  
            %Whether +pi/4 or not depends on experimental observations of the crystal system; 
            %For coumarin, yes; for mannitol, no; 
            epsilon = R*epsilon_c*R.';
            L  = (sqrt(epsilon(1,1))-sqrt(epsilon(2,2)))*2*pi*d/Lam;
            Lp = (sqrt((epsilon(1,1)+epsilon(1,2)+epsilon(2,1)+epsilon(2,2))./2)-...
                sqrt((epsilon(1,1)-epsilon(1,2)-epsilon(2,1)+epsilon(2,2))./2))*2*pi*d/Lam; 
            LB = real(L);
            LD = imag(L);
            LBp = real(Lp);
            LDp = imag(Lp);
            Mat = exp(-sqrt(LD^2+LDp^2))*expm([0,-LD,-LDp,0;-LD,0,0,LBp;-LDp,0,0,-LB;0,-LBp,LB,0]); % calculating Mueller, no change
            M = Mat * M; 
        end
        %R = MMrot(theta);
        Zdata(:,:,X,Y) = M; % Mueller matrix, output
        Zdatadiff(:,:,X,Y)= real(logm(Zdata(:,:,X,Y)));
        DI(:,X,Y) = sqrt(sum(sum((Zdata(:,:,X,Y)./Zdata(1,1,X,Y)).^2))-1)/sqrt(3);
    end
end



 %   function R = MMrot(t)
 %       R =
 %       [1,0,0,0;0,cos(2*t),-sin(2*t),0;0,sin(2*t),cos(2*t),0;0,0,0,1]; %Rotated matrix, only applicable in polar coordinates
  %  end

    function matrix = rotx(angle)
        matrix = [1,0,0;0,cos(angle),-sin(angle);0,sin(angle),cos(angle)]; % rotation matrix on x
    end

    function matrix = rotz(angle)
        matrix = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1]; % rotation matrix on z
    end

% Plotting
scrsz = get(0,'screensize');
ratio = 1.4;
if 0.8*scrsz(3)/ratio > scrsz(4)
    figPos = [50 50 ratio*scrsz(4)*0.8 scrsz(4)*0.8];
else
    figPos = [50 50 scrsz(3)*0.8 scrsz(3)*0.8/ratio];
end
figure('position',figPos)
pHandles=zeros(16,1);
for j = 1:4
    for k = 1:4
        vect=[(8+(k-1)*0.25*figPos(3))/figPos(3), (8+(4-j)*0.25*figPos(4))/figPos(4), 0.23, 0.23];
        hand = subplot('Position',vect);
        if max(max(Zdata(j,k,:,:))) - min(min(Zdata(j,k,:,:))) < 0.0001
            clim = [-1 1];
        else
            clim = [min(min(Zdata(j,k,:,:))),max(max(Zdata(j,k,:,:)))];
        end
        imagesc(squeeze(Zdata(j,k,:,:)),clim)
        colormap(hand,makeColormap(clim,true,'HotCold Bright'))
        set(hand,'nextplot','replacechildren');
        axis('off')
        colorbar
        pHandles(j+4*(k-1)) = hand;
    end
end

end
