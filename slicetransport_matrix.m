% Shiying Li and Caroline Moosmueller, 2024. 
% stochastsic slice transport - matrix version (using Radon slices at theta and its orthogonal angle) 

function [I,tp0,R0e1,R0e2,R1e1,R1e2] =slicetransport_matrix(I0,I1,theta,gamma)


%%% input
% gamma - step size 
% theta - an angle between 0 to 89 
% I0, I1 - input images - now only for square images 
% I0 is mu and I1 is sigma in the paper

% compute slice projections from I0 -->I1 to form h then morph I1 accordly 
%I =|detJ_h| I_1(h)


I0 = I0./sum(sum(I0));  I1 = I1./sum(sum(I1)); %  normalize

theta_perp = theta + 90; % orthogonal angle
% Radon transform in two direction
[R0,tp0] = radon(I0,[theta,theta_perp]);   % Radon slices on associated directions
[R1,~] = radon(I1,[theta,theta_perp]); % assume I0 and I1 are of the same size

R0e1 = R0(:,1)'; R0e2 = R0(:,2)'; % row vector, two slices for I0
R1e1 = R1(:,1)'; R1e2 = R1(:,2)';

epsR =(1e-10)* max([max(R0e1),max(R0e2),max(R1e1),max(R1e2)]);
epsR = max(eps, epsR);
R0e1 = R0e1 + epsR; R0e1 = R0e1/sum(R0e1); % make positive and normalized
R0e2 = R0e2 + epsR; R0e2 = R0e2/sum(R0e2); 
R1e1 = R1e1 + epsR; R1e1 = R1e1/sum(R1e1);
R1e2 = R1e2 + epsR; R1e2 = R1e2/sum(R1e2);

%%% caculate 1D optimal transport between slices 

rad = theta/180*pi; 
e1 = [cos(rad); sin(rad)]; % unit vector whose angle is theta
e2 = [-sin(rad); cos(rad)]; % orthorgonal vector

S0e1 = cumsum(R0e1); % cumulative sum of R0e1
S0e2 = cumsum(R0e2);

S1e1 = cumsum(R1e1); % cumulative sum of R1e1
S1e2 = cumsum(R1e2);



len0 = size(I0,1);
hlen0 = len0/2;
x_splim1 = -hlen0; x_splim2 = hlen0;
y_splim1 = x_splim1; y_splim2 = x_splim2;

lr0 = length(tp0);
xi = linspace(0,1,lr0); % interpolation points to S^{-1}

% OT for e1

S0e1_inv = interp1(S0e1, tp0, xi,'pchip');
tau_e1 = interp1(xi(2:end-1),S0e1_inv(2:end-1),S1e1,'pchip'); % tau_e1(tp0) = S0^{-1}_e1 \circ S1e1 (tp0) - ot from I1e1 to I0e1

tau_e1_gamma = (1-gamma)*tp0' + gamma*tau_e1; % stochastic combination

% OT for e2

S0e2_inv = interp1(S0e2,tp0,xi,'pchip');
tau_e2 = interp1(xi(2:end-1),S0e2_inv(2:end-1),S1e2,'pchip'); % tau_e2 (tp0) = S0^{-1}_e2 \circ S1e2 (tp0) - ot from I1e2 to I0e2
tau_e2_gamma = (1-gamma)*tp0' + gamma*tau_e2; % stochastic combination

% inverse of one-dimensinoal maps
g1 = interp1(tau_e1_gamma,tp0,tp0,'linear');
g2 = interp1(tau_e2_gamma,tp0,tp0,'linear');

inan = find(isnan(g1));
ivals = setdiff(1:length(g1),inan);
inanl = find(inan<ivals(1));
inanr = find(inan>ivals(end));
il = inan(inanl);
ir = inan(inanr);
g1(il) = g1(ivals(1))*ones(length(il),1);
g1(ir)= g1(ivals(end))*ones(length(ir),1);

inan = find(isnan(g2));
ivals = setdiff(1:length(g2),inan);
inanl = find(inan<ivals(1));
inanr = find(inan>ivals(end));
il = inan(inanl);
ir = inan(inanr);
g2(il) = g2(ivals(1))*ones(length(il),1);
g2(ir)= g2(ivals(end))*ones(length(ir),1);

% deratives
tp0_step = tp0(2)-tp0(1);
der_g1= gradient(g1,tp0_step);
der_g2 = gradient(g2,tp0_step);


%%% Construct 2D vector field

xg = -hlen0+.5+[0:len0-1];
yg =flip(xg); % grid evaluation points -pixels at center of each sub-block
[Xg,Yg] = meshgrid(xg,yg);
ng = size(Xg,1)*size(Xg,2);
Xv = reshape(Xg,1,ng); % reshape grid of x-coord in to a long row vector
Yv = reshape(Yg,1,ng); % reshape grid of y-coord in to a long row vector
XYv = [Xv;Yv]; % each column is a (x,y) coord of a point; 2 rows 
P = [e1,e2]; % orthogonal matrix

projXY = P'*XYv; % projected coord in directions e1 and e2
e1XY = projXY(1,:); % dot product of e1 with grid points, in a row vector
e2XY = projXY(2,:); % dot product of e2 with grid points, in a row vector

T = P*[interp1(tp0,g1,e1XY); interp1(tp0,g2,e2XY)]; % vector field T on gird points, 2 by n matrix, where n is # grid pts
% T(x) = P[tau_e1(e1^Tx); tau_e2(e2^Tx)]
J = interp1(tp0,der_g1,e1XY).*interp1(tp0,der_g2,e2XY); % Jacobian of T at a row vector of n grid pts

% displacement
% transform coordinates from unit square to image intrinsic coordinates
Crd1 = imref2d(size(I1),[x_splim1,x_splim2],[y_splim1,y_splim2]);
[T1Intrin, T2Intrin] = worldToIntrinsic(Crd1,T(1,:),T(2,:));
[XvIntrin, YvIntrin] = worldToIntrinsic(Crd1,Xv,Yv);
[ny,nx] = size(I1);
D = zeros(ny,nx,2); 
D(:,:,1) = reshape(T1Intrin-XvIntrin,ny,nx); % T1Intrin is a row vector, need to reshape to grid format
D(:,:,2) = -reshape(T2Intrin-YvIntrin,ny,nx);  % note intrinsic y-axis is opposite to the radon y-axis


% warp I1
JT = reshape(J,ny,nx); % reshape back to grid
I = JT.*imwarp(I1,D);












