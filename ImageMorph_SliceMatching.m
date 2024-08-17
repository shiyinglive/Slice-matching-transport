% Shiying Li and Caroline Moosmueller. Measure transfer via stochastic
% slicing and matching, 2024.  

% this scripts shows how a source image is "morphed" to a target image via the
% slice-matching iterative scheme 
% Two versions, at each iteration:
%%% (a) matrix slice -  a random orthogonal matrix is used
%%% (b) single slice -  a random angle is used 


close all; clear all;

ntrials = 15; % number of different trials

%%% (a) matrix slice version
n=20; % number of iteration in each trial
% %%% (b) single slice version
% n = 30;

sw_dist_all = zeros(ntrials,n+1);

for j = 1:ntrials

vis_const = 2.2; % this constant is for enhancing visualization. 
I1 = imread('I1.png');  % source image
I0 = imread('I0.png');  % target image

% generate slicing directions for all iterations

%%% matrix slice version - choose n random angles; note their orthogonal angles will be used to formulate the orthogonal matrix at each iteration
angs = randi([0,89],1,n); % the angles along with 

% %%% single slice version - choose n random angles from [0,pi)
% angs = randi([0,179],1,n);


% normalize to make these images probability density functions (pdfs)
I0 = double(I0);   I1 = double(I1);
I0 = I0./sum(sum(I0));  I1 = I1./sum(sum(I1)); %  pdfs
epsadd = 1e-7; % ensure positivity 
I0 = I0+epsadd; I1 = I1+epsadd;
I0 =  I0./sum(sum(I0)); I1 = I1./sum(sum(I1)); %  pdfs



% initialization
sw_dist = zeros(1,n+1);
sw_dist(1) = SW(I1,I0);

%%%%%
step_size = zeros(1,n);
for i = 1:n
    % step size
    gamma = (1+log2(i))/i; % step size satisfying the theorem assumptions

    step_size(i)=gamma;
    max0 = max(max(I0));  max1 = max(max(I1)); 
    theta = angs(i); 

    %%% matrix slice version
    [I1w,tp0,R0e1,R0e2,R1e1,R1e2] = slicetransport_matrix(I0,I1,theta,gamma); % transport I1 to I1w using sliced matching in theta and its orthogonal angle
    % %%% single slice version - comment out the previous line 
    % [I1w,tp0,R0e1,R1e1] = slicetransport_theta(I0,I1,theta,gamma);

    sw_dist(i+1) = SW(I1w,I0); % distances from morphed image I1w to target I0 at all iterations at the j-th trial
    I1 = I1w; 
end
sw_dist_all(j,:)= sw_dist; % collect distances for all trials 


end 
plot_errors(sw_dist_all)



