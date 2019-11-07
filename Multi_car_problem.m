% This script specifies the planning problem
nstep = 30; % number of sampling points on the path
dim = 4; % dimention of the problem, default 2
dt = 1; % sampling time (s)
start_1 = [0;0]; % starting point
ending_1 = [20;20]; % ending point
start_2 = [0;20];
ending_2 = [20;0];
margin = 0.25; % safety margin

I_2 = [eye(2), -eye(2);
    -eye(2), eye(2)];

%% Initialize the path by linear interpolation
path_1 = [];
path_2 = [];
refpath = [];
for i = 1:nstep
    path_1(:,i) = (nstep-i)/(nstep-1)*start_1+(i-1)/(nstep-1)*ending_1;
    path_2(:,i) = (nstep-i)/(nstep-1)*start_2+(i-1)/(nstep-1)*ending_2;
    refpath = [refpath; path_1(:, i); path_2(:, i)];
end

oripath = refpath;
%% Initialize the obstacles
nobj = 5;
obs = {};

obs{1}.v = [0;0]; % velocity of the obstacle
obs{1}.poly = [2 4.3 4 3;3 3 -1 -1]; % shape of the obstacle
obs{2}.v = [0;0];
obs{2}.poly = [6 7 7.5 6;0 1 -3 -3];
obs{3}.v = [0;0];
obs{3}.poly = [1 2 1.5 0;-2 -3 -0.1 -1];
obs{4}.v = [0;0];
obs{4}.poly = [1 6 1 5; 1 1 2 2];
obs{5}.v = [0;0];
obs{5}.poly = [3 7 7 3;-3 -3 -2.1 -2.1];

%% The cost function
% The distance metric between the original path and the new path
% The dimension of X would be [4*n, 1]. At each step, the x1, y1 and x2, y2
% form a 4x1 block, the whole vector X are constructed by n blocks.

% dim_low = dim / 2; %Using 2-dimension version matrix, stack them together to form the 4-dimension version matrix.
% Q1 = eye(nstep*dim_low);
% % The velocity
% Vdiff = eye(nstep*dim_low)-diag(ones(1,(nstep-1)*dim_low),dim_low);
% Q2 = Vdiff(1:(nstep-1)*dim_low,:)'*Q1(1+dim_low:end,1+dim_low:end)*Vdiff(1:(nstep-1)*dim_low,:);
% % The accelaration
% Adiff = Vdiff-diag(ones(1,(nstep-1)*dim_low),dim_low)+diag(ones(1,(nstep-2)*dim_low),dim_low*2);
% Q3 = Adiff(1:(nstep-2)*dim_low,:)'*Adiff(1:(nstep-2)*dim_low,:);
% % The weight
% cref = [0,0,0];
% cabs = [0,0,1];
% % The total costj
% Qref_temp = Q1*cref(1)/nstep+Q2*cref(2)*(nstep-1)+Q3*cref(3)*(nstep-1)^4/(nstep-2);
% Qabs_temp = Q1*cabs(1)/nstep+Q2*cabs(2)*(nstep-1)+Q3*cabs(3)*(nstep-1)^4/(nstep-2);
% 
% Qref = [Qref_temp, zeros(size(Qref_temp));
%     zeros(size(Qref_temp)), Qref_temp];
% Qabs = [Qabs_temp, zeros(size(Qabs_temp));
%     zeros(size(Qabs_temp)), Qabs_temp];

%% The cost function - version 2
% The distance metric between the original path and the new path
% The dimension of X would be [4*n, 1]. At each step, the x1, y1 and x2, y2
% form a 4x1 block, the whole vector X are constructed by n blocks.

dim = dim; %Using 2-dimension version matrix, stack them together to form the 4-dimension version matrix.
Q1 = eye(nstep*dim);
% The velocity
Vdiff = eye(nstep*dim)-diag(ones(1,(nstep-1)*dim),dim);
Q2 = Vdiff(1:(nstep-1)*dim,:)'*Vdiff(1:(nstep-1)*dim,:);
% The accelaration
Adiff = Vdiff-diag(ones(1,(nstep-1)*dim),dim)+diag(ones(1,(nstep-2)*dim),dim*2);
Q3 = Adiff(1:(nstep-2)*dim,:)'*Adiff(1:(nstep-2)*dim,:);
% The weight
cref = [0.05,2,0];
cabs = [0,0,10];
% The total costj
% Qref = Q1*cref(1)/nstep+Q2*cref(2)*(nstep-1)+Q3*cref(3)*(nstep-1)^4/(nstep-2);
% Qabs = Q1*cabs(1)/nstep+Q2*cabs(2)*(nstep-1)+Q3*cabs(3)*(nstep-1)^4/(nstep-2);
Qref = Q1*cref(1)+Q2*cref(2)+Q3*cref(3);
Qabs = Q1*cabs(1)+Q2*cabs(2)+Q3*cabs(3);

%% The boundary constraint
% Making sure that the start point and the end point are fixed.

% Aeq = zeros(2*dim,nstep*dim+nstep);
% Aeq(0*dim+1:1*dim,1:dim) = eye(dim);
% Aeq(1*dim+1:2*dim,(nstep-1)*dim+1:nstep*dim) = eye(dim);
% Aeq(2*dim+1:3*dim,1:2*dim) = [-eye(dim) eye(dim)];
% Aeq(3*dim+1:4*dim,(nstep-2)*dim+1:nstep*dim) = [-eye(dim) eye(dim)];
% beq = [path(:,1);path(:,end);path(:,2)-path(:,1);path(:,end)-path(:,end-1)];
Aeq = zeros(2*dim,nstep*dim);
Aeq(0*dim+1:1*dim,1:dim) = eye(dim);
Aeq(1*dim+1:2*dim,(nstep-1)*dim+1:nstep*dim) = eye(dim);
beq = [oripath(1:4);oripath(end-3:end)];

