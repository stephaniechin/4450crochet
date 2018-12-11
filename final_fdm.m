function [DOS_initlen,DOS_finlen,DOS_angle] = final_fdm(M,N,type,export) 
close all
%clear variables

%% helper functions
new_dist_angle = @(pt, a, l) pt + [cos(a)*l; sin(a)*l; 0];
euclid_len = @(list_vec) sqrt(sum(list_vec.^2,2));
euclid_dist = @(s_pts,t_pts) sqrt(sum((s_pts-t_pts).^2,2)); %euclid_dist = @(list_pts,s,t) sqrt(sum((list_pts(s,:)-list_pts(t,:)).^2,2)); %default use Euclidean distance, but could modfiy distance formula to adjust the initial "lengths" e.g. for dynamic relaxation 
custom_dist = euclid_dist;
pt_index = @(i,j,N) 1+(i-1)*N+j; 

DOS = @(std,mean,u) std.^2 + (mean-u);
deform = @(xarray_i, xarray_f) sqrt(sum((xarray_i-xarray_f).^2,2)); %sum along rows for each row = an edge, each column = a coordinate dimension

%% constants and design variables 
A = 2*pi; %constant
L = 10; %constant
%M = 1; %design variable
%N = 3; %design variable %% need to update calculation of pt_index expression if change N!

if strcmp(type,'uniform')
    filepre = strcat('uniform1000\uniform',num2str(M),num2str(N))
    depths = L/M * ones(1,M);
    angles = A/N * ones(1,N);
elseif strcmp(type,'random')
    filepre = strcat('initial\random',num2str(M),num2str(N))
    s = 0; rng(s);%randomization seed
    j_ = rand(1,M);
    depths = L/sum(j_)*j_
    j_ = rand(1,N);
    angles = A/sum(j_)*j_ 
elseif strcmp(type,'randomM')
    filepre = strcat('initial\randomM',num2str(M),num2str(N))
    s = 0; rng(s);%randomization seed
    j_ = rand(1,M);
    depths = L/sum(j_)*j_
    angles = A/N * ones(1,N);
elseif strcmp(type,'randomN')
    filepre = strcat('initial\randomN',num2str(M),num2str(N))
    depths = L/M * ones(1,M);
    s = 0; rng(s);%randomization seed
    j_ = rand(1,N);
    angles = A/sum(j_)*j_ 
%else if strcmp(type,'custom') %prompt for distance and angles    
end
%{
%for evenly spaced depths and angles
uniform_depths = L/M * ones(1,M);
uniform_angles = A/N * ones(1,N); 
%for randomly spaced depths and angles
s = 0; rng(s);%randomization seed
j_ = rand(1,M);
random_depths = L/sum(j_)*j_;
j_ = rand(1,N);
random_angles = A/sum(j_)*j_; 

depths = uniform_depths; 
angles = uniform_angles; 
%}
gangles = angles; %global cumulative angles 
for j = 2:N
    gangles(j) = gangles(j) + gangles(j-1);
end
disp(gangles)
%% set up lists of points and edges

pts = zeros(3,N,M); % num points = 1 + N*M %index points as pt_index = @(i,j) 1+(i-1)*N+j; 
num_pts = 1+ N*M;
num_edges = 2*N*M;
list_pts = zeros(num_pts,3); %includes center point, which is not in pts
list_edges = zeros(num_edges,2); % num lines = 2*N*M 
%initlen = zeros(num_edges);

Cpt = [0;0;0];
list_pts(1,:) = Cpt;
%% calculate and plot points
figure(1)
plot(Cpt(1),Cpt(2),'o'); hold on; 

for i = 1:M
    for j = 1:N
        
        if i == 1
            currentpoint = Cpt;
        else 
            currentpoint = pts(:,j,i-1);
        end
        new_point = new_dist_angle(currentpoint, gangles(j), depths(i));
        scatter(new_point(1), new_point(2),...
            'or'); hold on; 

        pts(:,j,i) = new_point;  
        list_pts(pt_index(i,j,N),:) = new_point;       
    end
end

%% calculate and plot edges
edge_index = 1; 
for i = 1:M
    for j = 1:N
        %radial edges
        if i == 1
            list_edges(edge_index,:) = [1,pt_index(i,j,N)];
            s = Cpt; 
        else
            list_edges(edge_index,:) = [pt_index(i-1,j,N),pt_index(i,j,N)];
            s = pts(:,j,i-1);
        end
        t = pts(:,j,i);
        plot([s(1) ; t(1)],[s(2) ; t(2)],...
             '-','color',[0,0,edge_index/num_edges]); hold on; 
         
        %list_initlen(edge_index) = custom_dist(list_pts,s,t); 
        edge_index = edge_index + 1; 
            
        %ring edges
        if j == N; j_ =1; else; j_ = j+1; end
      
        s = pts(:,j,i);
        t = pts (:,j_,i); 
        
        plot([s(1) ; t(1)],[s(2) ; t(2)],...
             '-','color',[0,0,edge_index/num_edges]); hold on; 
         
        list_edges(edge_index,:) = [pt_index(i,j,N),pt_index(i,j_,N)];
        %list_initlen(edge_index) = custom_dist(list_pts,s,t); 
        edge_index = edge_index + 1; 
    end    
end
figure(1)
axis equal;
axis([-10 10 -10 10]);
xticks(-10:2:10);
yticks(-10:2:10);
grid on
if export
    fileout = strcat(filepre,'init.jpg');
    saveas(gcf,fileout);
end


%%
%future: format graph and color points/lines 
%%
%future: export relevant locations, distances, etc to csv or save as .dat


%% FDM

%% Variables definition
% Initial positions (only those of fixed points matter)
x0 = list_pts(:,1);
y0 = list_pts(:,2);
z0 = list_pts(:,3);


% Fixed points indices
fixed_points = zeros(N,1);
threshold = 0.1; 
k = 1; %counter
for i = 1:length(x0)
    d = sqrt(sum([x0(i) y0(i) z0(i)].^2));
    if abs(d - L) < threshold
        fixed_points(k) = i;
        k = k+1; 
    end
end

num_f = length(fixed_points); %check that this equals N (design variable)
num_n = length(x0) - num_f;


%% Connectivity matrix
C = zeros(num_edges,num_pts);

for i = 1:num_edges
    s = list_edges(i,1);
    t = list_edges(i,2); 
    
    C(i,s) = 1;
    C(i,t) = -1; 
end

% Force densities
q = ones(num_edges, 1)*1000000;

%% Uniform load 
% applied to free points
p = [0, 0, -1];

%% Solve
% Prepare intermediate matrices
new_points = setdiff(1:size(x0), fixed_points); %aka free points
Q = diag(q);
Cn = C(:, new_points); %choose columns of C associated with free points
Cf = C(:, fixed_points);%choose columns of C associated with fixed points
Dn = Cn' * Q * Cn;
Df = Cn' * Q * Cf;
% Solve
xf = x0(fixed_points);
yf = y0(fixed_points);
zf = z0(fixed_points);
xn = linsolve(Dn, p(1) - Df * xf);
yn = linsolve(Dn, p(2) - Df * yf);
zn = linsolve(Dn, p(3) - Df * zf);

%% Post-process
% Recompose coordinates
coords = zeros(size(x0, 1), 3);
coords(fixed_points, :) = [xf, yf, zf];
coords(new_points, :) = [xn, yn, zn];
% Get member forces
l = sum((C * coords).^2, 2) .^ 0.5;
f = Q * l;



%% Graph
%%{
% Loop over each member and plot it
figure(2);
set(gcf, 'color', 'w');
for i = 1:num_edges
    % Get coordinates of nodes at both ends of that member
    c = coords(C(i, :) ~= 0, :);
    % Change the settings of the plot3 command to get variations in color,
    % line thickness, etc.  Here, I am making the structure blue and
    % thickening the lines based on their relative internal forces.
    fnorm = f ./ max(f);
    plot3(c(:, 1), c(:, 2), c(:, 3), 'b', 'LineWidth', fnorm(i));
    hold all; 
end
figure(2)
%axis equal;
%axis([-10 10 -10 10]);
xticks(-10:2:10);
yticks(-10:2:10);
%{
% label final positions of all points.
labels = cell(num_pts,1);
for i = 1:num_pts
    labels{i} = mat2str(round(coords(i,:)));
end
text(coords(:,1),coords(:,2),coords(:,3),labels)
axis equal;

%}
if export
    fileout = strcat(filepre,'final.jpg');
    saveas(gcf,fileout);
end
%%
initvec = list_pts(list_edges(:,2),:) - list_pts(list_edges(:,1),:); 
finvec =  coords(list_edges(:,2),:) - coords(list_edges(:,1),:); 
diffvec = finvec-initvec;
%initlen = custom_dist(list_pts(list_edges(:,1),:),list_pts(list_edges(:,2),:))
%initlen = custom_dist(initial start points, initial end points) 
initlen = euclid_len(initvec);                        
%finlen = custom_dist(coords(list_edges(:,1),:),coords(list_edges(:,2),:))
%finlen = custom_dist(final start points, final end points) 
%finlen = euclid_len(finvec);                       
finlen = l; 
difflen = euclid_len(diffvec);                    
%difflen = custom_dist(list_initlen, list_finlen)

diff_nodal = coords - list_pts; %final - initial coordinates
dist_nodal = euclid_len(diff_nodal); 

fin_ptangle1 = atan(coords(:,2)./coords(:,1)); % - angles(j)
diff_edgeangle2 = sec((finlen.^2 + initlen.^2 - difflen.^2)...
                    ./(2*finlen.*initlen)); %based on law of cosines to find the angle between the initial vector and the final vector for each edge
finlen2D = euclid_len(finvec(:,1:2));
initlen2D = euclid_len(initvec(:,1:2));
difflen2D = euclid_len(finvec(:,1:2)-initvec(:,1:2));

diff_edgeangle3 = sec(((finlen2D).^2 + (initlen2D).^2 - (difflen2D).^2) ...
                    ./(2*finlen2D.*initlen2D)); %
%%
DOS_initlen = DOS(std(initlen),mean(initlen),L/M);
DOS_finlen = DOS(std(finlen),mean(finlen),L/M);
DOS_angle = DOS(std(fin_ptangle1),mean(fin_ptangle1),A/N);
figure(3)
[sort_initlen, I_initlen] = sort(DOS_initlen);
plot(sort_initlen)
figure(4)
[sort_finlen, I_finlen] = sort(DOS_finlen);
plot(sort_finlen)
figure(5)
[sort_angle, I_angle] = sort(DOS_angle);
plot(sort_angle)
%%
%%{
if export
    fileout = strcat(filepre,'_edges.csv')
    csvwrite(fileout,[list_edges,initvec,finvec,difflen,difflen2D,diff_edgeangle2,diff_edgeangle3])

    fileout = strcat(filepre,'_nodes.csv')
    csvwrite(fileout,[list_pts,coords,diff_nodal,dist_nodal,fin_ptangle1])
end
%}
end 