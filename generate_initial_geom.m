%{
Given initial design variables:
INPUTS: 
    M = number of rings
    N = number of rays

Assumes following design constants:
    A = 2pi
    L = 10
    default type is uniform (unless otherwise specified)

Steps
- Calculates spacing of rings and rays
- Generates initial geometry 
- Plots initial node locations (OPT: export .jpg)
- Returns initial geometry

OUTPUTS: 
    dv = cell array of relevant design variables (M,N,A,L,depths,angles,filepre)
    coords_i = 1+N*M x 3 matrix 
    list_edges = 2*N*M x 2 matrix
    pts = NxMx3 matrix
%}

function [dv, coords_i, list_edges, pts] = generate_initial_geom(M,N,type,export)
    
    %% constants and design variables 
    A = 2*pi; %constant
    L = 10; %constant
    %M = 1; %design variable
    %N = 3; %design variable 

    % future: switch to handle optional input arguments
    subfolder = 'uniform1000';
    filepre = strcat(subfolder,'\',type,num2str(M),num2str(N))

    %% spacing between rings and rays
    switch type
        case 'uniform'
            depths = L/M * ones(1,M);
            angles = A/N * ones(1,N);
        case 'random'
            s = 0; rng(s);%randomization seed
            j_ = rand(1,M);
            depths = L/sum(j_)*j_
            j_ = rand(1,N);
            angles = A/sum(j_)*j_     
        case 'randomM'
            s = 0; rng(s); j_ = rand(1,M);%randomization seed
            depths = L/sum(j_)*j_
            angles = A/N * ones(1,N);
        case 'randomN'    
            depths = L/M * ones(1,M);
            s = 0; rng(s); j_ = rand(1,N);%randomization seed
            angles = A/sum(j_)*j_ 
        %case 'custom' 
        otherwise
            disp('default: uniform')
            depths = L/M * ones(1,M);
            angles = A/N * ones(1,N);        
    end

    gangles = angles; %global cumulative angles 
    for j = 2:N
        gangles(j) = gangles(j) + gangles(j-1);
    end
    disp(gangles)
    %% set up lists of points and edges

    pts = zeros(3,N,M); % num points = 1 + N*M %index points as pt_index = @(i,j) 1+(i-1)*N+j; 
    num_pts = 1+ N*M;
    num_edges = 2*N*M;
    coords_i = zeros(num_pts,3); %includes center point, which is not in pts
    list_edges = zeros(num_edges,2); % num lines = 2*N*M 
    %initlen = zeros(num_edges);

    Cpt = [0;0;0];
    coords_i(1,:) = Cpt;
    %% calculate and plot points
    figure
    plot(Cpt(1),Cpt(2),'o'); hold on; 

    axis equal;
    axis([-10 10 -10 10]);
    xticks(-10:2:10);
    yticks(-10:2:10);
    grid on
    
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
            coords_i(pt_index(i,j,N),:) = new_point;       
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

            %list_initlen(edge_index) = custom_dist(coords_i,s,t); 
            edge_index = edge_index + 1; 

            %ring edges
            if j == N; j_ =1; else; j_ = j+1; end

            s = pts(:,j,i);
            t = pts (:,j_,i); 

            plot([s(1) ; t(1)],[s(2) ; t(2)],...
                 '-','color',[0,0,edge_index/num_edges]); hold on; 

            list_edges(edge_index,:) = [pt_index(i,j,N),pt_index(i,j_,N)];
            %list_initlen(edge_index) = custom_dist(coords_i,s,t); 
            edge_index = edge_index + 1; 
        end    
    end
    if export
        fileout = strcat(filepre,'init.jpg');
        saveas(gcf,fileout);
    end

    %export design variables, in case needed in the future 
    dv = {M,N,L,A,depths,angles,filepre};
end

