%{
Given initial geometry (see generate_initial_geom.m) and other inputs:
INPUTS: 
    dv = cell array of relevant design variables (M,N,A,L,depths,angles,filepre)
    coords_i = 1+N*M x 3 matrix 
    list_edges = 2*N*M x 2 matrix
    q = 1+NxM vector of force densities
    p = 1x3 vector of applied force(free points)

Assumes: 
    fixed points are all points on outermost ring

Steps
- Uses FDM to calculate final geometry
- Plots final node locations in 3D (OPT: export .jpg)
- Returns final geometry, lengths, and forces

OUTPUTS: 
    dv = cell array of relevant design variables (M,N,A,L,depths,angles,filepre)
    coords_i = 1+N*M x 3 matrix 
    list_edges = 2*N*M x 2 matrix
    pts = NxMx3 matrix
%}

function [dv, coords, f, l] = fdm_plot(dv, coords_i, list_edges, q, p, export)
    M = dv{1}; N = dv{2}; L = dv{3}; A = dv{4};filepre = dv{7} %subfolder = 'uniform1000';filepre = strcat(subfolder,'\',type,num2str(M),num2str(N))

    num_pts = length(coords_i);%1+ N*M;
    num_edges = length(list_edges);%2*N*M;

    
    
    %% Variables definition
    % Initial positions (only those of fixed points matter)
    x0 = coords_i(:,1);
    y0 = coords_i(:,2);
    z0 = coords_i(:,3);


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
    %q = ones(num_edges, 1)*100000;

    %% Uniform load 
    % applied to free points
    %p = [0, -10, -1];

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
    figure
    xticks(-10:2:10);
    yticks(-10:2:10);
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



end
