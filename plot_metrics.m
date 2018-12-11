%{
Given initial and final geometry and performance(see generate_initial_geom.m, fdm_plot.m):
INPUTS: 
    dv = current design parameters
    coords_i = 1+N*M x 3 matrix 
    list_edges = 2*N*M x 2 matrix
    q = 1+NxM vector of force densities
    p = 1x3 vector of applied force(free points)

Steps
- Calculates various performance metrics
- Plots final node locations in 3D (OPT: export .jpg)
- Returns final geometry, lengths, and forces

OUTPUTS: 
    dv = cell array of relevant design variables (M,N,A,L,depths,angles,filepre)
    pm = cell array of performance metrics
    DOS
%}


function [dv, pm, DOS_init,DOS_fin,DOS_fina] = plot_metrics(dv, coords_i,coords,list_edges,q,p,f,l) 

    %% helper functions
    euclid_len = @(list_vec) sqrt(sum(list_vec.^2,2));
    euclid_dist = @(s_pts,t_pts) sqrt(sum((s_pts-t_pts).^2,2)); %euclid_dist = @(list_pts,s,t) sqrt(sum((list_pts(s,:)-list_pts(t,:)).^2,2)); %default use Euclidean distance, but could modfiy distance formula to adjust the initial "lengths" e.g. for dynamic relaxation 
    custom_dist = euclid_dist;
    
    DOS = @(std,mean,u) std.^2 + (mean-u);
    deform = @(xarray_i, xarray_f) sqrt(sum((xarray_i-xarray_f).^2,2)); %sum along rows for each row = an edge, each column = a coordinate dimension


    %% import design values
    M = dv{1}; 
    N = dv{2}; 
    L = dv{3}; 
    A = dv{4};
    filepre = dv{7} %subfolder = 'uniform1000';filepre = strcat(subfolder,'\',type,num2str(M),num2str(N))

    %% plot distributions 
    
    %% calculate final performance metrics
    T = sum(q.*l);  
    T2 = sum(f);
    
    rho_c = M/L;
    rho_l = N/A; 

    
    
    initvec = coords_i(list_edges(:,2),:) - coords_i(list_edges(:,1),:); 
    finvec =  coords(list_edges(:,2),:) - coords(list_edges(:,1),:); 
    diffvec = finvec-initvec;
    initlen = euclid_len(initvec);    
    
    diff_nodal = coords - coords_i; %final - initial coordinates
    dist_nodal = euclid_len(diff_nodal); 

    fin_ptangle1 = atan(coords(:,2)./coords(:,1)); % - angles(j)
    diff_edgeangle2 = sec((l.^2 + initlen.^2 - difflen.^2)...
                        ./(2*l.*initlen)); %based on law of cosines to find the angle between the initial vector and the final vector for each edge
    finlen2D = euclid_len(finvec(:,1:2));
    initlen2D = euclid_len(initvec(:,1:2));
    difflen2D = euclid_len(finvec(:,1:2)-initvec(:,1:2));

    diff_edgeangle3 = sec(((finlen2D).^2 + (initlen2D).^2 - (difflen2D).^2) ...
                        ./(2*finlen2D.*initlen2D)); %
                
    [sort_init, I_init] = sort(initlen);
    [sort_fin, I_fin] = sort(l);
    [sort_fina, I_fina] = sort(fin_ptangle1);
    
    DOS_init = DOS(std(initlen),mean(initlen),M);
    DOS_fin = DOS(std(l),mean(l),M);
    DOS_fina = DOS(std(fin_ptangle1),mean(fin_ptangle1),M);
    
    

    pm = {T, rho_c, rho_l, DOS_init, DOS_fin, DOS_fina, delta, F}; %could also express as vector since these should all be single values
    
    
end
