
%% helper functions
%new_dist_angle = @(pt, a, l) pt + [cos(a)*l; sin(a)*l; 0];
euclid_len = @(list_vec) sqrt(sum(list_vec.^2,2));
euclid_dist = @(s_pts,t_pts) sqrt(sum((s_pts-t_pts).^2,2)); %euclid_dist = @(list_pts,s,t) sqrt(sum((list_pts(s,:)-list_pts(t,:)).^2,2)); %default use Euclidean distance, but could modfiy distance formula to adjust the initial "lengths" e.g. for dynamic relaxation 
custom_dist = euclid_dist;
%pt_index = @(i,j,N) 1+(i-1)*N+j; 

DOS = @(std,mean,M) std.^2 + (mean-10./M);
deform = @(xarray_i, xarray_f) sqrt(sum((xarray_i-xarray_f).^2,2)); %sum along rows for each row = an edge, each column = a coordinate dimension


%%

nvec = 3:4%10;
mvec = 1:2%5; 

DOS_initlens = zeros(length(mvec),length(nvec));
DOS_finlens = zeros(length(mvec),length(nvec));
DOS_angle = zeros(length(mvec),length(nvec));

for n = nvec
    for m = mvec
        export = false;
        %[DOS_initlens(m,n),DOS_finlens(m,n),DOS_angle(m,n)]=final_fdm(m,n,'uniform',true);
        [DOS_initlens(m,n),DOS_finlens(m,n),DOS_angle(m,n)]=final_fdm_program(m,n,'uniform',false);
    end
end
%%
fileout='initial\gridsample_DOSlenresults.csv';
csvwrite(fileout,[DOS_initlens, DOS_finlens]);
fileout='initial\gridsample_DOSangleresults.csv';
csvwrite(fileout,DOS_angle);
%%
figure
plot(DOS_initlens(mvec,nvec),'ok'); hold on;
figure
plot(DOS_finlens(mvec,nvec),'ok'); hold on;
figure
plot(DOS_angle(mvec,nvec),'ok'); hold on;

%%
    M = 3; N = 30; type = 'uniform';
    [dv, coords_i, list_edges, pts] = generate_initial_geom(M,N,type,false);
%%
    
    num_pts = length(coords_i);%1+ N*M;
    num_edges = length(list_edges);%2*N*M;
    
    q = ones(num_edges, 1)*100000;
    p = [0, -10, -1];
    [dv, coords, f, l] = fdm_plot(dv,coords_i,list_edges,q,p,true);
    
%%
    figure
    plot(sort(f))
    figure
    hist(f)
