clear;
clc;
close all;

path = "\\storage3.ad.scilifelab.se\monalisa\Big_MoNaLISA\2021project\Project\characterization\intensity\405\";
name = "405_filtered_projected.tif";

im = double(imread(path+name));
px = 100;
%%
block = 50;
I = floor((size(im,1) - block)/block);
J = floor((size(im,2) - block)/block);
k = 1;
nx = zeros(1, (I+1)*(J+1));
ny = zeros(1, (I+1)*(J+1));
p = zeros(1, (I+1)*(J+1));
rs = zeros(1, (I+1)*(J+1));

for i = 0:I
    for j = 0:J
    range_x = 1+block*i:block+block*i;
    range_y = 1+block*j:block+block*j;
    crop = im(range_x, range_y);
    proj_x = squeeze(mean(crop, 1));
    x = 0:size(proj_x,2)-1;
    x = x*px;
    proj_x = (proj_x - min(proj_x(:)))./(max(proj_x(:))-min(proj_x(:))) - 1/2;
    [fit, r] = twosin(x, proj_x);
    p_x = 2*pi/fit.b1;
    rs_x = r.rsquare;

    proj_y = squeeze(mean(crop, 2));
    proj_y = proj_y';
    y = 0:size(proj_x,2)-1;
    y = y*px;
    proj_y = (proj_y - min(proj_y(:)))./(max(proj_y(:))-min(proj_y(:))) - 1/2;
    [fit, r] = twosin(y, proj_y);
    p_y = 2*pi/fit.b1;
    rs_y = r.rsquare;

    nx(k) = i;
    ny(k) = j;
    p(k) = (p_x+p_y)/2;
    rs(k) = (rs_x+rs_y)/2;
    k = k+1;
    end
end
%%
lim = 0.7;
scatter((nx(rs>lim)*block+block)/10, (ny(rs>lim)*block+block)/10, 1000, p(rs>lim), '.');
set(gca, 'clim', [895 910]);
