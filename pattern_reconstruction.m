%% Reconstruct an image from a 3D stack of 2D stacks, only taking the beads in focus

%% 1. Clear and close
clear all
close all
%% 2. Selecting data folder
datadir = '\\testalab4\D\OneDrive - KTH\2021_full\2022-03-02\2d\';
listing = dir(datadir);
I = size(listing, 1);

%% 3. Load image
i = 4;
fn = split(listing(i).name, '.')
image = double(loadimage(datadir, fn, 1));

%% 4. Localize
im_1 = imgaussfilt(max(image,[],3),2);
threshold=10;
pxSize=0.1;
xyBeads = Localize(im_1, threshold, 4, 5,1, pxSize);

h1 = figure(1);
imagesc(im_1,[min(min(im_1)) max(max(im_1))/0.5]), axis image, colormap jet, colorbar
title(['Localized beads']);
hold on
plot(xyBeads(:, 1), xyBeads(:, 2),'xr')
%% 5. Reconstruct fitting image
Cx = -0.5;
Cy = -0.5;
nFrames = 324;
m = sqrt(nFrames);
xr = zeros(size(xyBeads,1), nFrames);
yr = zeros(size(xyBeads,1), nFrames);
value = zeros(size(xyBeads,1), nFrames);
[Ny, Nx] = size(im_1);
for n = 1:nFrames
    % For each frame localize bead and construct image
    image = double(loadimage(datadir, fn, n));
        for b = 1:size(xyBeads,1)
            x = xyBeads(b,1) + Cx*floor((n-1)/m);
            y = xyBeads(b,2) + Cy*mod(n-1, m);

            yr(b, n) = max(1,min(2*Ny,2*y));
            xr(b, n) = max(1,min(2*Nx,2*x));

            r = 10;
            crop = image(max(1,floor(y-r)):min(Ny, ceil(y+r)), max(1,floor(x-r)):min(Nx, ceil(x+r)));
            value(b, n) = mean(mean(crop));
         end
    end


%% 11. Create full image
x_b = 2*Cx*(0:m-1);
y_b = 2*Cy*(0:m-1);

rec_im = zeros(2*Nx, 2*Ny);
for b = 1:size(value,1)
    xx = xr(b,1) + x_b;
    yy = yr(b,1) + y_b;
    
    [xx,  yy] = meshgrid(xx, yy);
    [xxx, yyy] = meshgrid(1:2*Nx, 1:2*Ny);
    
    square = reshape(value(b,:), [m, m]);
    rec_im = rec_im + interp2(xx, yy, square, xxx, yyy, 'linear', 0);
end

%% 12. Save image
savename = strcat(datadir, 'rec.hdf5');
h5create(savename,'/data',[2*Ny, 2*Nx]);
h5write(savename, '/data', rec_im);