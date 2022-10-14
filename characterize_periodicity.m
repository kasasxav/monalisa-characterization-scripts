%% Characterize periodicity over FOV

%% 1. Clear workspace
clear;
clc;

%% 2. Load image and beads coordinates
path = '\\testalab4.ad.scilifelab.se\D\Data\2020_full\2020-09-17\analysis\488\';
imm = h5read([path '-output_interp-488.hdf5'],'/data');

load([path 'xyBeads.mat']);

%% 3. Line profiles
scan = 1 / 0.05;
N = size(xyBeads, 1);
p = zeros(N, 2);
r = zeros(N, 2);
n = 6;
Cx = -0.5558;
Cy = -0.5206;
for b = 1:N
    x = xyBeads(b, 1);
    y = xyBeads(b, 2);
    if (x-scan > 0) && (y-scan >0) && (x < size(imm,1)) && (y < size(imm,2))
        rec_square = imm(x-scan:x, y-scan:y, n);

        % Horizontal direction
        line_x = sum(rec_square, 2);
        line_x = (line_x - min(line_x(:)))/(max(line_x(:))-min(line_x(:))) - 1/2;
        l = (1:size(line_x,1));
        [fit, gof] = sinFit(l, line_x);
        coeff = coeffvalues(fit);
        p(b, 1) = coeff(2);
        r(b, 1) = gof.rsquare;

        % Horizontal direction
        line_y = sum(rec_square, 1);
        line_y = (line_y - min(line_y(:)))/(max(line_y(:))-min(line_y(:))) - 1/2;
        l = (1:size(line_y,2));
        [fit, gof] = sinFit(l, line_y);
        coeff = coeffvalues(fit);
        p(b, 2) = coeff(2);
        r(b, 2) = gof.rsquare;    
    end
end
%% 4. Plots
p_t = (p(:,1) + p(:,2))/2;
gamma = 0.45;
cond = (r(:,1)>gamma & r(:,2)>gamma);
figure;
%subplot(2,2,1); 
scatter((xyBeads(cond,1))/(10)-27.9,xyBeads(cond,2)/(10)-60.6, 1000, 25*2*pi./p_t(cond,1), '.'); 
%scatter((xyBeads(cond,1))/(10)-27.9,xyBeads(cond,2)/(10)-60.6, 1000, 25*2*pi./p_t(cond,1), '.'); 
set(gca,'FontSize',24)
%title('Periodicity - V', 'FontSize', 24);
%xlabel('X [\mum]', 'FontSize', 24);
xlim([0 130])
%ylabel('Y [\mum]', 'FontSize', 24);
ylim([0 130])
% colorbar;
% caxis([230, 245]);
% subplot(2,2,2); scatter(xyBeads(r(:,1)>gamma,1)/(10),xyBeads(r(:,1)>gamma,2)/(10), 1000,(r(r(:,1)>gamma,1)+r(r(:,2)>gamma,2)/2) , '.'); 
% set(gca,'FontSize',24)
% title('R value - V', 'FontSize', 24);
% xlabel('X [\mum]', 'FontSize', 24);
% xlim([0 200])
% ylabel('Y [\mum]', 'FontSize', 24);
% ylim([0 210])
% colorbar;
% subplot(2,2,3); scatter(xyBeads(r(:,2)>gamma,1)/(10),xyBeads(r(:,2)>gamma,2)/(10), 1000, 25*2*pi./p(r(:,2)>gamma,2), '.'); 
% set(gca,'FontSize',24)
% title('Periodicity - H', 'FontSize', 24);
% xlabel('X [\mum]', 'FontSize', 24);
% xlim([0 200])
% ylabel('Y [\mum]', 'FontSize', 24);
% ylim([0 210])
% colorbar;
% caxis([230, 245]);
% subplot(2,2,4); scatter(xyBeads(r(:,2)>gamma,1)/(10),xyBeads(r(:,2)>gamma,2)/(10), 1000, r(r(:,2)>gamma,2), '.'); 
% set(gca,'FontSize',24)
% title('R value - H', 'FontSize', 24);
% xlabel('X [\mum]', 'FontSize', 24);
% xlim([0 200])
% ylabel('Y [\mum]', 'FontSize', 24);
% ylim([0 210])
% colorbar;

set(gca, 'clim', [210 240]);
