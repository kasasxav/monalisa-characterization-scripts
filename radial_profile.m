function [Tics,Average]=radial_profile(x, y, data,radial_step, xc, yc)
%main axii cpecified:
X= x -xc;
Y= y -yc;
% coordinate grid:
%[X,Y]=meshgrid(x,y);
% creating circular layers
Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % illustrating the principle:
% % figure;imagesc(Z_integer.*data)
% very fast MatLab calculations:
Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
Average=accumarray(Z_integer(:),data(:),[],@mean);
end