function [ projected_data, w ] = PCA( X, d)
%X: data
%d: the desired dimension of the output data
%   Detailed explanation goes here
mu = mean(X);
X_std = X - mu;
cov_mat = cov(X_std);
[eigenVecs, eigenVals] = eig(cov_mat);

% plot PCA
% eV = diag(eigenVals);
% s = 0;
% for i = size(eigenVals,1):-1:1
%     s = s+eV(i);
%     ratio(i) = s/sum(eV);
% end
% figure, plot(ratio(end:-1:1), 'Linewidth', 3)
% xx = 1:4:784;
% yy = 0.85*ones(size(xx));
% hold on, plot(xx,yy, '.', 'Linewidth', 2)
% xlabel('number of eigenvalues')
% ylabel('fraction of total variance retained')
% title('Fraction of total variance retained vesus number of eigenvalues')

w = eigenVecs(:,(end-d+1):end);
projected_data = X_std * w;

end

