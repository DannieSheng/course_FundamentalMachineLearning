%% EEL 5840/EEL4930 Elements of Machine Intelligence

clear
close all
clc

%% Initialize system
music = importdata('music.txt');
mix = importdata('corrupted_speech.txt');
fs = importdata('fs.txt');
music = (music-mean(music))';
mix = (mix-mean(mix))';
% w = 50000;
% w_start = 5;%140000;%
% music = music(w_start:w_start+w);
% mix = mix(w_start:w_start+w);
%%
tic
M_list = 5:5:100;
lambda_list = 0:0.01:0.5;
N = length(music);
MSE = zeros(length(M_list),length(lambda_list));
NMSE = zeros(length(M_list),length(lambda_list));
W = cell(length(M_list),length(lambda_list));
for idx_M = 1:length(M_list)
    d_mix = mix(M_list(idx_M):(end-1));
    
    for idx_lambda = 1:length(lambda_list)
        [W{idx_M,idx_lambda}, speech{idx_M,idx_lambda}, MSE(idx_M,idx_lambda), NMSE(idx_M,idx_lambda), erle_wiener(idx_M,idx_lambda)] =...
            Wiener_Estimization(music,d_mix,M_list(idx_M),lambda_list(idx_lambda));

    end
    display(['Filter order: ',num2str(M_list(idx_M)),'/',num2str(M_list(end))]);
end
toc

%% Plots
% [a_min,b_min]=find(MSE==min(min(MSE)));
% display(['MSE is minimum for filter order=',num2str(M_list(a_min)),...
%     ' and regularization=',num2str(lambda_list(b_min))]);
% figure,
% surf(M_list,lambda_list,MSE');xlabel('Filter order M');zlabel('Mean Square Error (MSE)');
% ylabel('Regularization parameter \lambda');colorbar;
% title(['Minimum MSE for M=',num2str(M_list(a_min)),' and \lambda=',num2str(lambda_list(b_min))]);

%%
% [c_min,d_min]=find(NMSE==min(min(NMSE)));
% display(['NMSE is minimum for filter order=',num2str(M_list(c_min)),...
%     ' and regularization=',num2str(lambda_list(d_min))]);

% figure,
% surf(M_list,lambda_list,NMSE');xlabel('Filter order M');zlabel('Normalized Mean Square Error (NMSE)');
% ylabel('Regularization parameter \lambda');colorbar;
% title(['Minimum NMSE for M=',num2str(M_list(c_min)),' and \lambda=',num2str(lambda_list(d_min))]);
% 
% figure,stem(W{a_min,b_min})
% xlabel('Time lags/taps');
% ylabel('Weight Coefficients');
% ylim([-1 1])
% title(['Weights for Filter Order = ' num2str(M_list(c_min))])
% 
% figure,stem(2:96, W{a_min,b_min}(2:96))
% xlabel('Time lags/taps');
% ylabel('Weight Coefficients');
% ylim([-0.3 0.15])
% title('Zoom from 2 to 96')

[M_max, lambda_max] = find(erle_wiener == max(max(erle_wiener)));
display(['ERLE is maximum for filter order=',num2str(M_list(M_max)),...
    ' and regularization=',num2str(lambda_list(lambda_max))]);
% figure
% plot(M_list, erle_wiener(:,lambda_max)', 'Linewidth',4)
% hold on
% title('ERLE Curve as a Function of Filter Order - Wiener Filter')
% legend(['\lambda = ' num2str(lambda_list(lambda_max))])

figure
plot(mix)
hold on
plot(speech{M_max, lambda_max})
legend('Corrupted Speech', 'Recovered Speech')
title('Comparison of Corrupted Speech and Recovered Speech - Wiener Filter')