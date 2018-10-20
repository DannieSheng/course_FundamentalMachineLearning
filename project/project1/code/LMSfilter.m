close all
clc
dbstop if error
load('music.txt')
load('corrupted_speech.txt')
load('fs.txt')

music = (music - mean(music))';
corrupted_speech = (corrupted_speech - mean(corrupted_speech))';

%%
M_list = 5:5:100;
ita_list = [10^(-6) 10^(-5) 10^(-4) 10^(-3)];
% ita_list = []
N = length(music);
MSE = zeros(length(M_list), length(ita_list));
% Wk = cell(length(M_list) , length(ita_list));
tic

for iter = 1:100;%:max_iter
    for idx_M = 1:length(M_list)  
        d_mix = corrupted_speech(M_list(idx_M):(end-1));
        for idx_ita = 1:length(ita_list)
        ww{idx_M, idx_ita}(:,1) = zeros(M_list(idx_M),1);
%         [Wk{idx_M, idx_ita}, Ek{idx_M, idx_ita}, ~, ~,Xk{idx_M}] = LMS_prediction(music, d_mix, M_list(idx_M), ita_list(idx_ita), 1, 0);
            [Wk{idx_M, idx_ita},Ek{idx_M, idx_ita},~,Xk{idx_M}] = LMS_estimation(music,d_mix,M_list(idx_M),ita_list(idx_ita),1,ww{idx_M, idx_ita}(:,iter));
            speech{idx_M, idx_ita} = d_mix - Xk{idx_M}'* Wk{idx_M, idx_ita}(:,end);
            [erle{idx_M, idx_ita}] = ERLE(d_mix,speech{idx_M, idx_ita});
            ww{idx_M, idx_ita}(:,iter+1) = Wk{idx_M, idx_ita}(:,end);
        end 
        
    end
    display(['Order ',num2str( M_list(idx_M)),'/',num2str(ita_list(idx_ita)),' done!'])
end