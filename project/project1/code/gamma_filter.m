clear
close all;
clc;
dbstop if error
%% Initialize system
fs = importdata('fs.txt');
music = importdata('music.txt');
mix = importdata('corrupted_speech.txt');

music = (music-mean(music))';
mix = (mix-mean(mix))';

tic
M_list = 5:5:100;
eta_list = [10^(-5) 10^(-4) 5*10^(-4) 10^(-3) 5*10^(-3)];
N = length(mix);

mu_list = 0.2;
idx_mu = 1;

tic
display('Order | Step-size');
% for iter = 1:20
iter = 1;
for idx_M = 1:length(M_list)
    d_mix = mix(M_list(idx_M):end);
    for idx_eta = 1:length(eta_list)
        ww{idx_M, idx_eta}(:,1) = zeros(M_list(idx_M),1);
        [Wk{idx_M, idx_eta}, Ek{idx_M, idx_eta}, ~, Xk{idx_M}] = GAMMA_estimization(music,d_mix,M_list(idx_M),eta_list(idx_eta),1,mu_list(idx_mu), ww{idx_M, idx_eta}(:,iter));
%         [Wk{idx_M, idx_ita}, Ek{idx_M, idx_ita}, ~, Xk{idx_M}, MSE{idx_M, idx_ita}] = GAMMA_estimization(music,d_mix,M_list(idx_M),ita_list(idx_ita),1,mu_list(idx_mu), ww{idx_M, idx_ita}(:,iter));

        speech{idx_M, idx_eta} = d_mix - Xk{idx_M}'* Wk{idx_M, idx_eta}(:,end);
        [erle{idx_M, idx_eta}] = ERLE(d_mix,speech{idx_M, idx_eta});
        ww{idx_M, idx_eta}(:,iter+1) = Wk{idx_M, idx_eta}(:,end);
    end
    display(['Order ',num2str(M_list(idx_M)),'/',num2str(length(M_list)),' done!'])
end
% end
% save('Gamma_weight.mat', 'ww')
[erle{idx_M, idx_eta}] = ERLE(d_mix,speech{idx_M, idx_eta});
toc
erle = cell2mat(erle);
[idx_Mmax, idx_itamax] = find(erle == max(max(erle)));
display(['ERLE is maximum for filter order=',num2str(M_list(idx_Mmax)),' and step size=',num2str(eta_list(idx_itamax))]);

figure
plot(M_list,erle,'Linewidth',2)
legend(['\eta = ' num2str(eta_list(1))])
xlabel('Filter Order')
ylabel('Echo Return Loss Enhancement')
title(['ERLE curve as a function of the filter order - Gamma Filter with \mu = ' num2str(mu_list(idx_mu))])

figure
plot(mix)
hold on
plot(speech{idx_Mmax})
title(['Comparison of Corrupted Speech and Recovered Speech - Gamma Filter with \mu = ' num2str(mu_list(idx_mu))])
legend('Corrupted Speech', 'Recovered Speech')

