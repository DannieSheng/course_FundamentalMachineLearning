function [erle] = ERLE(d,e)
% This function implements SNR improvement in dB by the 
%            ERLE = 10*log(E{d^2}/E{e^2})
%
% INPUT
% d: desired signal
% e: error signal
%
% OUTPUT
% erle: SNR in dB
% 

D2 = mean(d.^2); %power of the desired signal
E2 = mean(e.^2); %power of the error signal

f = D2/E2; % ratio

erle = 10*log10(f); % dB
