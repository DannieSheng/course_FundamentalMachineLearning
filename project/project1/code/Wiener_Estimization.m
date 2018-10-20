function [W_optimal,E_wiener,MSE,NMSE,erle_wiener] = Wiener_Estimization(input,desire,order,reg)
% This function implements the Wiener solution for Echo Cancellation
%
% INPUT
% input: input signal
% desire: desired signal
% order: filter order
% reg: (optional) regularizer parameter. Default value is reg=0.
%
% OUTPUT
% W_optimal: analytical weights
% E_wiener: error signal
% MSE: mean squared error
% NMSE: normalized mean squared error 
% erle_wiener: ERLE 
%
%

if nargin < 4
    reg = 0;
end

Xmean=mean(input.^2);       %input power

% D = desire(order+1:end);      %desired response
D = desire;

% Construction of the input matrix
DataMatrix = zeros(order,length(input)-order);
for i = 1:(length(input)-order)
    DataMatrix(:,i)=input(i+order-1:-1:i);
end

% Computation of R and P
R = DataMatrix*DataMatrix'/length(DataMatrix);     %auto-correlation
P = DataMatrix*D/length(DataMatrix);      %cross-correlation

% Adding regularizer term to R
Rreg = R + reg.*eye(order);

% Optimal Weights
W_optimal = Rreg\P;       %weights

% Error
E_wiener = D - DataMatrix'*W_optimal;

% MSE
MSE = mean(E_wiener.^2);

% Normalized MSE
NMSE = mean(E_wiener.^2)/Xmean;

% ERLE
erle_wiener = ERLE(D,E_wiener);


