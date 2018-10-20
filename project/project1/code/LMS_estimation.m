function [Wk,Ek,Y,R,Xk] = LMS_estimation(X,Dk,order,mu,method,random)
% function [Wk,Ek,Y,Xk] = LMS_estimation(X,Dk,order,eta,method,W)
%This funtion performs adaptive LMS/NLMS filter given the input data, desired
%output data, order and gain constant. It returns the weight values, the
%error, the output, the NMSE and auto-correlation function.
%
%   X     : input signal
%   Dk    : desired signal
%   order : order of LMS filter
%   mu    : learning rate (gain constant)
%   method: if 1, uses regular LMS. If 2, uses NLMS. default NLMS
%   random: if 1, uses random initialization of W. If 0, uses zeros as
%   weight initialization.
%
%   Wk    : weights
%   Ek    : error
%   Y     : output
% NMSE    : MSE normalized by the input power
%   R     : Auto-Correlation Function
%  MSE    : Mean Squared Error
%
%


% if nargin < 5
%     method = 2; %NLMS
%     random = 1; %random initialization of the weights
% end

if nargin < 6 
    W = zeros(order,1); %random initialization of the weights
end

Samples = length(X); %number of samples

% Initialization
Xk = zeros(order,length(X)-order); 
Ek=zeros(1,Samples);
Y=zeros(1,Samples);

% Input-delayed Matrix
for i=1:length(X)-order 
    Xk(:,i)=X(i+order-1:-1:i);
end

% Auto-correlation Matrix
% R = Xk*Xk'./length(Xk);

% % Choose between random or zero initialization of the weights
% if random==1
%     W = randn(order,1); %random
% else
%     W=zeros(order,1); %zeros
% end

if method ==1 % LMS algorithm
for k=1:Samples-order
    Y(k)= Xk(:,k)'*W; %output
    E = Dk(k) - Y(k); %instantaneous error
    Ek(:,k)=E; 
    W = W + 2 * eta * E * Xk(:,k); %weight update equation
    Wk(:,k)=W;
    NMSE(:,k) = mean((Dk-Xk'*W).^2)/mean(X.^2); %local NMSE
    MSE(:,k) = mean((Dk-Xk'*W).^2); %local MSE
end
elseif method ==2 % NLMS algorithm
    reg = 10^-10; %regularization term
    for k=1:Samples-order
    Y(k)= Xk(:,k)'*W; %output
    E = Dk(k) - Y(k); %instantaneous error
    Ek(:,k)=E;
    W = W + 2 * eta * E * Xk(:,k) ./ (reg+(Xk(:,k)'*Xk(:,k))); %weight update equation
    Wk(:,k)=W;
    NMSE(:,k) = mean((Dk-Xk'*W).^2)/mean(X.^2); %local NMSE
    MSE(:,k) = mean((Dk-Xk'*W).^2); %local MSE
    end
end


