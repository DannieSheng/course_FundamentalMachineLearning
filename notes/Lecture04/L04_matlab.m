
% EEL5840/EEL4930: Elements of Machine Intelligence Lecture 4

clear all
close all
clc
%% Polymonial Regression - Training

N = 1000; % number of samples
M = 4; % order of the polynomial
noisep = 0.2; % variance of normally distributed noise

% First, lets generate some simulated x- data.
input = linspace(0,2*pi,N); % input without noise

samples = input + noisep.*randn(1,N); % input data corroborated with normally distributed noise 
output = sin(samples)'; % output data

t = sin(input)'; % desired vector

% Data matrix as in the notes
X = power(repmat(input',1,M+1),repmat(0:M,N,1));

R = transpose(X)*X; % or you can write X'*X
p = X'*t; % or you can write transpose(X)*t

w = inv(R)*p;

display(w)
figure, stem(w,'LineWidth',2); axis([0.5,M+0.5,min(w)-1,max(w)+1]);
xlabel('Model order M','FontSize',13);ylabel('Coefficients w','FontSize',13);
title('Vector w','FontSize',15);

% Polynomial Regression - Test

Ntest = 100;

xranget = linspace(0,2*pi,Ntest);

% test data matrix X (format as in the notes)
Xtest = power(repmat(xranget',1,M+1),repmat(0:M,Ntest,1));

ttest = sin(xranget)'; % desired vector 
esty = Xtest*w;       % estimatied polynomial 

figure,
plot(input,output,'ob');hold on;
plot(xranget,esty,'-r','LineWidth',2);hold on;
plot(xranget,ttest,'-g','LineWidth',2);hold off;
xlabel('Input x','FontSize',13);ylabel('Desired t','FontSize',13);
legend('Training Data','Estimated Polynomial','True function');
axis([-0.5,M+2,-1,1]);
title([num2str(M),'th-order Polynomial Regression'],'FontSize',15);


%% Analysis of the Error

inst_error = esty-ttest; %intantaneous error
sq_error = inst_error.^2; %squared error

% plot everything
figure, 
subplot(2,2,1); hist(inst_error,20);xlabel('Input x','FontSize',13);
ylabel('Instantaneous Error e = y-t','FontSize',13);
title('Histogram of the Instantaneous Error','FontSize',15);

subplot(2,2,3); ksdensity(inst_error);xlabel('Input x','FontSize',13);
ylabel('Instantaneous Error e = y-t','FontSize',13);
title('Estimated Density of the Instantaneous Error','FontSize',15);

subplot(2,2,2); hist(sq_error,20);xlabel('Input x','FontSize',13);
ylabel('Squared Error e^Te = (y-t)^T(y-t)','FontSize',13);
title('Histogram of the Squared Error','FontSize',15);

subplot(2,2,4); ksdensity(sq_error);xlabel('Input x','FontSize',13);
ylabel('Squared Error e^Te = (y-t)^T(y-t)','FontSize',13);
title('Estimated Density of the Squared Error','FontSize',15);

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
    'Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Analysis of the Instantaneous and Squared Error',...
    'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',15)

%% Auto- and Cross-Correlation Matrices

R = transpose(X)*X; % or you can write X'*X
p = X'*t; % or you can write transpose(X)*t

% eigendecomposition
[V,D] = eig(R);

% V contains the eigenvectors in each column
% D contains the associated eigenvalues in its diagonal
eigenvalues = diag(D); % eigenvalues are now displayed in a vector

display('Eigenvalues of the auto-correlation matrix:');
display(sort(eigenvalues,'descend'));

RtimesV1 = R*V(:,1); 
Lambda1V1 = eigenvalues(1)*V(:,1);

display(RtimesV1)
display(Lambda1V1)

%% Polynomial Curve Fitting and Ill-Conditioned Systems

M = 7;
N = 100;
noisep = .30;

% First, let's generate some simulated x- data.  
input = linspace(0,1,N);

% Next, let's generate some noise
e = noisep.*randn(1,N);

% Suppose the true function is a sine curve and add the noise
t = (sin(2*pi.*input) + e)';

% Then we can fit the data using the polynomial curve fitting method we derived
X = power(repmat(input',1,M+1),repmat(0:M,N,1));

R = X'*X; % auto-correlation matrix
Rnoisy = X'*X + 0.1*eye(M+1,M+1); %diagonally-loaded auto-correlation

w = inv(R)*transpose(X)*t; % coefficients w
wnoisy = inv(Rnoisy)*transpose(X)*t; %coefficients w with perturbed R matrix

display(w)
display(wnoisy)

figure,stem(w,'LineWidth',2);hold on;stem(wnoisy,'Color','red','LineWidth',2); 
axis([0.5,M+.5,min(w)-1,max(w)+1]);
legend('Weights','Weights with Perturbation in R=X^TX');
xlabel('Model order M','FontSize',13);ylabel('Coefficients w','FontSize',13);
title('Vector w','FontSize',15);

% Now let us use the weights in test and plot results
xrange = linspace(0,1,N);  %get equally spaced points in the xrange
y = sin(2*pi.*xrange)'; %compute the true function value
X = [ones(N,1), power(repmat(xrange',1,M),repmat(1:M,N,1))];
esty = X*w; %compute the predicted value
esty_w = X*wnoisy;

% plot everything
figure,
plot(xrange,y,'-g','LineWidth',2); hold on;
plot(input,t,'ob');hold on;
plot(xrange,esty,'-r','LineWidth',2);hold on;
plot(xrange,esty_w,'c','LineWidth',2);hold off;
legend('True Function','Training Data','Estimated Polynomial','Loaded Polynomial');
xlabel('Input x');ylabel('Desired t');

[V,D] = eig(R);
[Vnoisy,Dnoisy] = eig(Rnoisy);

display(['Condition number of autocorrelation matrix: ',...
    num2str(max(diag(D))/min(diag(D)))]); 
display(['Condition number of diagonally-loaded autocorrelation matrix: ', ...
    num2str(max(diag(Dnoisy))/min(diag(Dnoisy)))]);
display(' ')
display('Eigenspectrum of autocorrelation matrix: ')
display(sort(diag(D),'descend'))
display(' ')
display('Eigenspectrum of diagonally-loaded autocorrelation matrix: ');
display(sort(diag(Dnoisy),'descend'))


%% Gradient Vector and Hessian Matrices

syms x y z

f = x*y + 2*z*x;
gradient(f, [x,y,z])
hessian(f, [x,y,z])

g = x^2 + y^2 + z^2;
gradient(g, [x,y,z])
hessian(g,[x,y,z])

h = x^3 + y^3 - 3*x*y;
gradient(h, [x,y])
hessian(h,[x,y])







