clear all
close all
%% Create a signal to test
t=linspace(-pi,pi,1000);
%x=[randn(length(t),1)]'+[sin(10*t)]; %Sin+ a bit of normal noise
%x=randn(length(t),1); %Normal noise
%x= randraw('bern', 0.5, [1 length(t)]); %Bernouli Distribution
x = randraw('norm', [], [1 length(t)]); %Another normal dist algorithm
%% As a test apply a moving average filter to the signal to see how stable test is
for windowSize = 1:100
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y = filter(b,a,x);
%Try four different measures
Krt(windowSize)=krt(y); %Regular Kurtosis
RobKurt(windowSize)=robkrt(y); %Robust approximation
Neg1(windowSize)=neg1(y); %Negentropy using hyperbolic approx
Neg2(windowSize)=neg2(y); %Negentropy using exponential approx
Neg3(windowSize)=neg3(y); %Negentropy using smoothed hyperbolic approx
end
%% As a test also apply a threshold filter
for windowSize = 101:200
y=x;
y2=y(y>=0);
y1=y(y<0);
y2((y2)<((200-windowSize)/100))=1.5;
y1(abs(y1)>((200-windowSize)/100))=-1.5;
y=[y1,y2];
%Try four different measures
Krt(windowSize)=krt(y); %Regular Kurtosis
RobKurt(windowSize)=robkrt(y); %Robust approximation
Neg1(windowSize)=neg1(y); %Negentropy using hyperbolic approx
Neg2(windowSize)=neg2(y); %Negentropy using exponential approx
Neg3(windowSize)=neg3(y); %Negentropy using smoothed hyperbolic approx
%% Plot the signal
figure(1), plot(x);
hold on
figure(1), plot(y,'*');
legend('Input Data','Filtered Data');
hold off
figure(2), hist(y);
end
%% Plots showing smoothness of Gaussianity measures
figure(3)
subplot(2,1,1), plot(1./abs(Krt));
legend('Kurtosis');
subplot(2,1,2), plot(1./abs(RobKurt));
legend('Robust Kurtosis');
figure();
plot(Neg1);
hold on
plot(Neg2);
plot((Neg1+Neg2)./2);
legend('Negentropy 1','Negentropy 2', 'Mix of both');

%% Functions
function y=krt(x)
   y=(mean(x.^4)-3*(mean(x.^2))^2);
end

function y=robkrt(x)    % more robust kurtosis
    y=(1/12)*((mean(x.^3))^2)+(1/48)*(krt(x))^2;
end

function y=neg1(x)
    g=randn(length(x),1);
    G=mean(log(cosh(g)));
    E=mean(log(cosh(x)));
    y=(E-G)^2;
end

function y=neg3(x)
    g=randn(length(x),1);
    G=0.5*mean(log(cosh(2*g)));
    E=0.5*mean(log(cosh(2*x)));
    y=(E-G)^2;
end

function y=neg2(x)
    g=randn(length(x),1);
    G=mean(-exp(-.5*(g.^2)));
    E=mean(-exp(-.5*(x.^2)));
    y=(E-G)^2;
end
