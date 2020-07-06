
clear all;close all

load dSE;       % Death data from FHM 20200311 to 20200703
D=100;          % Days considered in simulation
d=dSE(1:D);     % Remove last days which have missing data
D1=20;          % Intervention day

% Parameters that fit data
N=10e6;             % Population
I0=7.7226e+04;      % Initially infected
gamma=0.0195;       % SIR parameter gamma
lambda0 = 0.2847;   % SIR parameter lambda before intervention
lambda1 = 0.2549;    % SIR parameter lambda after intervention
cfr=6.5867e-04;     % Case fatality ratio

% Initial states
S(1)=N; 
I(1)=I0;
R(1)=0;  

for t=1:D-1
    if t<=D1
        lambda=lambda0;
    else
        lambda=lambda1;
    end
    S(t+1)=S(t)-lambda*I(t)*S(t)/N;
    I(t+1)=I(t)+lambda*I(t)*S(t)/N - gamma*I(t);
    R(t+1)=R(t)+gamma*I(t);
end
dhat=[0; cfr*diff(R(1:D)')];

figure
plot(1:D,d,'b-',1:D,dhat,'r-',D1,dhat(D1),'or')
legend({'Actual deaths','SIR model','Break point'})
xlabel('Days after March 11')
set(gca,'xlim',[0 length(dhat)],'ylim',[0 100])

figure
semilogy(1:D,S,'b-',1:D,I,'r-',1:D,R,'g-')
legend({'S','I','R'})
