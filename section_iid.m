%John Tran 25999001 FYP 2018
%plot of smoothing kernels f_Q(theta)
clc
clear
close all

Qa = 5; %Number of "receive" antennas
Qa_t = (Qa-1)/2; %Qa bar
Qb = 11; %Number of "receive" antennas
Qb_t = (Qb-1)/2; %Qb bar
t0_a = -0.23; %theta_0, corresponds to the physical location of a physical scatterer
t0_b = 0.13; %theta_0

%Values used to plot the graphs
theta = linspace(-0.5,0.5,1000);
qa = linspace(-Qa_t,Qa_t,Qa);
qb = linspace(-Qb_t,Qb_t,Qb);

xa = theta-t0_a;
xqa = (qa./Qa)-t0_a; %uniform sampling of xa
xb = theta-t0_b;
xqb = (qb./Qb)-t0_b; %uniform sampling of xb

%Graph a
J = (1/(Qa))*(exp(-1i*2*pi*xa*Qa_t));
f_Qa=J.*(sin(pi*Qa*xa)./sin(pi*xa));

%Corresponds to samples at the virtual angle
f_qa = (1/(Qa)).*(exp(-1i.*2.*pi.*xqa.*Qa_t)).*(sin(pi.*Qa.*xqa)./sin(pi.*xqa));

%Graph b
K = (1/(Qb))*(exp(-1i*2*pi*xb*Qb_t));
f_Qb = K.*sin(pi*Qb*xb)./sin(pi*xb);

%Corresponds to samples at the virtual angle
f_qb = (1/(Qb)).*(exp(-1i.*2.*pi.*xqb.*Qb_t)).*(sin(pi.*Qb.*xqb)./sin(pi.*xqb));

%% Plots
%a
figure()
hold on
plot(theta, abs(f_Qa))
ylabel('|f_Q(\theta - \theta_0)|')
xlabel('\theta')
title('plot a, \theta_0 = -0.23, Q = 5')

plot((qa./Qa),abs(f_qa),'*')
xlim([-0.5 0.5])
set(gca,'XTick',-0.5:0.1:0.5)
legend('|f_Q(\theta - \theta_0)|','|f_Q(q/Q - \theta_0)|')

%b
figure()
hold on
plot(theta, abs(f_Qb))
ylabel('|f_Q(\theta - \theta_0)|')
xlabel('\theta')
title('plot b, \theta_0 = 0.13, Q = 11')

plot((qb./Qb),abs(f_qb),'*')
xlim([-0.5 0.5])
set(gca,'XTick',-0.5:0.1:0.5)
legend('|f_Q(\theta - \theta_0)|','|f_Q(q/Q - \theta_0)|')