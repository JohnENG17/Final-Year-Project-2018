%John Tran 25999001 FYP 2018
%plots of theta=alpha*sin(phi)
clc
clear
close all

%alpha is the normalized antenna spacing
a1 = 0.5; 
a2 = 1.9;
p = asin(0.5/a2);
p2 = asin((2-0.5)/a2);

%phi is the angle relative to the horizontal axis in fig 2 and 3
phi = linspace(-0.5*pi,0.5*pi);

phi_2 = linspace(-p,p);
phi_3 = linspace(p,p2);
phi_4 = linspace(-p2,-p);
phi_5 = linspace(p2,0.5*pi);
phi_6 = linspace(-0.5*pi,-p2);

theta_1 = a1*sin(phi);

theta_2 = a2*sin(phi_2); %-0.5<=theta<0.5 due to the periodicity of theta
theta_3 = a2*sin(phi_3);
theta_4 = a2*sin(phi_4);
theta_5 = a2*sin(phi_5);
theta_6 = a2*sin(phi_6);

%plot a
figure()
plot(phi, theta_1)
xlim([-0.5*pi 0.5*pi])
ylim([-0.5 0.5])
xlabel("\phi/\pi")
ylabel("\theta")
title("\alpha = 0.5")
set(gca,'XTick',-0.5*pi:pi/4:0.5*pi)
set(gca,'XTickLabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})

%plot b
figure()
plot(phi_2, theta_2, '.b')
xlim([-0.5*pi 0.5*pi])
ylim([-0.5 0.5])
xlabel("\phi/\pi")
ylabel("\theta")
title("\alpha = 1.9")
set(gca,'XTick',-0.5*pi:pi/4:0.5*pi)
set(gca,'XTickLabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})

hold on
%plot(phi_3, theta_3, '.r')
plot(phi_3, theta_3-1, '.b')
plot(phi_4, theta_4+1, '.b')
plot(phi_5, theta_5-2, '.b')
plot(phi_6, theta_6+2, '.b')

% %% Some working out 

theta = a2*sin(phi);
plot(phi,theta)
xlabel("\phi/\pi")
ylabel("\theta")
