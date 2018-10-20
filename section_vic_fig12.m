%John Tran 25999001 FYP 2018
%E[|H_v(q,p)|^2] plots for single cluster
clc
clear
close all

%Initialising all the variables
P = 11;
Q = 11;

x = linspace(-0.4,0.4,P);
y = linspace(-0.4,0.4,Q);

L = 200;

a = 0.5;

%% a
%we know the cluster is pi/4 X pi/4 wide, we use this to find the angular spreads
%phi (0,0)
Sr = zeros(1,2);
St = zeros(1,2);

Sr(1) = -pi/8;
Sr(2) = pi/8;

St(1) = -pi/8;
St(2) = pi/8;

%determining size of cluster
Qc = zeros(1,2);
Pc = zeros(1,2);

Qc(1) = a*Q*sin(Sr(1));
Qc(2) = a*Q*sin(Sr(2));
Qc = round(Qc,0);
Pc(1) = a*P*sin(St(1));
Pc(2) = a*P*sin(St(2));
Pc = round(Pc,0);

Hva = zeros(P,Q);
%Calculating the channel power E(|Hv(q,p)|^2)
for i = Qc(1):Qc(2)
    for j = Pc(1):Pc(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hva(i+6,j+6) = BLp;  
    end
end
contour(x,y,Hva)
title('a=0.5')
xlabel('\theta_T')
ylabel('\theta_R')

sigma = trace(Hva'*Hva);

%% b
% Identicle to part a, but instead alpha = 1
figure()
a_b = 1;

%determining size of cluster
Qc_b = zeros(1,2);
Pc_b = zeros(1,2);

Qc_b(1) = a_b*Q*sin(Sr(1));
Qc_b(2) = a_b*Q*sin(Sr(2));
Qc_b = round(Qc_b,0);
Pc_b(1) = a_b*P*sin(St(1));
Pc_b(2) = a_b*P*sin(St(2));
Pc_b = round(Pc_b,0);

Hvb = zeros(P,Q);
%Calculating the channel power E(|Hv(q,p)|^2)
for i = Qc_b(1):Qc_b(2)
    for j = Pc_b(1):Pc_b(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hvb(i+6,j+6) = BLp;  
    end
end
contour(x,y,Hvb)
title('a=1')
xlabel('\theta_T')
ylabel('\theta_R')

%% c
% Identicle to part a, but instead alpha = 1.31
figure()
a_c = 1.31;

%determining size of cluster
Qc_c = zeros(1,2);
Pc_c = zeros(1,2);

Qc_c(1) = a_c*Q*sin(Sr(1));
Qc_c(2) = a_c*Q*sin(Sr(2));
Qc_c = round(Qc_c,0);

Pc_c(1) = a_c*P*sin(St(1));
Pc_c(2) = a_c*P*sin(St(2));
Pc_c = round(Pc_c,0);

%The for loops keeps the submatrix within the bounds of Hv by setting the
%bounds within [-5,5]
for i= 1:2
if Qc_c(i) <=-6
    Qc_c(i) = -5;
end
if Pc_c(i) <=-6
    Pc_c(i) = -5;
end
end

for i= 1:2
if Qc_c(i) >=6
    Qc_c(i) = 5;
end
if Pc_c(i) >=6
    Pc_c(i) = 5;
end
end

Hvc = zeros(P,Q);
for i = Qc_c(1):Qc_c(2)
    for j = Pc_c(1):Pc_c(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hvc(i+6,j+6) = BLp;  
    end
end
contour(x,y,Hvc)
title('a=1.31')
xlabel('\theta_T')
ylabel('\theta_R')










