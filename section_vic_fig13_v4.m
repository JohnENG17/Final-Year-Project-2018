%John Tran 25999001 FYP 2018
%E[|H_v(q,p)|^2] plots for flour clusters
clc
clear
close all

P = 11;
Q = 11;

x = linspace(-0.4,0.4,P);
y = linspace(-0.4,0.4,Q);

L = 200;

a= 1; %change a to get different plots

s = 1/16;

%initialising coordinates
Qc1 = [-3/16-s,-3/16+s]*Q*a*2;
Qc1 = round(Qc1);
Pc1 = [-1/16-s,-1/16+s]*P*a*2;
Pc1 = round(Pc1);

Qc2 = [1/16-s,1/16+s]*Q*a*2;
Qc2 = round(Qc2);
Pc2 = [1/16-s,1/16+s]*P*a*2;
Pc2 = round(Pc2);

Qc3 = [-1/16-s,-1/16+s]*Q*a*2;
Qc3 = round(Qc3);
Pc3 = [3/16-s,3/16+s]*P*a*2;
Pc3 = round(Pc3);

Qc4 = [-1/16-s,-1/16+s]*Q*a*2;
Qc4 = round(Qc4);
Pc4 = [-3/16-s,-3/16+s]*P*a*2;
Pc4 = round(Pc4);
%% Keeping the submatrices within the bounds of Hv
%cluster 1
for i= 1:2
if Qc1(i) <=-6
    Qc1(i) = -5;
end
if Pc1(i) <=-6
    Pc1(i) = -5;
end
end
if Qc1 == -5
    Qc1 = [0,5];
end

%cluster 2
for i= 1:2
if Qc2(i) >=6
    Qc2(i) = 5;
end
if Pc2(i) >=6
    Pc2(i) = 5;
end
end

%cluster 3
for i= 1:2
if Qc3(i) <=-6
    Qc3(i) = -5;
end
if Pc3(i) >=6
    Pc3(i) = 5;
end
end
if Pc3 == 5
    Pc3 = [-5,0];
end

%cluster 4
for i= 1:2
if Qc4(i) <=-6
    Qc4(i) = -5;
end
if Pc4(i) <=-6
    Pc4(i) = -5;
end
end
if Pc4 == -5
    Pc4 = [0,5];
end

%% Generating the contour plot
Hv = zeros(Q,P);

for i = Qc1(1):Qc1(2)
    for j = Pc1(1):Pc1(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hv(i+6,j+6) = BLp;  
    end
end
for i = Qc2(1):Qc2(2)
    for j = Pc2(1):Pc2(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hv(i+6,j+6) = BLp;  
    end
end
for i = Qc3(1):Qc3(2)
    for j = Pc3(1):Pc3(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hv(i+6,j+6) = BLp;  
    end
end
for i = Qc4(1):Qc4(2)
    for j = Pc4(1):Pc4(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hv(i+6,j+6) = BLp;  
    end
end

contour(x,y,Hv)
title(['a = ',num2str(a)])
xlabel('\theta_T')
ylabel('\theta_R')