clc
clear
close all

P = 11;
Q = 11;

x = linspace(-0.4,0.4,P);
y = linspace(-0.4,0.4,Q);

L = 200;

%% a
%s = is the size of the clusters divided by 2
s = pi/16;

% Sr1 = [-3*pi/8-s,-3*pi/8+s];
% St1 = [-pi/8-s,-pi/8+s];
% 
% Sr2 = [pi/8-s,pi/8+s];
% St2 = [pi/8-s,pi/8+s];
% 
% Sr3 = [-pi/8-s,-pi/8+s];
% St3 = [3*pi/8-s,3*pi/8+s];
% 
% Sr4 = [-pi/8-s,-pi/8+s];
% St4 = [-3*pi/8-s,-3*pi/8+s];
%==========================================================================
Sr1 = [-3*pi/16-s,-3*pi/16+s];
St1 = [-pi/16-s,-pi/16+s];

Sr2 = [pi/16-s,pi/16+s];
St2 = [pi/16-s,pi/16+s];

Sr3 = [-pi/16-s,-pi/16+s];
St3 = [3*pi/16-s,3*pi/16+s];

Sr4 = [-pi/16-s,-pi/16+s];
St4 = [-3*pi/16-s,-3*pi/16+s];
%==========================================================================
% a = 0.25;
% Qc1(1) = a*Q*sin(Sr1(1))-1;
% Qc1(2) = a*Q*sin(Sr1(2))-1;
% Qc1 = round(Qc1,0);
% Pc1(1) = a*P*sin(St1(1))-1;
% Pc1(2) = a*P*sin(St1(2));
% Pc1 = round(Pc1,0);
% 
% Qc2(1) = a*Q*sin(Sr2(1));
% Qc2(2) = a*Q*sin(Sr2(2));
% Qc2 = round(Qc2,0);
% Pc2(1) = a*P*sin(St2(1));
% Pc2(2) = a*P*sin(St2(2))+1;
% Pc2 = round(Pc2,0);
% 
% Qc3(1) = a*Q*sin(Sr3(1));
% Qc3(2) = a*Q*sin(Sr3(2));
% Qc3 = round(Qc3,0);
% Pc3(1) = a*P*sin(St3(1))+1;
% Pc3(2) = a*P*sin(St3(2))+1;
% Pc3 = round(Pc3,0);
% 
% Qc4(1) = a*Q*sin(Sr4(1))-1;
% Qc4(2) = a*Q*sin(Sr4(2));
% Qc4 = round(Qc4,0);
% Pc4(1) = a*P*sin(St4(1))-1;
% Pc4(2) = a*P*sin(St4(2))-1;
% Pc4 = round(Pc4,0);
%==========================================================================
% %% b
a= 0.5;
Qc1(1) = a*Q*sin(Sr1(1))-1;
Qc1(2) = a*Q*sin(Sr1(2))-1;
Qc1 = round(Qc1,0);
Pc1(1) = a*P*sin(St1(1))-1;
Pc1(2) = a*P*sin(St1(2));
Pc1 = round(Pc1,0);

Qc2(1) = a*Q*sin(Sr2(1));
Qc2(2) = a*Q*sin(Sr2(2))+1;
Qc2 = round(Qc2,0);
Pc2(1) = a*P*sin(St2(1));
Pc2(2) = a*P*sin(St2(2))+1;
Pc2 = round(Pc2,0);

Qc3(1) = a*Q*sin(Sr3(1))-1;
Qc3(2) = a*Q*sin(Sr3(2));
Qc3 = round(Qc3,0);
Pc3(1) = a*P*sin(St3(1))+1;
Pc3(2) = a*P*sin(St3(2))+1;
Pc3 = round(Pc3,0);

Qc4(1) = a*Q*sin(Sr4(1))-1;
Qc4(2) = a*Q*sin(Sr4(2));
Qc4 = round(Qc4,0);
Pc4(1) = a*P*sin(St4(1))-1;
Pc4(2) = a*P*sin(St4(2))-1;
Pc4 = round(Pc4,0);
%==========================================================================
%% c
% a= 1;
% Qc1(1) = a*Q*sin(Sr1(1))+3;
% Qc1(2) = a*Q*sin(Sr1(2))+3;
% Qc1 = round(Qc1,0);
% Pc1(1) = a*P*sin(St1(1))-1;
% Pc1(2) = a*P*sin(St1(2))+1;
% Pc1 = round(Pc1,0);
% 
% Qc2(1) = a*Q*sin(Sr2(1))+1;
% Qc2(2) = a*Q*sin(Sr2(2))+1;
% Qc2 = round(Qc2,0);
% Pc2(1) = a*P*sin(St2(1))-1;
% Pc2(2) = a*P*sin(St2(2))+1;
% Pc2 = round(Pc2,0);
% 
% Qc3(1) = a*Q*sin(Sr3(1))-1;
% Qc3(2) = a*Q*sin(Sr3(2))+1;
% Qc3 = round(Qc3,0);
% Pc3(1) = a*P*sin(St3(1))-3;
% Pc3(2) = a*P*sin(St3(2))-3;
% Pc3 = round(Pc3,0);
% 
% Qc4(1) = a*Q*sin(Sr4(1))+4;
% Qc4(2) = a*Q*sin(Sr4(2))+5;
% Qc4 = round(Qc4,0);
% Pc4(1) = a*P*sin(St4(1))+3;
% Pc4(2) = a*P*sin(St4(2))+3;
% Pc4 = round(Pc4,0);
%==========================================================================
for i= 1:2
if Qc1(i) <=-6
    Qc1(i) = -5;
end
if Pc1(i) <=-6
    Pc1(i) = -5;
end
end

for i= 1:2
if Qc2(i) >=6
    Qc2(i) = 5;
end
if Pc2(i) >=6
    Pc2(i) = 5;
end
end

for i= 1:2
if Qc3(i) <=-6
    Qc3(i) = -5;
end
if Pc3(i) >=6
    Pc3(i) = 5;
end
end

for i= 1:2
if Qc4(i) <=-6
    Qc4(i) = -5;
end
if Pc4(i) <=-6
    Pc4(i) = -5;
end
end

Hva = zeros(P,Q);
% BL = sum(-1+(1+1)*rand(1,L));
% BLp = mean(abs(BL)^2);
% Hva(Qc1(1)+6:Qc1(2)+6,Pc1(1)+6:Pc1(2)+6) = BLp; %1
% Hva(Qc4(1)+6:Qc4(2)+6,Pc4(1)+6:Pc4(2)+6) = BLp; %2
% Hva(Qc2(1)+6:Qc2(2)+6,Pc2(1)+6:Pc2(2)+6) = BLp; %3
% Hva(Qc3(1)+6:Qc3(2)+6,Pc3(1)+6:Pc3(2)+6) = BLp; %4

for i = Qc1(1):Qc1(2)
    for j = Pc1(1):Pc1(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hva(i+6,j+6) = BLp;  
    end
end

for i = Qc2(1):Qc2(2)
    for j = Pc2(1):Pc2(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hva(i+6,j+6) = BLp;  
    end
end

for i = Qc3(1):Qc3(2)
    for j = Pc3(1):Pc3(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hva(i+6,j+6) = BLp;  
    end
end

for i = Qc4(1):Qc4(2)
    for j = Pc4(1):Pc4(2)
      BL = -1+(1+1)*rand(1,L);
      BLp = mean(abs(BL).^2);
      Hva(i+6,j+6) = BLp;  
    end
end

contour(x,y,Hva)




