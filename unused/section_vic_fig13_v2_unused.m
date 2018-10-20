clc
clear
close all

P = 11;
Q = 11;

x = linspace(-0.4,0.4,P);
y = linspace(-0.4,0.4,Q);

L = 200;

a = 0.5;

%% a
figure()
pq1 = P*[-1/16,-3/16];
pq2 = P*[1/16,1/16];
pq3 = P*[3/16,-1/16];
pq4 = P*[-3/16,-1/16];

pq1 = round(pq1,0);
pq2 = round(pq2,0);
pq3 = round(pq3,0);
pq4 = round(pq4,0);

Hva = zeros(P,Q);
BL = sum(-1+(1+1)*rand(1,L));
BLp = mean(abs(BL)^2);
Hva(pq1(2)+5:pq1(2)+6,pq1(1)+6:pq1(1)+7) = BLp;
Hva(pq2(2)+5:pq2(2)+6,pq2(1)+6:pq2(1)+7) = BLp;
Hva(pq3(2)+5:pq3(2)+6,pq3(1)+6:pq3(1)+7) = BLp;
Hva(pq4(2)+5:pq4(2)+6,pq4(1)+6:pq4(1)+7) = BLp;
contour(x,y,Hva)

%% b
% figure()
% s2 = 1/8;
% pq1b = P*2*[-1/16,-3/16];
% pq2b = P*2*[1/16,1/16];
% pq3b = P*2*[3/16,-1/16];
% pq4b = P*2*[-3/16,-1/16];
% 
% pq1b = round(pq1b,0);
% pq2b = round(pq2b,0);
% pq3b = round(pq3b,0);
% pq4b = round(pq4b,0);
% 
% Hvb = zeros(P,Q);
% BL = sum(-1+(1+1)*rand(1,L));
% BLp = mean(abs(BL)^2);
% Hvb(pq1b(2)+5:pq1b(2)+6,pq1b(1)+6:pq1b(1)+7) = BLp;
% Hvb(pq2b(2)+5:pq2b(2)+6,pq2b(1)+6:pq2b(1)+7) = BLp;
% Hvb(pq3b(2)+5:pq3b(2)+6,pq3b(1)+6:pq3b(1)+7) = BLp;
% Hvb(pq4b(2)+5:pq4b(2)+6,pq4b(1)+6:pq4b(1)+7) = BLp;
% contour(x,y,Hvb)

figure()
s = 1/16;

Sr1 = Q*2*[-3/8,-3/8];
Sr1 = Q*[-3/8-s,-3/8+s];
Sr1 = round(Sr1,0);

St1 = P*2*[-1/8,-1/8];
St1 = P*[-1/8-s,-1/8+s];
St1 = round(St1,0);


Sr2 = Q*[1/8-s,1/8+s];
Sr2 = round(Sr2,0);
St2 = P*[1/8-s,1/8+s];
St2 = round(St2,0);

Sr3 = Q*[-1/16-s,-1/16+s];
Sr3 = round(Sr3,0);
St3 = P*[3/16-s,3/16+s];
St3 = round(St3,0);

Sr4 = Q*[1/16-s,-1/16+s];
Sr4 = round(Sr4,0);
St4 = P*[-3/16-s,-3/16+s];
St4 = round(St4,0);

Hvb = zeros(P,Q);
BL = sum(-1+(1+1)*rand(1,L));
BLp = mean(abs(BL)^2);
Hvb(Sr1(1)+6:Sr1(2)+6,St1(1)+6:St1(2)+6) = BLp;
Hvb(Sr2(1)+6:Sr2(2)+6,St2(1)+6:St2(2)+6) = BLp;
Hvb(Sr3(1)+6:Sr3(2)+6,St3(1)+6:St3(2)+6) = BLp;
Hvb(Sr4(1)+6:Sr4(2)+6,St4(1)+6:St4(2)+6) = BLp;
contour(x,y,Hvb)
