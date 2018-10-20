clc
clear
close all

P = 11;
Q = 11;

x = linspace(-0.4,0.4,P);
y = linspace(-0.4,0.4,Q);
q = linspace(-(Q-1)/2,(Q-1)/2,Q)/Q;
p = linspace(-(P-1)/2,(P-1)/2,P)/P;

rho = 100;
SNR = 10*log10(rho);
k=10;
sf = sqrt( P^2/(P+k*(2*P-k-1)) );

L = 200;

a = 0.5;



%% constructing the cluster
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

%% plotting

Hva = zeros(P,Q); %for the contour plot
Hvac = zeros(P,Q); %for the capacity plot

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

%array steering and response vectors
for i = 1:P
    for j = 1:P
ar(i,j) = 1/sqrt(P)*( exp(-j*2*pi*p(j)*(P-12+i)) )';
at(i,j) = 1/sqrt(Q)*( exp(-j*2*pi*q(j)*(Q-12+i)) )';
    end
end

Cb = 0;
%capacity bounds
for v = 1:P
   Cb = Cb+log2(1+(rho/P)*chi2rnd(2*v)) ;
end

%calculating capacity
for i=1:1000
    
% for m = 1:Q
%    for n = 1:P
%        Hvm(m,n) = sum(-1+(1+1)*rand(1,L)); %maximally rich scattering
%    end
% end

Hvac(Qc(1)+6:Qc(2)+6,Pc(1)+6:Pc(2)+6) = sum(-1+(1+1)*rand(1,L));
Hk=Hvac*ar*at';
C(i) = log2( det( eye(P,Q)+(rho/P)*(Hk.*Hk') ) )/sf;
end

mean(abs(C))
%Hk = Hva*ar*at';
C = abs(C);
C = sort(C/sf,'descend');
cap = linspace(45,64,1000);

figure()
plot(cap,C)
