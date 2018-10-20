clc
clear
close all

P = 11;
Q = 11;

rho = 100;
SNR = 10*log10(rho);

L = 200;

a = 0.5;


k = linspace(0,P-1,11);
kd = 0;
kf = k(end);
%Scaling factor
%for k = 0 gives diagonal approximation
%for k = P-1 gives full matrix
SF_d = sqrt(P^2./( P+kd.*(2*P-kd-1) ));
SF_f = sqrt(P^2./( P+kf.*(2*P-kf-1) ));

%Channel capacity
%C = log2( det( I+(p*Hv*Hhv)/P ) )

%x = sqrt(rho)*Hs+w;
Hva = zeros(P,Q);
BL = sum(-1+(1+1)*rand(1,L));
BLp = mean(abs(BL)^2);
Hva(4:8,4:8) = BLp;

% for i = 4:8
%    
%     for j = 4:8
%        Hva(i,j) = mean(abs( sum(-1+(1+1)*rand(1,L)) )^2);
%     end
%     
% end

x = linspace(-0.4,0.4,P);
y = linspace(-0.4,0.4,Q);

q = linspace(-0.5*(Q-1),0.5*(Q-1),Q);
thet_rq = q/Q;
p = linspace(-0.5*(P-1),0.5*(P-1),P);
thet_tp = p/P;

%setting At and Ar for equation (15)
At = zeros(1,P);
at = zeros(1,P)';
for i = 1:P
   at(i,1) = (P^-0.5)*exp(-i*2*pi*thet_tp(i)*(i-1)); 
end
for j = 1:P
for i = 1:P
   At(j,i) = (P^-0.5)*exp(-i*2*pi*thet_tp(j)*(i-1)); 
end
end

Ar = zeros(1,Q)';
ar = zeros(1,Q)';
for i = 1:Q
   ar(i,1) = (Q^-0.5)*exp(-i*2*pi*thet_rq(i)*(i-1)); 
end
for j = 1:Q
for i = 1:Q
   Ar(i,j) = (Q^-0.5)*exp(-i*2*pi*thet_rq(j)*(i-1)); 
end
end

%equation (15)
H = Ar*Hva*At;

%setting bounds for (28)
qk = max(p(1),p-0);
m = min(p(P),p+0);

% equation (28)
Hk = 0;
for i = 1:P
    for j = 1:P
Hk =Hk+ Hva(j,i)*ar(j)*at(i)'/SF_d;
    end
end

Hvk = zeros(P,Q);
Hvk(4:8,4:8) = Hk; 

% equation (56)
CH = log2( det( eye(P,Q)+(rho/P)*Hva*Hva' ) );

contour(x,y,Hva) 

%%
figure()
% for k = 1:P
%     CH(1,k) = log2( 1+(rho/P)*chi2rnd(2*k) );
% end

xch = linspace(44,62,11);

plot(xch,CH)

Sr = [pi/16,5*pi/16];
St = [-5*pi/16,-pi/16];

Qm = a*Q*sin(Sr(1));
Qp = a*Q*sin(Sr(2));

Pm = a*P*sin(St(1));
Pp = a*P*sin(St(2));
