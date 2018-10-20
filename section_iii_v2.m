%John Tran 25999001 FYP 2018
%Contour plot of |H_v(q,p)|
clc
clear
close all
%%
%x-axis
P = 21;
Pb = (P-1)/2;
p = linspace(-Pb,Pb,P);

%y-axis
Q = 21;
Qb = (Q-1)/2;
q = linspace(-Qb,Qb,Q);

a = 0.5;
L = 200;

%[DELETE]
Srq = linspace( -1/(2*Q),1/(2*Q),200 );
Stp = linspace( -1/(2*P),1/(2*P),200 );

%% a
figure()
%phi_t @ (pi/8,-pi/8)
%phi_r @ (-pi/8,pi/8)

%Angular spreads for top left point scatterer
Sr = pi/8;
St = -pi/8;

%finding location of the top left point scatterers
Pa_c = P*a*sin(St);
Pa_c = round(Pa_c,0); 
Qa_c = Q*a*sin(Sr);
Qa_c = round(Qa_c,0); 

x = linspace(-0.4,0.4,P);
y = linspace(-0.4,0.4,Q);

z = zeros(P,Q);
%top left point scatterer 
Hva = z;
Hva2 = z;
Hva(11+Qa_c,11+Pa_c) = abs(sum(-1+(1+1)*rand(1,L)));
%bottom right point scatterer 
Hva2(11+Pa_c,11+Qa_c) = abs(sum(-1+(1+1)*rand(1,L)));

contour(x,y,Hva) 
hold on
contour(x,y,Hva2) 
xlabel('Transmit Angle')
ylabel('Receive Angle')
title('|H_V(q,p)|: Point Scatterers')

%% b
figure()
%angular spreads for top left
Srb = [pi/16,3*pi/16];
Stb = [-3*pi/16,-pi/16];

%Finding cluster size of top left
Qb_c = zeros(1,2);
Pb_c = zeros(1,2);

Qb_c(1) = a*Q*sin(Srb(1));
Qb_c(2) = a*Q*sin(Srb(2));
Qb_c = round(Qb_c,0);
Pb_c(1) = a*P*sin(Stb(1));
Pb_c(2) = a*P*sin(Stb(2));
Pb_c = round(Pb_c,0);

%cluster size [delete]
qb = linspace(Qb_c(1),Qb_c(2),5); %Q- = 2, Q+ = 6
pb = linspace(Pb_c(1),Pb_c(2),5); %P- = -6, P+ = -2

%rank
r = min( qb(end)-qb(1)+1,pb(end)-pb(1)+1 ); %rank of submatrix

Hvb = z;

%[delete]
% for i= 1:2
% if Qb_c(i) >= 11
%     Qb_c(i) = 10;
% end
% if Pb_c(i) <=-11
%     Pb_c(i) = -10;
% end
% end

%bottom right cluster
for i = Pb_c(1)+11:Pb_c(2)+11
   for j = Qb_c(1)+11:Qb_c(2)+11
      Hvb(i,j) = sum(-1+(1+1)*rand(1,L)); 
   end
end

%top left cluster
for i = Pb_c(1)+11:Pb_c(2)+11
   for j = Qb_c(1)+11:Qb_c(2)+11
      Hvb(j,i) = sum(-1+(1+1)*rand(1,L)); 
   end
end

contour(x,y,abs(Hvb)) 
xlabel('Transmit Angle')
ylabel('Receive Angle')
title('|H_V(q,p)|: Clustered Scatterers')
%rankb = rank(Hvb);

%% c
figure()
Hvc = zeros(21,21);
c = 1;

for i = 0:P-1
    Hvc(P-i,c) = sum(-1+(1+1)*rand(1,L)); 
    c = c + 1;
end

% Hvc(1,20) = Hvc(1,21);
% Hvc(2,21) = Hvc(1,21);
% Hvc(20,1) = Hvc(21,1);
% Hvc(21,2) = Hvc(21,1);

% Hvc(1,20) = sum(-1+(1+1)*rand(1,L)); 
% Hvc(2,21) = sum(-1+(1+1)*rand(1,L)); 
% Hvc(20,1) = sum(-1+(1+1)*rand(1,L)); 
% Hvc(21,2) = sum(-1+(1+1)*rand(1,L)); 

contour(x,y,abs(Hvc)) 
xlabel('Transmit Angle')
ylabel('Receive Angle')
title('|H_V(q,p)|: Diagonal Scatterering')
%% d
figure()
Hvd = zeros(21,21);

d = 0;
for i = 1:Q
   for j = 1:P
       Hvd(i,j) = sum(-1+(1+1)*rand(1,L));
   end
end

contour(x,y,abs(Hvd)) 
xlabel('Transmit Angle')
ylabel('Receive Angle')
title('|H_V(q,p)|: Maximally Rich Scatterering')







%% working out

%Transmit
phi_tx = pi/8;
phi_ty = -pi/8;
theta_tx =a*sin(phi_tx); 
theta_ty =a*sin(phi_ty); 

%Receive
phi_rx =-pi/8;
phi_ry = pi/8;
theta_rx =a*sin(phi_rx); 
theta_ry =a*sin(phi_ry); 

% N = 21;
% D = pi/8;
% Value = 0;
% H = phased.ULA(N,D);