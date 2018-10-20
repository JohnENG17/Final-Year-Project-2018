clc
clear
close all

%% to graph contour plot, need three vectors
% p, q and Hv(q,p)

Q = 21;
Q_bar = (Q-1)/2;
%receive angle - y-axis
q = linspace(-Q_bar,Q_bar,Q);
P = 21;
P_bar = (Q-1)/2;
%transmit angle - x-axis
p = linspace(-P_bar,P_bar,P); 

alpha_T = 0.5;
alpha_R = 0.5;
L = 200; %number of paths

S_rq = linspace( -1/(2*Q),1/(2*Q) );
S_tp = linspace( -1/(2*P),1/(2*P) );

%% part a
%==========================================================================
% %Bottom right cluster?
% phi_t1 = pi/8;
% phi_r1 = -pi/8;
% phi_1 = linspace(phi_t1,phi_r1);
% 
% %Top left cluster?
% phi_t2 =-pi/8;
% phi_r2 = pi/8;
% phi_2 = linspace(phi_t2,phi_r2);

% p_t = pi/8; q_t = -pi/8;
% phi_T = linspace(p_t,q_t);
% 
% p_r = -pi/8; q_r = pi/8;
% phi_R = linspace(p_r,q_r);
% 
% theta_r = alpha_R*sin(phi_R);
% theta_t = alpha_T*sin(phi_T);

% qq = q_t/Q
% pp = p_t/P

% syms k
% 
% G_b = symsum(dirac((q_t)/Q-mod(k,1))*dirac((p_t)/P-mod(k,1)), k, 1,L)

%H_v = G_b/(P*Q);
%contour(p,q,H_v);

% l = 1;
% G = @(tr,tt) dirac(tr-theta_r(l)).*dirac(tt-theta_t(l));
% 
% Ga = 0;
% for l = 1:L
%    
%     %Ga = Ga + G(theta_r(L),theta_t(L));
%     Ga = Ga + G(q,p);
%     idx = Ga == Inf;
%     Ga(idx) = 1;
%     
% end
% 
% Hv = Ga/(P*Q);


% G_bar = @(theta_r,theta_t) G(mod(theta_r,1),mod(theta_t,1));
% G_bar = @(theta_r,theta_t) dirac(theta_r-q/Q)*dirac(theta_t-p/P);
% H_v = @(q,p) G_bar(q./Q,p./P)./(P.*Q);
% 
% H_va = abs(H_v(q,p));

% for qi = -Q_bar:Q_bar 
%     for pi = -P_bar:P_bar
%     G_bar = dirac(theta_r-qi/Q).*dirac(theta_t-pi/P);
%     end
% end
%==========================================================================

Bl = rand(1,L);

%% part b
n = 21; %might refer to P = Q = 21;

%Angular spreads for the top left cluster
Sr_b = linspace(pi/16,3*pi/16);
St_b = linspace(-3*pi,-pi/16);
% Sr_b = [pi/16,3*pi/16]; %receive angle
% St_b = [-3*pi/16,-pi/16]; %transmit angle

Q_neg = alpha_R*Q*sin(Sr_b(1));
Q_pos = alpha_R*Q*sin(Sr_b(end));
q_b = linspace(Q_neg,Q_pos);

P_neg = alpha_T*P*sin(St_b(1));
P_pos = alpha_T*P*sin(St_b(end));
p_b = linspace(P_neg,P_pos);

r = min(Q_pos - Q_neg + 1, P_pos - P_neg + 1);

%(Q-,Q+) = (2,6)
%(P-,P+) = (-6,-2)
%r = 5
subm = ones(5,5);
%% part c
% kc = 0;
Sr_c = [-pi/2,pi/2];
St_c = [-pi/2,pi/2];
 
%% part d
% kd = P-1;

%% working out
x = linspace(-0.4,0.4);
y = linspace(-0.4,0.4);
z = [9 6 1;
    6 0 6;
    1 6 9];
z = ones(100,100);
z2 = z;

%Transmit
phi_tx = pi/8;
phi_ty = -pi/8;
theta_tx =alpha_T*sin(phi_tx); 
theta_ty =alpha_T*sin(phi_ty); 

%Receive
phi_rx =-pi/8;
phi_ry = pi/8;
theta_rx =alpha_R*sin(phi_rx); 
theta_ry =alpha_R*sin(phi_ry); 

[Y1,X1] = find(x==0.1913);
[Y2,X2] = find(y==-0.1913);

z(75,26) = rand(1);
z2(26,75) = rand(1);

contour(x,y,z) 
hold on
contour(x,y,z2)
xlabel('Transmit Angle')
ylabel('Receive Angle')


N = 21;
D = pi/8;
Value = 0;
H = phased.ULA(N,D);











