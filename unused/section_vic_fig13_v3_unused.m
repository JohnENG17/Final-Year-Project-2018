clc
clear
close all

P = 11;
Q = 11;

x = linspace(-0.4,0.4,P);
y = linspace(-0.4,0.4,Q);

L = 200;

a = 0.5;

Hva = zeros(P,Q);
BL = sum(-1+(1+1)*rand(1,L));

Hva(4,5) = BL; %1
Hva(7,7) = BL; %2
Hva(5,8) = BL; %3
Hva(6,4) = BL; %4

% for i = 0:2
%     for j = 0:2
% Hva(3+i,4+j) = sum(-1+(1+1)*rand(1,L)); %1
% Hva(6+i,6+j) = sum(-1+(1+1)*rand(1,L)); %2
% Hva(4+i,7+j) = sum(-1+(1+1)*rand(1,L)); %3
% Hva(5+i,3+j) = sum(-1+(1+1)*rand(1,L)); %4
%     end
% end

contour(x,y,Hva)