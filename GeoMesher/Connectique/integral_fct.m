function [Integ] = integral_fct(f, Q, Xe, prec)
% AUTEUR : Diallo Amadou, 28/09/2020
N = size(Q,1);
points = zeros(3,2);  
points(3,:) = Xe;
Integ = 0;

for i = 1:N-1
   points(1:2,:) = Q(i:i+1,:);
   [X,W] = Quad2D(prec, points');
   Integ = Integ + sum(W.*f(X(:,1),X(:,2)));
end
points(1,:) = Q(1,:); points(2,:) = Q(N,:);
[X,W] = Quad2D(prec,points');
Integ = Integ + sum(W.*f(X(:,1),X(:,2)));
    
end