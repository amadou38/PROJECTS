function [X,W] = Quad2D(prec,T)
% AUTEUR : Diallo Amadou, 28/09/2020

[points, weight] = quadrature(prec);
A = T(:,1); B = T(:,2); C = T(:,3);
X = A + [B-A, C-A]*points';
X = X';
W = abs(det([B-A, C-A]))*weight;
end