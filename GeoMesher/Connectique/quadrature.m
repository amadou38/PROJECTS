function [point, weight] = quadrature(i)
% AUTEUR : Diallo Amadou, 28/09/2020

point = zeros(i,2);
weight = zeros(i,1);
if ( i == 1) 
    a = 1/3;
    point = [a a];
    weight(1) = 1/2;
elseif ( i == 3)
    a = 2/3;
    b = 1/6;
    w = 1/6;
    point = [b b; a b; b a];
    weight = w*ones(3,1);
elseif ( i == 4)
    a = 1/5;
    b = 3/5;
    c = 1/3;
    w1 = 25/(24*4);
    w2 = -27/(24*4);
    point = [a a; b a; a b; c c];
    weight(1:i-1) = w1; weight(i) = w2; 
end








end