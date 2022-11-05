clear; clc; close all;
Verts = [0,0;1,0;1,1;0,1];
ne = length(Verts); np = 3;
Area = 0.5*abs(sum(Verts(:,1).*Verts([2:end,1],2) - Verts([2:end,1],1).*Verts(:,2)));
Xe = sum((Verts + Verts([2:end,1],:)) .* repmat(Verts(:,1).*Verts([2:end,1],2) - Verts([2:end,1],1).*Verts(:,2),1,2));
Xe = Xe/(6*Area);
he = 0;
for i = 1:ne-1
    for j = (i+1):ne
        he = max(he, norm(Verts(i, :) - Verts(j, :)));
    end
end
lambda = 1; nu = 1;
p = basis(Xe,he);
D = dof(np,p,ne,Verts,Xe,Area);
G = LHS(nu,lambda,np,p,ne,Verts,Xe,Area,he);
RHS = RHS_P(nu,lambda,np,p,ne,Verts,he);
P = G\RHS; G(1:3,:) = 0;
M0 = P'*G*P;
M1 = (eye(2*ne) - D*P)'*(eye(2*ne) - D*P);
alpha = 0.5*trace(M0);
M = M0 + alpha*M1
M0
M1
alpha = trace(M0)
G = LHS(nu,lambda,np,p,ne,Verts,Xe,Area,he)
RHS
D
P