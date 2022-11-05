function [e2neigh,Nod2e,e2edg,nedg,freeedg,freenod] = getIDs(inElt, vertices)
  
  nelt = length(inElt);
  nnod = size(vertices,1);
  e2nod = inElt;
  Nod2e = alg1(e2nod,nnod,nelt);
  e2neigh = alg2(e2nod,Nod2e,nelt);
  nedg = alg4(e2neigh,nelt);
  Verts = vertices;
  [e2edg,nedg,freeedg,freenod] = alg5(e2neigh,e2nod,Verts,nelt,nnod,nedg);
end
