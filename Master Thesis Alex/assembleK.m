function K = assembleK(Ke,nele,nnode)
% assembles the global stiffness matrix of bar
%   Ke: vector of element stiffness
%   nele: number of elements
%   nnode: number of nodes

K = sparse(nnode,nnode);

K = K + sparse(1:nele,1:nele,Ke,nnode,nnode);
K = K + sparse(1:nele,2:nele+1,-Ke,nnode,nnode);
K = K + sparse(2:nele+1,1:nele,-Ke,nnode,nnode);
K = K + sparse(2:nele+1,2:nele+1,Ke,nnode,nnode);


end