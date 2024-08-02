function u = febar( D, q, L, nele )
% returns the nodal displacements of a bar subjected to uniform axial loading
%  D: vector with axial resistances at each element
%  q: uniformly distributed load
%  L: length of bar
%  nele: number of elements

% number of nodes
nnode = nele+1;

% element length
l = L/nele;

% consistent nodal forces
F = zeros(nnode,size(q,2));

% assemble force vector
F(1:end-1,:) = l/2*q; 
F(2:end,:) = F(2:end,:) + l/2*q; 

% assemble stiffness matrix
Ke = D/l;
K = assembleK(Ke,nele,nnode);

% reduce stiffness matrix and force vector to apply fixed end boundary
% condition

Kred = K(2:nnode,2:nnode);
Fred = F(2:nnode,:);

% solution for nodal displacements
u = Kred\Fred; % This operator skip the inversing on the matrix, saving computational cost

u = [zeros(1,size(q,2));u];


end