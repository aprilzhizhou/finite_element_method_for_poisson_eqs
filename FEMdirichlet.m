function [ ] = FEM_poisson( f,g, Hmax)
% f(X) - RHS function, X = [x,y];
% g(x,y) -  derichlet boundary condition;
% Hmax - meshsize 
%% generate mesh 
model = createpde(1); 
geometryFromEdges(model,@circleg);
MESH=generateMesh(model,'Hmax',Hmax); % use PDE toolbox to generate triangulation
pdeplot(model)
set(gca,'fontsize',14)
axis equal

[p,e,t] = meshToPet(model.Mesh); % export mesh data
node = p'; % nodes 
elem = t(1:3,:)'; % elements  
bdy = e(1,:); % boundary nodes 
bdy(2,:) = 0; % input homogeneous dirichlet boundary condition
for j = 1:length(bdy(1,:))
    i = bdy(1,j);
    bdy(2,j)=g(node(i,1),node(i,2));
end
M = size(elem,1); % number of finite elements
N = size(node,1); % numbe of nodes

%% construct global stiffness matrix
% construct lobal stiffness matrix
locS = zeros(3*M,3);locB = zeros(3,M);
for l = 1:M
    p = [node(elem(l,1),1),node(elem(l,2),1),node(elem(l,3),1);
        node(elem(l,1),2),node(elem(l,2),2),node(elem(l,3),2)];
        % p = [x1,x2,x3;y1,y2,y3]; (xi,yi), i = 1,2,3, are the coordinates of
        % three vertices of the triangle l 
    
    % CONSTRUCT LOCAL STIFFNESS
    % we adopt Long Chen's algorithm to compute local stiffness matrices
    Bnode = [p(1,1)-p(1,3),p(2,1)-p(2,3); p(1,2)-p(1,3),p(2,2)-p(2,3)];
    G = [[1,0]',[0,1]',[-1,-1]'];
    area = 0.5*abs(det(Bnode));
    s = zeros(3,3);
    for i = 1:3
        for j = 1:3
            s(i,j) = area*((Bnode\G(:,i))'*(Bnode\G(:,j)));
        end
    end
    % store local stiffness matrices 
    locS(3*(l-1)+1:3*(l-1)+3,:) = s;
    
    % we use Gaussian quadrature to compute local load vectors 
    mid1 = ( node(elem(l,2),:)+node(elem(l,3),:) )/2;
    mid2 = ( node(elem(l,3),:)+node(elem(l,1),:) )/2;
    mid3 = ( node(elem(l,1),:)+node(elem(l,2),:) )/2;
    b1 = area.* (f(mid2)+f(mid3))/6;
    b2 = area.* (f(mid3)+f(mid1))/6;
    b3 = area.* (f(mid1)+f(mid2))/6;
    % store local stiffness matrices 
    locB(:,l) = [b1;b2;b3];
end
% assemble global stiffness matrix
S= sparse(N,N); B =zeros(N,1);
for l = 1:M
    for i = 1:3
        for j = 1:3
            S(elem(l,i),elem(l,j)) = S(elem(l,i),elem(l,j)) + locS((l-1)*3+i,j);
            % global assembly of global stiffness matrix 
        end
    end
    for k = 1:3
            B(elem(l,k)) = B(elem(l,k)) + locB(k,l);
            % global assembly of global load vector
    end 
end

%% implement dirichlet boundary condition
Xbdy = bdy(1,:);
Xint = setdiff(1:N,Xbdy); %separate boundary and interior nodes 
Uh = zeros(N,1);
Uh(Xbdy) = bdy(2,:)';
B = B - S*Uh; % insert nonzero dirichlet BCs
S2 = S(Xint,Xint);
B2 = B(Xint);
Uh(Xint) = S2\B2; % use backslash operator to solve sparse linear system

%% plot solution & FEM approximation
figure (1)
trisurf(elem,node(:,1)',node(:,2)',Uh); % plot numerical solution
title('FEM Approximation')
xlabel('x'); ylabel('y'); zlabel('u_h(x,y)')
set(gca,'fontsize',14)
axis equal
end

