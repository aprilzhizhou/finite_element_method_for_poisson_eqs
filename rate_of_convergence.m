%% problem setup 
% non-homogeneous dirichlet boundary condition 
u=@(x,y) -sin(pi*x).*cos(2*pi*y); % exact solution u(x,y)
f = @(X) -5*pi^2.*sin(pi*X(:,1)).*cos(2*pi*X(:,2)); % RHS function f(x,y)
gradu1 = @(x,y) -pi* cos(pi*x).*cos(2*pi*y); 
gradu2 = @(x,y) 2 *pi*sin(pi* x).* sin(2 *pi* y); % gradient of u(x,y)
%% generate mesh 
model = createpde(1);
% construct irregular geometry
% g = decsg(gd,sf,ns);
% geometryFromEdges(model,g);
% geometryFromEdges(model,@scatterg);
% geometryFromEdges(model,@squareg);
geometryFromEdges(model,@circleg);
% geometryFromEdges(model,@lshapeg);


hList = 1*0.5.^[1:5]; % hList contains different meshsizes for triangulation 
for kk = 1:length(hList)
    Hmax = hList(kk); 
MESH=generateMesh(model,'Hmax',Hmax); % use PDE toolbox to generate triangulation
[p,e,t] = meshToPet(model.Mesh); % export mesh data
node = p'; % nodes 
elem = t(1:3,:)'; % elements  
bdy = e(1,:); % boundary nodes 
bdy(2,:) = 0; % input inhomogeneous dirichlet boundary condition
for j = 1:length(bdy(1,:))
    i = bdy(1,j);
    bdy(2,j)=u(node(i,1),node(i,2));
end
M = size(elem,1); % number of finite elements
N = size(node,1); % numbe of nodes


%% construct global stiffness matrix
% construct lobal stiffness matrix & load vector
locS = zeros(3*M,3);locB = zeros(3,M);
for l = 1:M
    p = [node(elem(l,1),1),node(elem(l,2),1),node(elem(l,3),1);
        node(elem(l,1),2),node(elem(l,2),2),node(elem(l,3),2)]; 
    % p = [x1,x2,x3;y1,y2,y3]; (xi,yi), i = 1,2,3, are the coordinates of
    % three vertices of the triangle l 
    
    % CONSTRUCT LOCAL STIFFNESS
    % We adopt Long Chen's algorithm to compute local stiffness matrices
    Bnode = [p(1,1)-p(1,3),p(2,1)-p(2,3); p(1,2)-p(1,3),p(2,2)-p(2,3)];
    G = [[1,0]',[0,1]',[-1,-1]'];
    area = 0.5*abs(det(Bnode)); %  area of triangles  
    s = zeros(3,3);
    for i = 1:3
        for j = 1:3
            s(i,j) = area*((Bnode\G(:,i))'*(Bnode\G(:,j)));
        end
    end
    % store local stiffness matrices 
    locS(3*(l-1)+1:3*(l-1)+3,:) = s;
    
    % CONSTRUCT LOCAL LOAD VECTOR 
    % we use Gaussian quadratue to compute local load vector 
    mid1 = ( node(elem(l,2),:)+node(elem(l,3),:) )/2;
    mid2 = ( node(elem(l,3),:)+node(elem(l,1),:) )/2;
    mid3 = ( node(elem(l,1),:)+node(elem(l,2),:) )/2;
    b1 = area.* (f(mid2)+f(mid3))/6;
    b2 = area.* (f(mid3)+f(mid1))/6;
    b3 = area.* (f(mid1)+f(mid2))/6;
    % store local load vectors 
    locB(:,l) = [b1;b2;b3];
end
% assemble global stiffness matrix & load vector
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

%% implement non-homogeneous dirichlet boundary condition

Xbdy = bdy(1,:); 
Xint = setdiff(1:N,Xbdy); %separate boundary and interior nodes 
Uh = zeros(N,1);
Uh(Xbdy) = bdy(2,:)';
B = B - S*Uh; % implement nonzero derichlet BCs
S2 = S(Xint,Xint);
B2 = B(Xint);
Uh(Xint) = S2\B2; % use backslash operator to solve sparse linear system


%% plot solution & FEM approximation
Ue = zeros(N,1); % exact solution
for i = 1:N
    Ue(i) = u(node(i,1),node(i,2));
    % plot exact solution
end
% xlabel('x'); ylabel('y'); zlabel('u(x,y)')

%% H1 semi-norm 
localError1 = zeros(1,M);
localError2 = zeros(1,M);
for l = 1:M
    N1 = [node(elem(l,1),1),node(elem(l,1),2)]; % N1 = (x1,y1)
    N2 = [node(elem(l,2),1),node(elem(l,2),2)]; % N2 = (x2,y2)
    N3 = [node(elem(l,3),1),node(elem(l,3),2)]; % N3 = (x3,y3)
    A = det( [1,N1(1),N1(2); 1,N2(1),N2(2); 1,N3(1),N3(2)]  );
    gradPhi = 1/A*[-N3(2)+N2(2),N3(1)-N2(1); N3(2)-N1(2),-N3(1)+N1(1); -N2(2)+N1(2),N2(1)-N1(1)];
    % gradPhi are the gradients of phi1, phi2, and phi3.
    gradUh  = Uh(elem(l,1))*gradPhi(1,:) + Uh(elem(l,2))*gradPhi(2,:) + Uh(elem(l,3))*gradPhi(3,:);
    % gradUh is the gradiant of Uh.
    g1 = @(x,y) (gradu1(x,y)-gradUh(1))^2 + (gradu2(x,y)-gradUh(2))^2;
    % g(x,y) is the function (grad(U)(x,y) - grad(Uh))^2
    localError2(l) = Gauss(N1,N2,N3,g1);
    % Use Gaussian Quadrature to approximate the H1 norm of (grad(U)(x,y) - grad(Uh))^2
    p = [node(elem(l,1),1),node(elem(l,2),1),node(elem(l,3),1);
         node(elem(l,1),2),node(elem(l,2),2),node(elem(l,3),2)];
    xx = 1/3*sum(p(1,:)); yy = 1/3*sum(p(2,:));
    Phi = 1/(det([1,p(1,1),p(2,1); 1, p(1,2),p(2,2); 1,p(1,3),p(2,3)]))*[det([1,xx,yy; 1, p(1,2),p(2,2); 1,p(1,3),p(2,3)]);
        det([1, p(1,1),p(2,1);1,xx,yy; 1,p(1,3),p(2,3)]);
        det([1,p(1,1),p(2,1);1,p(1,2),p(2,2);1,xx,yy])];
    g2 =@(x,y) (u(x,y) - ( Uh(elem(l,1))*Phi(1)+ Uh(elem(l,2))*Phi(2)+Uh(elem(l,3))*Phi(3) ))^2;
    localError1(l) = Gauss(N1,N2,N3,g2);
end
globalError = sqrt(sum(localError2)+sum(localError1));
H1err(kk) = globalError; % global H1 error 

end


RoC = zeros(length(hList)-1,1); % Rate of Convergence
for i = 2:length(hList)
    RoC(i) = log(H1err(i-1)/H1err(i))/log(2);
end

T = table(hList',H1err',RoC,'VariableNames',{'h_max','H1Error','RateofConvergence'})

figure (1)
loglog(hList,H1err,'linewidth',2);
set(gca,'fontsize',15)
title('Loglog Plot','interpreter','latex','fontsize',18)
xlabel('mesh size $log(h)$','interpreter','latex','fontsize',20);
ylabel('$log(||u - u_h||_{H^1})$','interpreter','latex','fontsize',20);
grid on
Slope = abs(polyfit(log(hList),log(H1err),1));
fprintf('The slope of the loglog plot is approximately equal to %8.5f\n',Slope(1))

% CONSTRUCT GAUSSIAN QUADRATURE
% Gausian quadrature used to approxiamte integral of f(x,y) over triangles
% with vertices N1, N2, and N3.
function z = Gauss(N1,N2,N3,f)
A = 1/2*abs(det([N1(1),N2(1),N3(1); N1(2),N2(2),N3(2);1,1,1]));
% A is the area of the triangle with vertices N1 N2 N3
z = A*f((N1(1)+N2(1)+N3(1))/3,(N1(2)+N2(2)+N3(2))/3);
end
