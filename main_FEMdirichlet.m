%% problem setup 
% non-homogeneous dirichlet boundary condition 
u=@(x,y) -sin(pi*x).*cos(2*pi*y); % exact solution u(x,y)
f = @(X) -5*pi^2.*sin(pi*X(:,1)).*cos(2*pi*X(:,2)); % RHS function f(x,y)

% we use function FEMdirichlet to solve Possion equation \nabla^2 u = f for
% u. We input exaction solution u(x,y) in the second argument for dirichlet
% boundary condition g(x,y).
FEMdirichlet(f,u,0.05)


