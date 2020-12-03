%discretisation of Burgers' eqn using AJR's holistic approach to order
%gamma^2. Variable diffusion of epsilon. 
%Implements Roberts, Applied Numerical Mathematics 2001.

%input:numx - number of lattice points (excluding boundaries).
%h - mesh spacing
%gamma - holistic coupling coefficient (set to 1 to simulate full pde)
%epsilon - diffusivity coefficient for burgers' pde
function f = burgDisc(numx,h,gamma,epsilon)
%%Set up difference operators.
ovec = ones(numx-1,1);
D = -2*eye(numx)+diag(ovec,1)+diag(ovec,-1); %D/h^2 evaluates the centred second deriv, truncated near bnds
D2 = D^2; %centred fourth order deriv
C = (diag(ovec,1) + diag(-ovec,-1))/2; %C/h is the centred first deriv, truncated near bnds


g=gamma; g2 = gamma^2; %gamma=1 for discretization; uses AJR holistic approach

%set up functions to be used at both boundaries
%a(1) is a boundary condition. a(2) is h^2 * da/dt
%diffusion terms
funD =@(a) [g+g2/6,    -g/12-g2/45;...
    -g2/12,  g2/90] *a;

funC = @(a,u) [g2*u(1)/2,  -g2*u(1)/24; 0 0]*a;

%evaluate at boundaries
%left boundary
lBnd =@(t,u,al) epsilon/h^2*funD(al) + 1/h*funC(al,u);

%right boundary
rBnd = @(t,u,ar) flipud( epsilon/h^2*funD(ar) - 1/h*funC(ar,flipud(u)) );

%global
bnd = @(t,u,al,ar) [lBnd(t,u,al); zeros(numx-4,1); rBnd(t,u,ar)];


f=@(t,u,al,ar) g*epsilon/h^2 *  D*u ... 
    -g2/h * (C*u).*u ...
    -g2*epsilon/12/h^2 * D2*u...
    +bnd(t,u,al,ar);
