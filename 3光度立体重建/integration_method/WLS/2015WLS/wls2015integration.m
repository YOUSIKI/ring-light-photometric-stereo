function [u,weights] = wls2015integration(p,q,mask,lambda,solver,gamma_int,u0,tol_precond,tol_linsolve,maxit)
%	Minimises : 
%	F(u) = \iint_{mask} w(x,y) ||\nabla u(x,y) - [p(x,y),q(x,y)]||^2 + \lambda (u(x)-u0(x))^2 dx dy
%
%	Required : 
%	p : estimate derivatives of u in the x-direction (towards bottom)
%	q : estimate derivatives of u in the y-direction (towards right)
%	
%	Optional : 
%	mask : binary mask (should at least be connected)	
%	lambda : fidelity weight (default is 1e-5)
%	solver : 'cholesky', 'gmres' or 'bicgstab' (default is cholesky)
%	gamma_int : weighting parameter for integrability (default is 10)
%	u0 : initial estimate (default is simchony(p,q))
%	tol_precond : droptolerance for incomplete LU preconditioning, for use with gmres/bicgstab (default is 1e-5)
%	tol_linsolve : stopping criterion for iterative linear solvers (default is 1e-6) 
%	maxit : stopping criterion for iterative linear solvers (default is 50)
%
%	Author : Yvain Queau - University of Toulouse - IRIT, UMR CNRS 5505
%	yvain.queau@enseeiht.fr
%
	
	if(nargin<2)
		disp('estimated gradient (p,q) must be provided')
		return;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ÌÝ¶È»¥»»fanhao
    temp = p; p=q; q = temp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (~exist('mask','var')|isempty(mask)) mask=ones(size(p));end
	if (~exist('lambda','var')|isempty(lambda)) lambda=1e-5;end
	if (~exist('solver','var')|isempty(solver)) solver='cholesky';end
	
	if (~exist('gamma_int','var')|isempty(gamma_int)) gamma_int=10;end	
	if (~exist('u0','var')|isempty(u0)) u0=simchony(p,q);end
	
	if (~exist('tol_precond','var')|isempty(tol_precond)) tol_precond=1e-3;end	
	if (~exist('tol_linsolve','var')|isempty(tol_linsolve)) tol_linsolve=1e-6;end	
	if (~exist('maxit','var')|isempty(maxit)) maxit=50;end	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calculate some usefuk masks : 								   %	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Omega = zeros(size(mask,1),size(mask,2),4);
	Omega_padded = padarray(mask,[1 1],0);
	Omega(:,:,1) = Omega_padded(3:end,2:end-1).*mask;	% In mask, neigbor "under" in mask 
	Omega(:,:,2) = Omega_padded(1:end-2,2:end-1).*mask; % In mask, neigbor "over" in mask 
	Omega(:,:,3) = Omega_padded(2:end-1,3:end).*mask;	% In mask, neigbor "right" in mask 
	Omega(:,:,4) = Omega_padded(2:end-1,1:end-2).*mask; % In mask, neigbor "left" in mask 
	clear Omega_padded 	
	clear uc_padded neighbors_control
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Compute integrability-based weightings:			 							   %	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
	p_y = 0.5*(p(:,[2:end end])-p(:,[1 1:end-1]));
	q_x = 0.5*(q([2:end end],:)-q([1 1:end-1],:));
	int = zeros(size(p));
	% p forward, q forward
	%~ p_y = p(:,[2:end end])-p;
	%~ q_x = q([2:end end],:)-q;
% 	int = int+abs(p_y-q_x);
% 	weights = exp(-gamma_int*int);
    weights= ones(size(p));
	clear p_y q_x

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calculate A and B such that the problem amounts to solve Au =B   %
	% A is the modified (weighted) Laplacian (5 points stencil) 	   %
	% B is the modified (weighted) Divergence of [p,q]				   %	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	indices_mask = find(mask>0);
	mapping_matrix = zeros(size(p));
	mapping_matrix(indices_mask)=1:length(indices_mask); 	% mapping_matrix(i,j) indicates the corresponding index of pixel (i,j) in the VECTOR u = u(find(mask>0))
	pbar = 0.5*(p+p([2:end end],:));
	qbar = 0.5*(q+q(:,[2:end end]));
	
	% Compute the system Ax = B, with size(A,1)=size(A,2)=length(B)=length(find(mask>0))	
	[A,B] = calculate_system(pbar,qbar,mask,Omega,lambda,u0,weights,mapping_matrix);
	clear Omega
	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	u = u0;	% Initialization for defining size	
% 	u(find(mask==0))=NaN;
	indices_mask = find(mask>0);
	
	% Solution 1 : normal equation by Cholesky
	if(strcmp(solver,'cholesky'))
		AtB=A'*B;
		AtA=A'*A;		
		spparms('spumoni',0)	
		sol_chol=AtA\AtB;
		u(indices_mask)=sol_chol;
	elseif(strcmp(solver,'gmres'))	
		[L,U] = ilu(A,struct('type','ilutp','droptol',tol_precond));		
		[sol,fl1,rr1,it1,rv1] = gmres(A,B,[],tol_linsolve,maxit,L,U,u(indices_mask));		
		u(indices_mask)=sol;
	elseif(strcmp(solver,'bicgstab'))	
		[L,U] = ilu(A,struct('type','ilutp','droptol',tol_precond));	
		[sol,fl1,rr1,it1,rv1] = bicgstab(A,B,tol_linsolve,maxit,L,U,u(indices_mask));	
		u(indices_mask)=sol;
	end

	%%% That's all folks ! %%%
	
end



function [A,B] = calculate_system(pbar,qbar,mask,Omega,lambda,u0,weights,mapping_matrix)

	I = [];
	J = [];
	K = [];
	indices_mask = find(mask>0);
	B = zeros(length(indices_mask),1);

	% All px in mask
	[X,Y]=find(mask>0);
	indices_center=sub2ind(size(mask),X,Y);
	I_centre = mapping_matrix(indices_center);
	A_centre = -lambda *ones(length(I_centre),1);	
	I=[I;I_centre]; % Which lines are we talking about ? 
	J=[J;I_centre]; % Which columns ?
	K=[K;A_centre]; % Which values ? 	
	B(I_centre)=-lambda*u0(indices_mask);	

	% In mask, right neighbor in mask
	set = Omega(:,:,3);
	[X,Y]=find(set>0);
	indices_centre=sub2ind(size(mask),X,Y);
	I_centre = mapping_matrix(indices_centre);
	indices_voisins=sub2ind(size(mask),X,Y+1);
	I_voisins = mapping_matrix(indices_voisins);
	A_centre = -(weights(indices_centre));
	A_voisin =weights(indices_centre);
	
	K=[K;A_centre(:);A_voisin(:)]; % Which values ?
	I=[I;I_centre(:);I_centre(:)]; % Which values ?
	J=[J;I_centre(:);I_voisins(:)]; % Which values ?
	B(I_centre) = B(I_centre)+weights(indices_centre).*qbar(indices_centre);
	
	
	% In mask, left neighbor in mask
	set = Omega(:,:,4);
	[X,Y]=find(set>0);
	indices_centre=sub2ind(size(mask),X,Y);
	I_centre = mapping_matrix(indices_centre);
	indices_voisins=sub2ind(size(mask),X,Y-1);
	I_voisins = mapping_matrix(indices_voisins);
	A_centre = -(weights(indices_voisins));
	A_voisin =weights(indices_voisins);
	K=[K;A_centre(:);A_voisin(:)]; % Which values ?
	I=[I;I_centre(:);I_centre(:)]; % Which values ?
	J=[J;I_centre(:);I_voisins(:)]; % Which values ?
	B(I_centre) = B(I_centre)-weights(indices_voisins).*qbar(indices_voisins);
	
	% In mask, top neighbor in mask
	set = Omega(:,:,2);
	[X,Y]=find(set>0);
	indices_centre=sub2ind(size(mask),X,Y);
	I_centre = mapping_matrix(indices_centre);
	indices_voisins=sub2ind(size(mask),X-1,Y);
	I_voisins = mapping_matrix(indices_voisins);
	A_centre = -(weights(indices_voisins));
	A_voisin =weights(indices_voisins);
	K=[K;A_centre(:);A_voisin(:)]; % Which values ?
	I=[I;I_centre(:);I_centre(:)]; % Which values ?
	J=[J;I_centre(:);I_voisins(:)]; % Which values ?
	B(I_centre) = B(I_centre)-weights(indices_voisins).*pbar(indices_voisins);	
		
	% In mask, bottom neighbor in mask
	set = Omega(:,:,1);
	[X,Y]=find(set>0);
	indices_centre=sub2ind(size(mask),X,Y);
	I_centre = mapping_matrix(indices_centre);
	indices_voisins=sub2ind(size(mask),X+1,Y);
	I_voisins = mapping_matrix(indices_voisins);
	A_centre = -(weights(indices_centre));
	A_voisin =weights(indices_centre);
	K=[K;A_centre(:);A_voisin(:)]; % Which values ?
	I=[I;I_centre(:);I_centre(:)]; % Which values ?
	J=[J;I_centre(:);I_voisins(:)]; % Which values ?
	B(I_centre) = B(I_centre)+weights(indices_centre).*pbar(indices_centre);	
	
	% Construction de A : 
	A = sparse(I,J,K);	
	
end


function [U,rmse] = simchony(p,q,gt)
% An implementation of the method from Simchony et Al for integration 
% of a normal field with natural boundary condition
% Ref : Direct Analytical Methods for Solving Poisson Equation in 
% Computer Vision problems - PAMI 1990
%
% example :
% p = ones(100,100); 
% q = ones(100,100); 
% u = simchony(p,q); 
%
% This performs the least square solution to \nabla u = [p,q], i.e. :
% min \int_\Omega \| \nablua U - [p,q] \|^2
% where \Omega is square and the natural Neumann boundary condition 
% \mu \cdotp (\nabla u -[p,q]) = 0 is used (2nd order derivatives 
% are neglected), where \mu is the outer 
% normal to \Omega. 
% This boundary condition is the one given by the 
% calculus of variations.
%
% % Axis : O->y
%        |
%        x
%
% The problem comes to solve a linear system Ax=b, where A is a bloc 
% Toplitz matrix. Fast solution is provided by Discrete Cosine Transform
%
% Implementation : Yvain Queau
% Universite de Toulouse, IRIT, UMR CNRS 5505
% yvain.queau@enseeiht.fr
% See other codes at http://ubee.enseeiht.fr/photometricstereo/

%~ p = 0.5*(p([2:end end],:)+p);
%~ q = 0.5*(q(:,[2:end end])+q);

% Compute div(p,q)
px = 0.5*(p([2:end end],:)-p([1 1:end-1],:));
qy = 0.5*(q(:,[2:end end])-q(:,[1 1:end-1]));

% Div(p,q) + Boundary Condition
f = px+qy;
f(1,2:end-1) = 0.5*(p(1,2:end-1)+p(2,2:end-1));
f(end,2:end-1) = 0.5*(-p(end,2:end-1)-p(end-1,2:end-1));
f(2:end-1,1) = 0.5*(q(2:end-1,1)+q(2:end-1,2));
f(2:end-1,end) = 0.5*(-q(2:end-1,end)-q(2:end-1,end-1));

f(1,1)=0.5*(p(1,1)+p(2,1)+q(1,1)+q(1,2));
f(end,1)=0.5*(-p(end,1)-p(end-1,1)+q(end,1)+q(end,2));
f(1,end)=0.5*(p(1,end)+p(2,end)-q(1,end)-q(1,end-1));
f(end,end)=0.5*(-p(end,end)-p(end-1,end)-q(end,end)-q(end,end-1));

% Sine transform of f
fsin=dct2(f);

% Denominator
[x,y] = meshgrid(0:size(p,2)-1,0:size(p,1)-1);
denom = (2*cos(pi*x/(size(p,2)))-2) + (2*cos(pi*y/(size(p,1))) - 2);
Z = fsin./(denom);
Z(1,1)=0.5*Z(1,2)+0.5*Z(2,1); %Or whatever...

% Inverse Sine transform :
U=idct2(Z);


	if(nargin>2)% Ground truth available
		moyenne_ecarts=mean(U(:)-gt(:));
		U=U-moyenne_ecarts;
		npix=size(p,1)*size(p,2);
		rmse=sqrt((sum((U(:)-gt(:)).^2))/npix);
	else
		U=U-min(U(:));
	end

end



