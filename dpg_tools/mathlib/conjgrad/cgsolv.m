function  [result,P,R_cg] = cgsolv( A, b, x, k_max, verbose)
% CGSOLV   solve Ax = b using the Conjugate Gradient algorithm
%
% [result,P,R_cg] = cgsolv( A, b, x, k_max, verbose)
%
% taken from p.529 of Golub and Van Loan, "Matrix Computations," 3rd ed.
% S. 10.2 The Conjugate Gradient Method
% S. 10.2.6 Some Practical Details.
%
% See also CMPR, PRACT_SD, PRACT_CG

if nargin < 3, x = randn(size(A,2),1); end;
if nargin < 4, k_max = 20; end;
if nargin < 5, verbose = 0; end;

k = 0;
r = b - A*x;

p = 0;
xi = x;

% clear P R_cg X_cg A_cg;
if verbose, X_cg(:,1) = x; end;

while ( k < k_max ) 

  rho1 = sqrt(r'*r);         
  if verbose; R_cg(k+1) = rho1; end;
  k = k+1;
  if rho1 < eps*prod(size(b)), break; end;
  if k == 1, 
    beta = 0; 
  else, 
    beta = (rho1 / rho0)^2 ; 
  end;
  p = r + beta*p;
    
  w = A*p;                              % linear operation on p
  
  alpha = rho1^2 / (p'*w);
  x = x + alpha*p;
  if verbose, 
    X_cg(:,k+1) = x; A_cg(:,k) = alpha; P(:,k) = r; 
  end;
  
  rho0 = rho1;
  r = r - alpha*w;
end;

if verbose; keyboard; end;
result = x;
% Copyright 2001 William Scott Hoge (shoge@ece.neu.edu or shoge@ieee.org) 
% All rights reserved. Licensed according to the GNU General Public Licence
% (see http://www.gnu.org)

