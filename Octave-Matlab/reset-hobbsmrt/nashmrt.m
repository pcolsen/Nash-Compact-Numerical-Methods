% NASHMRT.M -- Nash Marquart in MATLAB
% This is explicitly for Hobbs problem.
% initiation
echo off
echo hobbsf off
hobbs0
b = (input('parameters in row form :'))';
disp(b)
[f, g, r, J, H] = hobbsf(b,y);
disp('Initial function value');
disp(f);
% check for non-computable function
if f < 0, 
  error('Non computable function at initial point');
%  break;
  
end;
% now set up iteration
E = 1000*ones(n,1);
lambda = 0.0001;
phi = 1;
% the offset vector
x = -1E+24*ones(n,1); 
JCALL=1;
% initial "best parameters" set so while does not stop first time through
while any((x+E)~=(b+E)),
  f0 = f;
  fprintf('After %g, Evals, function value = %g\n',JCALL, f);
  x = b;
  A = J'*J;
  g0 = g;
  while ((f>=f0) & (lambda>=0)),
% set lambda negative if we converge the parameters
    fprintf('LAMDA=%g\n',lambda);
    C = A + lambda.*(diag(diag(A)) + phi.*eye(n));
    d = -2.*C\g0; % we solve with built-in solver
    b = x + d;
%    fprintf('Trial b '); disp(b');
    if all((x+E)==(b+E)),
       lambda=-1;
       fprintf('Parameters not changed \n');
%      pause;
%  parameters not moved, so converged
     else
% parameters have changed, so evaluate function
      [f, g, r, J, H] = hobbsf(b,y);
      JCALL=JCALL+1;
      if f<0,
         if lambda<1E-20,
           lambda=eps;
         end;
         lambda=10*lambda;
         fprintf('NOT COMPUTABLE: new lambda =%g\n',lambda);
         f=f0+1E+20;
% ensure we have a big function to keep going
      else
         if f<f0,
           lambda=0.4*lambda;
           fprintf('SUCCESS: new lambda =%g\n',lambda);
%           pause;
% success
         else
           if lambda<1E-20,
             lambda=eps;
           end;
           lambda=10*lambda;
%           fprintf('FAILURE: new lambda =%g\n',lambda);
%           pause;
         end;   % of if f<f0
      end; % of if f<0
     end;
% of if (x+E)=(b+E) 
   end; 
% of while (f>=f0)
end; 
% while parameters not converged
disp('CONVERGED');
fprintf('Final parameters \n'); disp(x);
fprintf('Function value = %g\n',f0);
fprintf('gradient \n'); disp(g);
ss=input('Hit [cr] to continue ','s');
fprintf('Hessian  \n'); disp(H);
fprintf('H_inverse\n'); disp(inv(H));
fprintf(' Eigengvalues of H \n');
disp((eig(H))');
fprintf(' J''*J    \n'); disp(A);
fprintf('J''*J_inv \n'); disp(inv(J'*J));

