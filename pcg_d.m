function x = pcg_d(a,b,x,maxit,tol,U,Lambda)

n=length(b);
beta=0; 
r=b-a*x; ap=0*r;
gamma=r.'*r;
res0=norm(r);
res=res0;
mem=1;
k=0;
while (res/res0 > tol & k<maxit)

  z= r*U*(U'./diag(Lambda));
  k=k+1;
  gamma=r.'*z;
  if (k==1), p=z;else, beta=gamma/gamma0;p=z+beta*p;end
  ap=a*p;
  delta=p.'*ap;
  alfa = gamma/delta;
  x = x + alfa*p;
  r = r - alfa*ap;
  gamma0=gamma;
  res=norm(r);  
   mem=[mem,res/res0];
%disp([k,res,res/res0])

end	
 fprintf('CG its %d  rel.res %d\n',k,res/res0)

%semilogy(mem)
