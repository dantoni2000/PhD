function [sz,Z]=Nystrom_RK_Lyap(A,B1,B2,params,X)

[n,~]=size(A);
lN_0=10; lN_max=n/2; qN=5; tolN_e=1e-4; tolN_r=1e-4; miN=1;
[UN,Lhat] = Adaptive_Nystrom(-A,lN_0,lN_max,qN,tolN_e,tolN_r,miN);
[l,~]=size(Lhat);
%[UN,Lhat] = Nystrom(-A,l);

tic
   m=params.m;
   tol=params.tol;
   s1=params.smin;
   emax=params.smax;
   ch=params.ch;
   period=params.period;
   
   dim1 = [];

[n,n]=size(A);
B1=full(B1);
B2=full(B2);
p=size(B1,2);
I=speye(p);O=0*I;
In=speye(n);
Il=speye(l);

csize = size(B1,2);

[V1,rr1]=qr(B1,0); 
[V2,rr2]=qr(B2,0);
nrmb=norm((rr1),'fro')*norm((rr2),'fro'); 
beta1=V1'*B1; Beta1=beta1*beta1';
beta2=V2'*B2; Beta2=beta2*beta2';
errtot=[];

VV1=V1;
VV2=V2;

H1=sparse(p*(m+2),p*(m+1));
H2=sparse(p*(m+2),p*(m+1));
nrmrestotnew=[];
%nrma=norm(A,'fro');

if (norm(A-A',1)<1e-14), symm=1;
fprintf('The matrix A is symmetric \n')
else symm=0;
fprintf('The matrix A is nonsymmetric \n')
end

fprintf('\n')

fprintf('     # iter      rel.res.\n')

 newAv1=A*V1;
 newAv2=A*V2;
 %newAv=A'*V;
 K1=full(V1'*newAv1); 
 K2=full(V2'*newAv2);
 s=s1(1);
 z=s1(1);
 eH1=eig(K1);
 eH2=eig(K2);
 eH1points = sort([s1(:)',emax]);
 eH2points = sort([s1(:)',emax]);
 snew=newpolei(eH1points,eH1,s1(1)*ones(p,1));
 znew=newpolei(eH2points,eH2,s1(1)*ones(p,1));
 if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end
 if real(znew)<0, znew=-real(znew)+sqrt(-1)*imag(znew);end
 s=[s,snew];
 z=[z,znew];
% additional steps
cmplxflag=0;
itsinner=0;


i=0;

itcheck = i;
while i < m


  i=i+1;

  paired=0;
  itp=1;
  c1flag = 0;
  c2flag = 0;
  V1wrk = V1;
  V2wrk = V2;

  while (paired==0),

    i1=i+1; it = 0; t=0.;
    w1=V1wrk;
    w2=V2wrk;
 
    
     its_cg=0; res_cg=0;
     ANyinv_s = UN*(UN'./diag(-sparse(Lhat)-snew*Il)) - 1/snew*(In-UN*UN');
     snew
     aaaaa = norm ((A-snew*In)^-1 - ANyinv_s)
     % ANyinv_s = (A-snew*In)^-1;
     ANyinv_z = UN*(UN'./diag(-sparse(Lhat)-znew*Il)) - 1/znew*(In-UN*UN');
     % ANyinv_z = (A-znew*In)^-1;

     w1rk1 = ANyinv_s*w1; 
     w2rk1 = ANyinv_z*w2;
     
     %%%% All real basis implementation for RKSM %%%%%
     
     if imag(w1rk1) ~= 0 & c1flag == 0
         w1rk = real(w1rk1);
         c1flag = 1;
     elseif imag(w1rk1) ~= 0 & c1flag == 1
         w1rk = imag(w1rk1);
     else
         w1rk = w1rk1;
     end

     if imag(w2rk1) ~= 0 & c2flag == 0
         w2rk = real(w2rk1);
         c2flag = 1;
     elseif imag(w2rk1) ~= 0 & c2flag == 1
         w2rk = imag(w2rk1);
     else
         w2rk = w2rk1;
     end
     
     
% Gram-Schmidt step
    jms=(i-1)*p+1;j1s=(i+1)*p;js=i*p;js1=js+1;
    B1m(jms:js,1:csize)= V1'*B1;
    B2m(jms:js,1:csize)= V2'*B2;

    for it=1:2,
      for kk=1:i
        k1=(kk-1)*p+1; k2=kk*p; 
        ww1=VV1(1:n,k1:k2);
        ww2=VV2(1:n,k1:k2);
        gamma1=ww1'*w1rk;
        gamma2=ww2'*w2rk;
        H1(k1:k2,jms:js) = H1(k1:k2,jms:js)+ gamma1;
        H2(k1:k2,jms:js) = H2(k1:k2,jms:js)+ gamma2;
        w1rk = w1rk - ww1*gamma1;
        w2rk = w2rk - ww2*gamma2;
      end
    end
    
    [V1,h1inv]=qr(w1rk,0); H1(js1:j1s,jms:js)=h1inv; %hinv = inv(hinv);
    [V2,h2inv]=qr(w2rk,0); H2(js1:j1s,jms:js)=h2inv; %hinv = inv(hinv);

    if (cmplxflag), 
        snew=conj(snew); s=[s,snew];
        znew=conj(znew); z=[z,znew];cmplxflag=0;
        newAv1=A*V1;
        newAv2=A*V2;
        %newAv=A'*V;
        D1 = kron(spdiags(s(2:end)),I);
        D2 = kron(spdiags(z(2:end)),I);
        g = VV1'*newAv1;
        g1 = g; 
        g2 = V1'*A*VV1; g3 = V1'*A*V1;
        %g2 = V'*A'*VV; g3 = V'*A'*V;
        K1 = [K1 g1; g2, g3];
        VV1=[VV1,V1];

        q = VV2'*newAv2;
        q1 = q; 
        q2 = V2'*A*VV2; q3 = V2'*A*V2;
        %q2 = V2'*A'*VV2; q3 = V2'*A'*V2;
        K2 = [K2 q1; q2, q3];
        VV2=[VV2,V2];
        i=i+1; itp=itp+1;
    else, 
        paired=1; 
    end
  end

  VV1new = [VV1, V1];
  VV2new = [VV2, V2];
  ih1=i1; ih=i;
  newAv1=A*V1;
  newAv2=A*V2;
    %newAv=A'*V;
  D1 = kron(spdiags(s(2:end)),I);
  D2 = kron(spdiags(z(2:end)),I);
  g = VV1'*newAv1;
  q = VV2'*newAv2;
    
    
  if (symm), K1=(K1+K1')/2; K2=(K2+K2')/2; end

  if (rem(i,period)==0)

% Solve the projected problem
      Y=lyap(K1,K2,B1m*B2m'); 

      err(i)=norm(X - VV1(:,1:size(Y,1))*Y*VV2(:,1:size(Y,1))');
      % 
      % ee=X - VV1(:,1:size(Y,1))*Y*VV2(:,1:size(Y,1))';
      % energy(i)=trace(ee'*A*ee + ee'*ee*A);
      
% computed residual   (exact, in exact arithmetic) cheaper computation possible
     u1=newAv1-VV1*g;
     d1=-VV1*(Y*(H1(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H1(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
     U1=full([-V1*s(end),  d1 u1 ]);
     rr1=qr(U1,0); rr1=triu(rr1(1:size(rr1,2),:));

     u2=newAv2-VV2*q;
     d2=-VV2*(Y*(H2(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H2(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
     U2=full([-V2*z(end),  d2 u2 ]);
     rr2=qr(U2,0); rr2=triu(rr2(1:size(rr2,2),:));
     

%backward error
     nrmres = norm(rr1*sparse([O I O; I O I; O I O ])*rr2','fro');
     % anrmres = norm(A*VV1(:,1:size(Y,1))*Y*VV2(:,1:size(Y,1))'+VV1(:,1:size(Y,1))*Y*VV2(:,1:size(Y,1))'*A + B1*B2','fro');
     % why=anrmres-nrmres
% relative residual norm
     nrmresnew = (nrmres)/nrmb;
 
     nrmrestotnew = [nrmrestotnew, nrmresnew];

     dim = size(VV1,2);
     dim1 = [dim1,dim];

     disp([i,nrmresnew])

     if (nrmresnew<tol), 
        break
     end
  end 

% New poles and zeros
  eH1=sort(eig(K1));
  eH1orig=eH1;
  eH2=sort(eig(K2));
  eH2orig=eH2;
  %eK=sort(eig(full(K)));

  if (ch)                     % Complex poles. Compute set for next complex pole of r_m

    if (any(imag(eH1)) ~=0 & max(abs(imag(eH1)))>1e-5 & length(eH1)>2) % Roots lambdas come from convex hull too
     eH1=[eH1;-emax];
      ij=convhull(real(eH1),imag(eH1)); eH1=eH1(ij);
      ieH1=length(eH1); missing1=ih*p-ieH1;
      while missing1>0,                         % include enough points from the border
        neweH1=(eH1(1:ieH1-1)+eH1(2:ieH1))/2;missing1=ih*p-length(eH1);
        eH1=[eH1;neweH1];
      end
    % eH1=eH1(1:ih);
      eH1points=-eH1;
      eH1=eH1orig;
    else                                  % if all real eigs, no convex hull possible
      eH1points = sort([s1; emax.';-real(eH1)]);
    end

    if (any(imag(eH2)) ~=0 & max(abs(imag(eH2)))>1e-5 & length(eH2)>2) % Roots lambdas come from convex hull too
     eH2=[eH2;-emax];
      iij=convhull(real(eH2),imag(eH2)); eH2=eH2(iij);
      ieH2=length(eH2); missing2=ih*p-ieH2;
      while missing2>0,                         % include enough points from the border
        neweH2=(eH2(1:ieH2-1)+eH2(2:ieH2))/2;missing2=ih*p-length(eH2);
        eH2=[eH2;neweH2];
      end
    % eH2=eH2(1:ih);
      eH2points=-eH2;
      eH2=eH2orig;
    else                                  % if all real eigs, no convex hull possible
      eH2points = sort([s1; emax.';-real(eH2)]);
    end


  else   % Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
     
     if (any(imag(eH1)) ~=0 & length(eH1)>2)    % Roots lambdas come from convex hull too
       eH1=[eH1;-s1;-emax.'];
       ij=convhull(real(eH1),imag(eH1)); eH1=eH1(ij);
       ieH1=length(eH1); missing1=ih*p-ieH1;
       while missing1>0, % include enough points from the border
         neweH1=(eH1(1:ieH1-1)+eH1(2:ieH1))/2;
         eH1=[eH1;neweH1];
         missing1=ih*p-length(eH1);
       end
       eH1=eH1(1:ih*p);
     end
      eH1points = sort([s1; emax.';-real(eH1)]);
      eH1=eH1orig;

     if (any(imag(eH2)) ~=0 & length(eH2)>2)    % Roots lambdas come from convex hull too
       eH2=[eH2;-s1;-emax.'];
       ij=convhull(real(eH2),imag(eH2)); eH2=eH2(ij);
       ieH2=length(eH2); missing2=ih*p-ieH2;
       while missing2>0, % include enough points from the border
         neweH2=(eH2(1:ieH2-1)+eH2(2:ieH2))/2;
         eH2=[eH2;neweH2];
         missing2=ih*p-length(eH2);
       end
       eH2=eH2(1:ih*p);
     end
      eH2points = sort([s1; emax.';-real(eH2)]);
      eH2=eH2orig;
  end


  gs=kron(s(2:end),ones(1,p))';
  qs=kron(z(2:end),ones(1,p))';
  snew = newpolei(eH1points,eH1,gs);
  znew = newpolei(eH2points,eH2,qs);
  if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end
  if real(znew)<0, znew=-real(znew)+sqrt(-1)*imag(znew);end

% If pole is complex, include its conjugate

  if (imag(snew) ~=0), cmplxflag=1;end
  s=[s,snew];

  g1 = g; 
  g2 = V1'*A*VV1; g3 = V1'*A*V1;
  %g2 = V1'*A'*VV1; g3 = V1'*A'*V1;
  K1 = [K1 g1; g2, g3];
  VV1=[VV1,V1];
    
    
  if (imag(znew) ~=0), cmplxflag=1;end
  z=[z,znew];

  q1 = q; 
  q2 = V2'*A*VV2; q3 = V2'*A*V2;
  %q2 = V2'*A'*VV2; q3 = V2'*A'*V2;
  K2 = [K2 q1; q2, q3];
  VV2=[VV2,V2];
    
    
end

% factored solution 

% [uY,sY]=eig(Y); [sY,id]=sort(diag(sY));
% sY=flipud(sY); uY=uY(:,id(end:-1:1));
% is=sum(abs(sY)>1e-8);
% Y0 = uY(:,1:is)*diag(sqrt(sY(1:is))); 
% Z1 = VV1(:,1:size(Y0,1))*Y0; 
% Z2 = VV2(:,1:size(Y0,1))*Y0;

[uY,sY,dY]=svd(Y); [sY,id]=sort(diag(sY));
sY=flipud(sY); uY=uY(:,id(end:-1:1)); dY=dY(id(end:-1:1),:);
is=sum(abs(sY)>1e-08);

Y0=uY(:,1:is)*diag(sY(1:is))*dY(:,1:is)';
Z=VV1(:,1:size(Y0,1))*Y0*VV2(:,1:size(Y0,1))';
final_rank=is;
sz=size(VV1,2);
t2=toc;
 fprintf(' \n')
fprintf('Space dim %d  Solution rank %d time %d\n',sz,is,t2);
 
%semilogy(nrmrestotnew,'-gd')
semilogy(err,'-bd')

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfun(x,eH,s)
 
for j=1:length(x)
r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfuni(x,eH,s)

for j=1:length(x)
r(j)=1./abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew=newpolei(eHpoints,eH,s)

for j=1:length(eHpoints)-1
%   snew(j) = fminbnd( @(x) ratfuni(x,eH,s), eHpoints(j),eHpoints(j+1),optimset('TolX',1e-3));
    sval=linspace(eHpoints(j),eHpoints(j+1),20);
   %sval=linspace(eHpoints(j),eHpoints(j+1),100);
    [sf,jx] = max (abs(ratfun(sval,eH,s)));
    snew(j)=sval(jx);
end
[sn,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function RA=compute_RA(A,eH,s);
%
%ne=length(eH);
%s=[1e20,s];
%I=speye(size(A));
%RA=I;
%for k=1:ne,
%   RA = (A-eH(k)*I)/(A-s(k)*I)*RA;
%end
%return
