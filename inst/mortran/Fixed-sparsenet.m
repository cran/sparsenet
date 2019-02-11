"
c
c                          SparseNet (1/24/12)
c                                                    
c
c call sparsenet(no,ni,x,y,w,jd,ne,nx,ngam,nlam,fxgam,flmin,
c            par,igrid,istart,thr,maxit,lmu,a0,ca,ia,nin,rsq,nlp,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   x(no,ni) = predictor data matrix flat file (overwritten)
c   y(no) = response vector (overwritten)
c   w(no)= observation weights (overwritten)
c   jd(jd(1)+1) = predictor variable deletion flag
c      jd(1) = 0  => use all variables
c      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
c   ne = maximum number of variables allowed to enter largest model
c        (stopping criterion)
c   nx = maximum number of variables allowed to enter all models
c        along path (memory allocation, nx > ne).
c   ngam = number of gamma values
c   nlam = (maximum) number of lambda values
c   fxgam = largest gamma value less than infinity
c   flmin = user control of lamda values
c        minimum lamda = flmin*(largest lamda value)
c   par(2,ngam,nlam) = (gamma,lambda) values for solutions
c   igrid = 0/1 => don't/do use user supplied par values
c      igrid = 0 => default values are used and returned in the par array
c   istart = warm starting flag
c      istart=1 => lambda warm starting only
c      istart=2 => gamma warm starting only
c      istart=3 => both. Solution with smallest criterion value used.
c   thr = convergence threshold for each lambda solution.
c      iterations stop when the maximum reduction in the criterion value
c      as a result of each parameter update over a single pass
c      is less than thr times the null criterion value.
c      (suggested value, thr=1.0e-5)
c   maxit = maximum allowed number of passes over the data for all lambda
c      values (suggested values, maxit = 1000000)
c
c output:
c
c   lmu = actual number of lamda values used
c   a0(ngam,lmu) = intercept values for each (gamma,lambda) solution
c   ca(nx,ngam,lmu) = compressed coefficient values for each solution
c   ia(nx) = pointers to compressed coefficients
c   nin(ngam,lmu) = number of compressed coefficients for each solution
c   rsq(ngam,lmu) = R**2 values for each solution
c   nlp = actual number of passes over the data for all gamma-lamda values
c   jerr = error flag:
c      jerr  = 0 => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations for some gamma value.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value for some gamma value.
c
c
c
c least-squares utility routines:
c
c
c uncompress coefficient vectors for all solutions:
c
c call solns(ni,nx,ngam,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx,ngam = input to sparsenet
c    lmu,ca,ia,nin = output from sparsenet
c
c output:
c
c    b(ni,ngam,lmu) = all sparsenet returned solutions in uncompressed format
c
c
c uncompress coefficient vector for particular solution:
c
c call uncomp(ni,ca,ia,nin,a)
c
c input:
c
c    ni = total number of predictor variables
c    ca(nx) = compressed coefficient values for the solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for the solution
c
c output:
c
c    a(ni) =  uncompressed coefficient vector
c             referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call modval(a0,ca,ia,nin,n,x,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(n) = model predictions
c
c
c
"
subroutine sparsenet(no,ni,x,y,w,jd,ne,nx,ngam,nlam,fxgam,flmin,
   par,igrid,istart,thr,maxit,lmu,a0,ca,ia,nin,rsq,nlp,jerr);
implicit double precision(a-h,o-z);
double precision x(no,ni),y(no),w(no),ca(nx,ngam,nlam),rsq(ngam,nlam);
double precision a0(ngam,nlam),par(2,ngam,nlam);
integer jd(*),ia(nx),nin(ngam,nlam);
%fortran
      double precision, dimension (:), allocatable :: xm,xs
      integer, dimension (:), allocatable :: ju
%mortran
allocate(xm(1:ni),stat=jerr); if(jerr.ne.0) return;
allocate(xs(1:ni),stat=jerr); if(jerr.ne.0) return;
allocate(ju(1:ni),stat=jerr); if(jerr.ne.0) return;
call chkvars(no,ni,x,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; return;>
call standard1(no,ni,x,y,w,ju,xm,xs,ym,ys,jerr);
if(jerr.ne.0) return;
call sparsenet2(ni,ju,y,no,ne,nx,x,ngam,nlam,fxgam,flmin,par,igrid,istart,
   thr,maxit,lmu,ca,ia,nin,rsq,nlp,jerr);
if(jerr.gt.0) return;
<k=1,lmu; <i=1,ngam; nk=nin(i,k);
   <l=1,nk; ca(l,i,k)=ys*ca(l,i,k)/xs(ia(l));>
   a0(i,k)=ym-dot_product(ca(1:nk,i,k),xm(ia(1:nk)));
>>
deallocate(xm,xs,ju);
return;
end;
subroutine standard1 (no,ni,x,y,w,ju,xm,xs,ym,ys,jerr);
implicit double precision(a-h,o-z);
double precision x(no,ni),y(no),w(no),xm(ni),xs(ni); integer ju(ni);
%fortran
      double precision, dimension (:), allocatable :: v
%mortran
allocate(v(1:no),stat=jerr); if(jerr.ne.0) return;
w=w/sum(w); v=sqrt(w);
<j=1,ni; if(ju(j).eq.0) next;
   xm(j)=dot_product(w,x(:,j)); x(:,j)=v*(x(:,j)-xm(j));
   xs(j)=sqrt(dot_product(x(:,j),x(:,j))); x(:,j)=x(:,j)/xs(j);
>
ym=dot_product(w,y); y=v*(y-ym); ys=sqrt(dot_product(y,y)); y=y/ys;
deallocate(v);
return;
end;
subroutine sparsenet2(ni,ju,y,no,ne,nx,x,ngam,nlam,fxgam,flmin,par,igrid,
   istart,thr,maxit,lmu,ao,ia,kin,rsqo,nlp,jerr);
implicit double precision(a-h,o-z);
parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, rsqmax=0.999);
double precision y(no),x(no,ni),ao(nx,ngam,nlam),rsqo(ngam,nlam);
double precision par(2,ngam,nlam);
integer ju(ni),ia(nx),kin(ngam,nlam);
%fortran
      double precision, dimension (:), allocatable :: a,g,rsqs,a1
      double precision, dimension (:,:), allocatable :: ys
      integer, dimension (:), allocatable :: mm
%mortran
allocate(a(1:ni),stat=jerr); if(jerr.ne.0) return;
allocate(a1(1:ni),stat=jerr); if(jerr.ne.0) return;
allocate(mm(1:ni),stat=jerr); if(jerr.ne.0) return;
allocate(g(1:ni),stat=jerr); if(jerr.ne.0) return;
%fortran
       allocate(ys(1:no,1:ngam),stat=ierr)
%mortran
jerr=jerr+ierr;
allocate(rsqs(1:ngam),stat=jerr); if(jerr.ne.0) return;
mm=0; /nlp,nin/=0; g=0.0;
<j=1,ni; if(ju(j).eq.0) next; g(j)=abs(dot_product(y,x(:,j)));>
if igrid.eq.0 < call pargrid(ni,g,ju,fxgam,flmin,ngam,nlam,par,alm0);>
else < alm0=maxval(g);>
<m=1,nlam;
   <n=1,ngam; gam=par(1,n,m); alm=par(2,n,m);
      if m.eq.1 <
         if n.eq.1 < a=0.0; rsq=0.0;
            call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,
               rsq,nlp,nin,maxit,thr,jerr);
            if(jerr.ne.0) go to :done:;
            ys(:,n)=y; rsqs(n)=rsq;                     
         >
         else < y=ys(:,n-1); rsq=rsqs(n-1);
            call uncomp(ni,ao(:,n-1,m),ia,kin(n-1,m),a);
            call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,
               rsq,nlp,nin,maxit,thr,jerr);
            if(jerr.ne.0) go to :done:;
            ys(:,n)=y; rsqs(n)=rsq;
         >
      >
      elseif n.eq.1 < y=ys(:,n); rsq=rsqs(n);
            call uncomp(ni,ao(:,n,m-1),ia,kin(n,m-1),a);
            call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,
               rsq,nlp,nin,maxit,thr,jerr);
            if(jerr.ne.0) go to :done:;
            ys(:,n)=y; rsqs(n)=rsq;
      >
      else <  if(istart.eq.2) go to :gamma:; y=ys(:,n); rsq=rsqs(n);
         call uncomp(ni,ao(:,n,m-1),ia,kin(n,m-1),a);
         call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,
            rsq,nlp,nin,maxit,thr,jerr);
         if(jerr.ne.0) go to :done:;
         ys(:,n)=y; rsqs(n)=rsq; if(istart.eq.1) go to:out:;
         cri=0.5*(1.0-rsq)+penalty(ni,a,gam,alm);
         :gamma:y=ys(:,n-1); rsq=rsqs(n-1);
         call uncomp(ni,ao(:,n-1,m),ia,kin(n-1,m),a1);
         call soln(m,gam,alm,no,ni,x,y,a1,g,ia,mm,ju,nx,
            rsq,nlp,nin,maxit,thr,jerr);
         if(jerr.ne.0) go to :done:; 
         if istart.eq.2 < ys(:,n)=y; rsqs(n)=rsq; a=a1;>
         elseif 0.5*(1.0-rsq)+penalty(ni,a1,gam,alm).lt.cri <
            ys(:,n)=y; rsqs(n)=rsq; a=a1;
         >
      > 
      :out:if(nin.gt.0) ao(1:nin,n,m)=a(ia(1:nin)); kin(n,m)=nin;
      rsqo(n,m)=rsqs(n); if(n.gt.1) next;
      me=0; <j=1,nin; if(ao(j,n,m).ne.0.0) me=me+1;>
      if (me.lt.ne) next; if(rsqs(n).le.rsqmax) next;
      if m.gt.1 < if(rsqs(n)-rsqo(n,m-1).ge.sml*rsqo(n,m-1)) next;>
      go to :done:;
   >
   lmu=m; alm0=alm;
>
:done:deallocate(a,mm,g,rsqs,a1,ys);
return;
end;
subroutine pargrid(ni,g,ju,fxgam,flmin,ngam,nlam,par,alm0);
implicit double precision(a-h,o-z);
parameter(eps=1.0e-6,big=9.9e35);
double precision g(ni),par(2,ngam,nlam); integer ju(ni);
eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));
gaf=(1.0/fxgam)**(1.0/(ngam-2)); onepeps=1.0+eps;
<m=1,nlam;
   if m.gt.1 < alm=alm*alf;>
   else < alm0=0.0;
      <j=1,ni; if(ju(j).gt.0) alm0=max(alm0,g(j));>
      alm=alf*alm0;
   >
   <n=1,ngam;
      if n.eq.1 < gam=big;>
      elseif n.eq.2 < gam=fxgam;>
      else < gam=gaf*gam;>
      par(1,n,m)=max(gam,onepeps); par(2,n,m)=alm;
   >
>
return;
end;
subroutine soln
   (m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,rsq,nlp,nin,maxit,thr,jerr);
implicit double precision(a-h,o-z);
double precision y(no),x(no,ni),a(ni),g(ni);
integer ia(*),mm(*),ju(ni);
fgam=gam/(gam-1.0); almglm=gam*alm;
loop < 
   nlp=nlp+1; dlx=0.0;
   <k=1,ni; gk=dot_product(y,x(:,k));
      ak=a(k); u=ak+gk; au=abs(u); a(k)=0.0;
      if au.gt.alm <
         if au.gt.almglm < a(k)=u;>
         else < a(k)=sign((au-alm)*fgam,u);>
      >
      if(a(k).eq.ak) next;
      if mm(k).eq.0 < nin=nin+1;
         if nin.gt.nx < jerr=-10000-m;  return;>
         mm(k)=nin; ia(nin)=k;
      >
      del=a(k)-ak; rsq=rsq+del*(2.0*gk-del);
      y=y-del*x(:,k); dlx=max(del**2,dlx);         
   >
   if nlp.gt.maxit < jerr=-m; return;> 
   if(dlx.lt.thr) exit; 
   loop < nlp=nlp+1; dlx=0.0;
      <l=1,nin; k=ia(l); gk=dot_product(y,x(:,k));
         ak=a(k); u=gk+ak; au=abs(u); a(k)=0.0;
         if au.gt.alm <
            if au.gt.almglm < a(k)=u;>
            else < a(k)=sign((au-alm)*fgam,u);>
         >               
         if(a(k).eq.ak) next;
         del=a(k)-ak; rsq=rsq+del*(2.0*gk-del);
         y=y-del*x(:,k); dlx=max(del**2,dlx);
      >
      if(dlx.lt.thr) exit;
      if nlp.gt.maxit < jerr=-m; return;>
   >
>
return;
end;
function penalty(ni,a,gam,alm);
implicit double precision(a-h,o-z);
double precision a(ni);
ga=gam*alm; tga=2.0*ga; hga=0.5*ga; pen=0.0;
<j=1,ni; if(a(j).eq.0.0) next;
   if abs(a(j)).gt.ga < pen=pen+hga;>
   else < pen=pen+abs(a(j))-a(j)**2/tga;>
>
penalty=pen;
return;
end;
subroutine chkvars(no,ni,x,ju);
implicit double precision(a-h,o-z);
double precision x(no,ni); integer ju(ni);
<j=1,ni; ju(j)=0; t=x(1,j);
   <i=2,no; if(x(i,j).eq.t) next; ju(j)=1; exit;>
>
return;
end;
subroutine uncomp(ni,ca,ia,nin,a);
implicit double precision(a-h,o-z);
double precision ca(*),a(ni); integer ia(*);
a=0.0; if(nin.gt.0) a(ia(1:nin))=ca(1:nin);
return;
end;
subroutine modval(a0,ca,ia,nin,n,x,f);
implicit double precision(a-h,o-z);
double precision ca(nin),x(n,*),f(n); integer ia(nin);
f=a0; if(nin.le.0) return;
<i=1,n; f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)));>
return;
end;
subroutine solns(ni,nx,ngam,lmu,a,ia,nin,b);
implicit double precision(a-h,o-z);
double precision a(nx,ngam,lmu),b(ni,ngam,lmu); integer ia(nx),nin(ngam,lmu);
<lgam=1,ngam; <lam=1,lmu;
   call uncomp(ni,a(:,lgam,lam),ia,nin(lgam,lam),b(:,lgam,lam));
>>
return;
end;
%%
