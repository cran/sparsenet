c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))             
      subroutine sparsenet(no,ni,x,y,w,jd,ne,nx,ngam,nlam,fxgam,flmin,  
     *par,igrid,istart,thr,maxit,lmu,a0,ca,ia,nin,rsq,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),ca(nx,ngam,nlam),rsq(ngam,nl
     *am)
      double precision a0(ngam,nlam),par(2,ngam,nlam)                   
      integer jd(*),ia(nx),nin(ngam,nlam)                               
      double precision, dimension (:), allocatable :: xm,xs             
      integer, dimension (:), allocatable :: ju                         
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 10021                                   
      jerr=7777                                                         
      return                                                            
10021 continue                                                          
      call standard1(no,ni,x,y,w,ju,xm,xs,ym,ys,jerr)                   
      if(jerr.ne.0) return                                              
      call sparsenet2(ni,ju,y,no,ne,nx,x,ngam,nlam,fxgam,flmin,par,igrid
     *,istart,  thr,maxit,lmu,ca,ia,nin,rsq,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 10031 k=1,lmu                                                  
      do 10041 i=1,ngam                                                 
      nk=nin(i,k)                                                       
      do 10051 l=1,nk                                                   
      ca(l,i,k)=ys*ca(l,i,k)/xs(ia(l))                                  
10051 continue                                                          
      continue                                                          
      a0(i,k)=ym-dot_product(ca(1:nk,i,k),xm(ia(1:nk)))                 
10041 continue                                                          
      continue                                                          
10031 continue                                                          
      continue                                                          
      deallocate(xm,xs,ju)                                              
      return                                                            
      end                                                               
      subroutine standard1 (no,ni,x,y,w,ju,xm,xs,ym,ys,jerr)            
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),xm(ni),xs(ni)               
      integer ju(ni)                                                    
      double precision, dimension (:), allocatable :: v                 
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=w/sum(w)                                                        
      v=sqrt(w)                                                         
      do 10061 j=1,ni                                                   
      if(ju(j).eq.0)goto 10061                                          
      xm(j)=dot_product(w,x(:,j))                                       
      x(:,j)=v*(x(:,j)-xm(j))                                           
      xs(j)=sqrt(dot_product(x(:,j),x(:,j)))                            
      x(:,j)=x(:,j)/xs(j)                                               
10061 continue                                                          
      continue                                                          
      ym=dot_product(w,y)                                               
      y=v*(y-ym)                                                        
      ys=sqrt(dot_product(y,y))                                         
      y=y/ys                                                            
      deallocate(v)                                                     
      return                                                            
      end                                                               
      subroutine sparsenet2(ni,ju,y,no,ne,nx,x,ngam,nlam,fxgam,flmin,par
     *,igrid,  istart,thr,maxit,lmu,ao,ia,kin,rsqo,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, rsqmax=0.999)       
      double precision y(no),x(no,ni),ao(nx,ngam,nlam),rsqo(ngam,nlam)  
      double precision par(2,ngam,nlam)                                 
      integer ju(ni),ia(nx),kin(ngam,nlam)                              
      double precision, dimension (:), allocatable :: a,g,rsqs,a1       
      double precision, dimension (:,:), allocatable :: ys              
      integer, dimension (:), allocatable :: mm                         
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(a1(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(g(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
       allocate(ys(1:no,1:ngam),stat=ierr)                              
      jerr=jerr+ierr                                                    
      allocate(rsqs(1:ngam),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      g=0.0                                                             
      do 10071 j=1,ni                                                   
      if(ju(j).eq.0)goto 10071                                          
      g(j)=abs(dot_product(y,x(:,j)))                                   
10071 continue                                                          
      continue                                                          
      if(igrid .ne. 0)goto 10091                                        
      call pargrid(ni,g,ju,fxgam,flmin,ngam,nlam,par,alm0)              
      goto 10101                                                        
10091 continue                                                          
      alm0=maxval(g)                                                    
10101 continue                                                          
      continue                                                          
      do 10111 m=1,nlam                                                 
      do 10121 n=1,ngam                                                 
      gam=par(1,n,m)                                                    
      alm=par(2,n,m)                                                    
      if(m .ne. 1)goto 10141                                            
      if(n .ne. 1)goto 10161                                            
      a=0.0                                                             
      rsq=0.0                                                           
      call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,  rsq,nlp,nin,maxit,
     *thr,jerr)
      if(jerr.ne.0) go to 10170                                         
      ys(:,n)=y                                                         
      rsqs(n)=rsq                                                       
      goto 10181                                                        
10161 continue                                                          
      y=ys(:,n-1)                                                       
      rsq=rsqs(n-1)                                                     
      call uncomp(ni,ao(:,n-1,m),ia,kin(n-1,m),a)                       
      call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,  rsq,nlp,nin,maxit,
     *thr,jerr)
      if(jerr.ne.0) go to 10170                                         
      ys(:,n)=y                                                         
      rsqs(n)=rsq                                                       
10181 continue                                                          
      continue                                                          
      goto 10131                                                        
10141 if(n .ne. 1)goto 10191                                            
      y=ys(:,n)                                                         
      rsq=rsqs(n)                                                       
      call uncomp(ni,ao(:,n,m-1),ia,kin(n,m-1),a)                       
      call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,  rsq,nlp,nin,maxit,
     *thr,jerr)
      if(jerr.ne.0) go to 10170                                         
      ys(:,n)=y                                                         
      rsqs(n)=rsq                                                       
      goto 10201                                                        
10191 continue                                                          
      if(istart.eq.2) go to 10210                                       
      y=ys(:,n)                                                         
      rsq=rsqs(n)                                                       
      call uncomp(ni,ao(:,n,m-1),ia,kin(n,m-1),a)                       
      call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,  rsq,nlp,nin,maxit,
     *thr,jerr)
      if(jerr.ne.0) go to 10170                                         
      ys(:,n)=y                                                         
      rsqs(n)=rsq                                                       
      if(istart.eq.1) go to10220                                        
      cri=0.5*(1.0-rsq)+penalty(ni,a,gam,alm)                           
10210 continue                                                          
      y=ys(:,n-1)                                                       
      rsq=rsqs(n-1)                                                     
      call uncomp(ni,ao(:,n-1,m),ia,kin(n-1,m),a1)                      
      call soln(m,gam,alm,no,ni,x,y,a1,g,ia,mm,ju,nx,  rsq,nlp,nin,maxit
     *,thr,jerr)
      if(jerr.ne.0) go to 10170                                         
      if(istart .ne. 2)goto 10241                                       
      ys(:,n)=y                                                         
      rsqs(n)=rsq                                                       
      a=a1                                                              
      goto 10231                                                        
10241 if(0.5*(1.0-rsq)+penalty(ni,a1,gam,alm) .ge. cri)goto 10251       
      ys(:,n)=y                                                         
      rsqs(n)=rsq                                                       
      a=a1                                                              
10251 continue                                                          
10231 continue                                                          
10201 continue                                                          
10131 continue                                                          
10220 continue                                                          
      if(nin.gt.0) ao(1:nin,n,m)=a(ia(1:nin))                           
      kin(n,m)=nin                                                      
      rsqo(n,m)=rsqs(n)                                                 
      if(n.gt.1)goto 10121                                              
      me=0                                                              
      do 10261 j=1,nin                                                  
      if(ao(j,n,m).ne.0.0) me=me+1                                      
10261 continue                                                          
      continue                                                          
      if(me.lt.ne)goto 10121                                            
      if(rsqs(n).le.rsqmax)goto 10121                                   
      if(m .le. 1)goto 10281                                            
      if(rsqs(n)-rsqo(n,m-1).ge.sml*rsqo(n,m-1))goto 10121              
10281 continue                                                          
      go to 10170                                                       
10121 continue                                                          
      continue                                                          
      lmu=m                                                             
      alm0=alm                                                          
10111 continue                                                          
      continue                                                          
10170 continue                                                          
      deallocate(a,mm,g,rsqs,a1,ys)                                     
      return                                                            
      end                                                               
      subroutine pargrid(ni,g,ju,fxgam,flmin,ngam,nlam,par,alm0)        
      implicit double precision(a-h,o-z)                                
      parameter(eps=1.0e-6,big=9.9e35)                                  
      double precision g(ni),par(2,ngam,nlam)                           
      integer ju(ni)                                                    
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
      gaf=(1.0/fxgam)**(1.0/(ngam-2))                                   
      onepeps=1.0+eps                                                   
      do 10291 m=1,nlam                                                 
      if(m .le. 1)goto 10311                                            
      alm=alm*alf                                                       
      goto 10321                                                        
10311 continue                                                          
      alm0=0.0                                                          
      do 10331 j=1,ni                                                   
      if(ju(j).gt.0) alm0=max(alm0,g(j))                                
10331 continue                                                          
      continue                                                          
      alm=alf*alm0                                                      
10321 continue                                                          
      continue                                                          
      do 10341 n=1,ngam                                                 
      if(n .ne. 1)goto 10361                                            
      gam=big                                                           
      goto 10351                                                        
10361 if(n .ne. 2)goto 10371                                            
      gam=fxgam                                                         
      goto 10381                                                        
10371 continue                                                          
      gam=gaf*gam                                                       
10381 continue                                                          
10351 continue                                                          
      par(1,n,m)=max(gam,onepeps)                                       
      par(2,n,m)=alm                                                    
10341 continue                                                          
      continue                                                          
10291 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine soln (m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,rsq,nlp,nin,m
     *axit,thr,jerr)
      implicit double precision(a-h,o-z)                                
      double precision y(no),x(no,ni),a(ni),g(ni)                       
      integer ia(*),mm(*),ju(ni)                                        
      fgam=gam/(gam-1.0)                                                
      almglm=gam*alm                                                    
      continue                                                          
10391 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10401 k=1,ni                                                   
      gk=dot_product(y,x(:,k))                                          
      ak=a(k)                                                           
      u=ak+gk                                                           
      au=abs(u)                                                         
      a(k)=0.0                                                          
      if(au .le. alm)goto 10421                                         
      if(au .le. almglm)goto 10441                                      
      a(k)=u                                                            
      goto 10451                                                        
10441 continue                                                          
      a(k)=sign((au-alm)*fgam,u)                                        
10451 continue                                                          
      continue                                                          
10421 continue                                                          
      if(a(k).eq.ak)goto 10401                                          
      if(mm(k) .ne. 0)goto 10471                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 10491                                         
      jerr=-10000-m                                                     
      return                                                            
10491 continue                                                          
      mm(k)=nin                                                         
      ia(nin)=k                                                         
10471 continue                                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del)                                          
      y=y-del*x(:,k)                                                    
      dlx=max(del**2,dlx)                                               
10401 continue                                                          
      continue                                                          
      if(nlp .le. maxit)goto 10511                                      
      jerr=-m                                                           
      return                                                            
10511 continue                                                          
      if(dlx.lt.thr)goto 10392                                          
      continue                                                          
10521 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10531 l=1,nin                                                  
      k=ia(l)                                                           
      gk=dot_product(y,x(:,k))                                          
      ak=a(k)                                                           
      u=gk+ak                                                           
      au=abs(u)                                                         
      a(k)=0.0                                                          
      if(au .le. alm)goto 10551                                         
      if(au .le. almglm)goto 10571                                      
      a(k)=u                                                            
      goto 10581                                                        
10571 continue                                                          
      a(k)=sign((au-alm)*fgam,u)                                        
10581 continue                                                          
      continue                                                          
10551 continue                                                          
      if(a(k).eq.ak)goto 10531                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del)                                          
      y=y-del*x(:,k)                                                    
      dlx=max(del**2,dlx)                                               
10531 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 10522                                          
      if(nlp .le. maxit)goto 10601                                      
      jerr=-m                                                           
      return                                                            
10601 continue                                                          
      goto 10521                                                        
10522 continue                                                          
      goto 10391                                                        
10392 continue                                                          
      return                                                            
      end                                                               
      function penalty(ni,a,gam,alm)                                    
      implicit double precision(a-h,o-z)                                
      double precision a(ni)                                            
      ga=gam*alm                                                        
      tga=2.0*ga                                                        
      hga=0.5*ga                                                        
      pen=0.0                                                           
      do 10611 j=1,ni                                                   
      if(a(j).eq.0.0)goto 10611                                         
      if(abs(a(j)) .le. ga)goto 10631                                   
      pen=pen+hga                                                       
      goto 10641                                                        
10631 continue                                                          
      pen=pen+abs(a(j))-a(j)**2/tga                                     
10641 continue                                                          
      continue                                                          
10611 continue                                                          
      continue                                                          
      penalty=pen                                                       
      return                                                            
      end                                                               
      subroutine chkvars(no,ni,x,ju)                                    
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni)                                         
      integer ju(ni)                                                    
      do 10651 j=1,ni                                                   
      ju(j)=0                                                           
      t=x(1,j)                                                          
      do 10661 i=2,no                                                   
      if(x(i,j).eq.t)goto 10661                                         
      ju(j)=1                                                           
      goto 10662                                                        
10661 continue                                                          
10662 continue                                                          
10651 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine uncomp(ni,ca,ia,nin,a)                                 
      implicit double precision(a-h,o-z)                                
      double precision ca(*),a(ni)                                      
      integer ia(*)                                                     
      a=0.0                                                             
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                               
      return                                                            
      end                                                               
      subroutine modval(a0,ca,ia,nin,n,x,f)                             
      implicit double precision(a-h,o-z)                                
      double precision ca(nin),x(n,*),f(n)                              
      integer ia(nin)                                                   
      f=a0                                                              
      if(nin.le.0) return                                               
      do 10671 i=1,n                                                    
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                   
10671 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine solns(ni,nx,ngam,lmu,a,ia,nin,b)                       
      implicit double precision(a-h,o-z)                                
      double precision a(nx,ngam,lmu),b(ni,ngam,lmu)                    
      integer ia(nx),nin(ngam,lmu)                                      
      do 10681 lgam=1,ngam                                              
      do 10691 lam=1,lmu                                                
      call uncomp(ni,a(:,lgam,lam),ia,nin(lgam,lam),b(:,lgam,lam))      
10691 continue                                                          
      continue                                                          
10681 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
