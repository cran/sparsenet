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
      subroutine sparsenet(no,ni,x,y,w,jd,ne,nx,ngam,nlam,fxgam,flmin,      209 
     *par,igrid,istart,thr,maxit,lmu,a0,ca,ia,nin,rsq,nlp,jerr)
      real x(no,ni),y(no),w(no),ca(nx,ngam,nlam),rsq(ngam,nlam)             210
      real a0(ngam,nlam),par(2,ngam,nlam)                                   211
      integer jd(*),ia(nx),nin(ngam,nlam)                                   212
      real, dimension (:), allocatable :: xm,xs                                 
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          217
      allocate(xs(1:ni),stat=ierr)                                          217
      jerr=jerr+ierr                                                        218
      allocate(ju(1:ni),stat=ierr)                                          218
      jerr=jerr+ierr                                                        219
      if(jerr.ne.0) return                                                  220
      call chkvars(no,ni,x,ju)                                              221
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  222
      if(maxval(ju) .gt. 0)goto 10021                                       222
      jerr=7777                                                             222
      return                                                                222
10021 continue                                                              223
      call standard1(no,ni,x,y,w,ju,xm,xs,ym,ys,jerr)                       224
      if(jerr.ne.0) return                                                  225
      call sparsenet2(ni,ju,y,no,ne,nx,x,ngam,nlam,fxgam,flmin,par,igrid    227 
     *,istart,  thr,maxit,lmu,ca,ia,nin,rsq,nlp,jerr)
      if(jerr.gt.0) return                                                  228
10030 do 10031 k=1,lmu                                                      228
10040 do 10041 i=1,ngam                                                     228
      nk=nin(i,k)                                                           229
10050 do 10051 l=1,nk                                                       229
      ca(l,i,k)=ys*ca(l,i,k)/xs(ia(l))                                      229
10051 continue                                                              230
10052 continue                                                              230
      a0(i,k)=ym-dot_product(ca(1:nk,i,k),xm(ia(1:nk)))                     231
10041 continue                                                              231
10042 continue                                                              231
10031 continue                                                              232
10032 continue                                                              232
      deallocate(xm,xs,ju)                                                  233
      return                                                                234
      end                                                                   235
      subroutine standard1 (no,ni,x,y,w,ju,xm,xs,ym,ys,jerr)                236
      real x(no,ni),y(no),w(no),xm(ni),xs(ni)                               236
      integer ju(ni)                                                        237
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           240
      if(jerr.ne.0) return                                                  241
      w=w/sum(w)                                                            241
      v=sqrt(w)                                                             242
10060 do 10061 j=1,ni                                                       242
      if(ju(j).eq.0)goto 10061                                              243
      xm(j)=dot_product(w,x(:,j))                                           243
      x(:,j)=v*(x(:,j)-xm(j))                                               244
      xs(j)=sqrt(dot_product(x(:,j),x(:,j)))                                244
      x(:,j)=x(:,j)/xs(j)                                                   245
10061 continue                                                              246
10062 continue                                                              246
      ym=dot_product(w,y)                                                   246
      y=v*(y-ym)                                                            246
      ys=sqrt(dot_product(y,y))                                             246
      y=y/ys                                                                247
      deallocate(v)                                                         248
      return                                                                249
      end                                                                   250
      subroutine sparsenet2(ni,ju,y,no,ne,nx,x,ngam,nlam,fxgam,flmin,par    252 
     *,igrid,  istart,thr,maxit,lmu,ao,ia,kin,rsqo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, rsqmax=0.999)           253
      real y(no),x(no,ni),ao(nx,ngam,nlam),rsqo(ngam,nlam),par(2,ngam,nl    254 
     *am)
      integer ju(ni),ia(nx),kin(ngam,nlam)                                  255
      real, dimension (:), allocatable :: a,g,rsqs,a1                           
      real, dimension (:,:), allocatable :: ys                                  
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                           261
      allocate(a1(1:ni),stat=ierr)                                          261
      jerr=jerr+ierr                                                        262
      allocate(mm(1:ni),stat=ierr)                                          262
      jerr=jerr+ierr                                                        263
      allocate(g(1:ni),stat=ierr)                                           263
      jerr=jerr+ierr                                                        264
       allocate(ys(1:no,1:ngam),stat=ierr)                                      
      jerr=jerr+ierr                                                        268
      allocate(rsqs(1:ngam),stat=ierr)                                      268
      jerr=jerr+ierr                                                        269
      if(jerr.ne.0) return                                                  270
      mm=0                                                                  270
      nlp=0                                                                 270
      nin=nlp                                                               270
      g=0.0                                                                 271
10070 do 10071 j=1,ni                                                       271
      if(ju(j).eq.0)goto 10071                                              271
      g(j)=abs(dot_product(y,x(:,j)))                                       271
10071 continue                                                              272
10072 continue                                                              272
      if(igrid .ne. 0)goto 10091                                            272
      call pargrid(ni,g,ju,fxgam,flmin,ngam,nlam,par,alm0)                  272
      goto 10101                                                            273
10091 continue                                                              273
      alm0=maxval(g)                                                        273
10101 continue                                                              274
10081 continue                                                              274
10110 do 10111 m=1,nlam                                                     275
10120 do 10121 n=1,ngam                                                     275
      gam=par(1,n,m)                                                        275
      alm=par(2,n,m)                                                        276
      if(m .ne. 1)goto 10141                                                277
      if(n .ne. 1)goto 10161                                                277
      a=0.0                                                                 277
      rsq=0.0                                                               278
      call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,  rsq,nlp,nin,maxit,    280 
     *thr,jerr)
      if(jerr.ne.0) go to 10170                                             281
      ys(:,n)=y                                                             281
      rsqs(n)=rsq                                                           282
      goto 10181                                                            283
10161 continue                                                              283
      y=ys(:,n-1)                                                           283
      rsq=rsqs(n-1)                                                         284
      call uncomp(ni,ao(:,n-1,m),ia,kin(n-1,m),a)                           285
      call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,  rsq,nlp,nin,maxit,    287 
     *thr,jerr)
      if(jerr.ne.0) go to 10170                                             288
      ys(:,n)=y                                                             288
      rsqs(n)=rsq                                                           289
10181 continue                                                              290
10151 continue                                                              290
      goto 10131                                                            291
10141 if(n .ne. 1)goto 10191                                                291
      y=ys(:,n)                                                             291
      rsq=rsqs(n)                                                           292
      call uncomp(ni,ao(:,n,m-1),ia,kin(n,m-1),a)                           293
      call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,  rsq,nlp,nin,maxit,    295 
     *thr,jerr)
      if(jerr.ne.0) go to 10170                                             296
      ys(:,n)=y                                                             296
      rsqs(n)=rsq                                                           297
      goto 10201                                                            298
10191 continue                                                              298
      if(istart.eq.2) go to 10210                                           298
      y=ys(:,n)                                                             298
      rsq=rsqs(n)                                                           299
      call uncomp(ni,ao(:,n,m-1),ia,kin(n,m-1),a)                           300
      call soln(m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,  rsq,nlp,nin,maxit,    302 
     *thr,jerr)
      if(jerr.ne.0) go to 10170                                             303
      ys(:,n)=y                                                             303
      rsqs(n)=rsq                                                           303
      if(istart.eq.1) go to10220                                            304
      cri=0.5*(1.0-rsq)+penalty(ni,a,gam,alm)                               305
10210 continue                                                              305
      y=ys(:,n-1)                                                           305
      rsq=rsqs(n-1)                                                         306
      call uncomp(ni,ao(:,n-1,m),ia,kin(n-1,m),a1)                          307
      call soln(m,gam,alm,no,ni,x,y,a1,g,ia,mm,ju,nx,  rsq,nlp,nin,maxit    309 
     *,thr,jerr)
      if(jerr.ne.0) go to 10170                                             310
      if(istart .ne. 2)goto 10241                                           310
      ys(:,n)=y                                                             310
      rsqs(n)=rsq                                                           310
      a=a1                                                                  310
      goto 10231                                                            311
10241 if(0.5*(1.0-rsq)+penalty(ni,a1,gam,alm) .ge. cri)goto 10251           312
      ys(:,n)=y                                                             312
      rsqs(n)=rsq                                                           312
      a=a1                                                                  313
10251 continue                                                              314
10231 continue                                                              314
10201 continue                                                              315
10131 continue                                                              315
10220 continue                                                              315
      if(nin.gt.0) ao(1:nin,n,m)=a(ia(1:nin))                               315
      kin(n,m)=nin                                                          316
      rsqo(n,m)=rsqs(n)                                                     316
      if(n.gt.1)goto 10121                                                  317
      me=0                                                                  317
10260 do 10261 j=1,nin                                                      317
      if(ao(j,n,m).ne.0.0) me=me+1                                          317
10261 continue                                                              318
10262 continue                                                              318
      if(me.lt.ne)goto 10121                                                318
      if(rsqs(n).le.rsqmax)goto 10121                                       319
      if(m .le. 1)goto 10281                                                319
      if(rsqs(n)-rsqo(n,m-1).ge.sml*rsqo(n,m-1))goto 10121                  319
10281 continue                                                              320
      go to 10170                                                           321
10121 continue                                                              322
10122 continue                                                              322
      lmu=m                                                                 322
      alm0=alm                                                              323
10111 continue                                                              324
10112 continue                                                              324
10170 continue                                                              324
      deallocate(a,mm,g,rsqs,a1,ys)                                         325
      return                                                                326
      end                                                                   327
      subroutine pargrid(ni,g,ju,fxgam,flmin,ngam,nlam,par,alm0)            328
      parameter(eps=1.0e-6,big=9.9e35)                                      329
      real g(ni),par(2,ngam,nlam)                                           329
      integer ju(ni)                                                        330
      eqs=max(eps,flmin)                                                    330
      alf=eqs**(1.0/(nlam-1))                                               331
      gaf=(1.0/fxgam)**(1.0/(ngam-2))                                       331
      onepeps=1.0+eps                                                       332
10290 do 10291 m=1,nlam                                                     333
      if(m .le. 1)goto 10311                                                333
      alm=alm*alf                                                           333
      goto 10321                                                            334
10311 continue                                                              334
      alm0=0.0                                                              335
10330 do 10331 j=1,ni                                                       335
      if(ju(j).gt.0) alm0=max(alm0,g(j))                                    335
10331 continue                                                              336
10332 continue                                                              336
      alm=alf*alm0                                                          337
10321 continue                                                              338
10301 continue                                                              338
10340 do 10341 n=1,ngam                                                     339
      if(n .ne. 1)goto 10361                                                339
      gam=big                                                               339
      goto 10351                                                            340
10361 if(n .ne. 2)goto 10371                                                340
      gam=fxgam                                                             340
      goto 10381                                                            341
10371 continue                                                              341
      gam=gaf*gam                                                           341
10381 continue                                                              342
10351 continue                                                              342
      par(1,n,m)=max(gam,onepeps)                                           342
      par(2,n,m)=alm                                                        343
10341 continue                                                              344
10342 continue                                                              344
10291 continue                                                              345
10292 continue                                                              345
      return                                                                346
      end                                                                   347
      subroutine soln (m,gam,alm,no,ni,x,y,a,g,ia,mm,ju,nx,rsq,nlp,nin,m    349 
     *axit,thr,jerr)
      real y(no),x(no,ni),a(ni),g(ni)                                       350
      integer ia(*),mm(*),ju(ni)                                            351
      fgam=gam/(gam-1.0)                                                    351
      almglm=gam*alm                                                        352
10390 continue                                                              352
10391 continue                                                              352
      nlp=nlp+1                                                             353
      dlx=0.0                                                               354
10400 do 10401 k=1,ni                                                       354
      gk=dot_product(y,x(:,k))                                              355
      ak=a(k)                                                               355
      u=ak+gk                                                               355
      au=abs(u)                                                             355
      a(k)=0.0                                                              356
      if(au .le. alm)goto 10421                                             357
      if(au .le. almglm)goto 10441                                          357
      a(k)=u                                                                357
      goto 10451                                                            358
10441 continue                                                              358
      a(k)=sign((au-alm)*fgam,u)                                            358
10451 continue                                                              359
10431 continue                                                              359
10421 continue                                                              360
      if(a(k).eq.ak)goto 10401                                              361
      if(mm(k) .ne. 0)goto 10471                                            361
      nin=nin+1                                                             362
      if(nin .le. nx)goto 10491                                             362
      jerr=-10000-m                                                         362
      return                                                                362
10491 continue                                                              363
      mm(k)=nin                                                             363
      ia(nin)=k                                                             364
10471 continue                                                              365
      del=a(k)-ak                                                           365
      rsq=rsq+del*(2.0*gk-del)                                              366
      y=y-del*x(:,k)                                                        366
      dlx=max(del**2,dlx)                                                   367
10401 continue                                                              368
10402 continue                                                              368
      if(nlp .le. maxit)goto 10511                                          368
      jerr=-m                                                               368
      return                                                                368
10511 continue                                                              369
      if(dlx.lt.thr)goto 10392                                              370
10520 continue                                                              370
10521 continue                                                              370
      nlp=nlp+1                                                             370
      dlx=0.0                                                               371
10530 do 10531 l=1,nin                                                      371
      k=ia(l)                                                               371
      gk=dot_product(y,x(:,k))                                              372
      ak=a(k)                                                               372
      u=gk+ak                                                               372
      au=abs(u)                                                             372
      a(k)=0.0                                                              373
      if(au .le. alm)goto 10551                                             374
      if(au .le. almglm)goto 10571                                          374
      a(k)=u                                                                374
      goto 10581                                                            375
10571 continue                                                              375
      a(k)=sign((au-alm)*fgam,u)                                            375
10581 continue                                                              376
10561 continue                                                              376
10551 continue                                                              377
      if(a(k).eq.ak)goto 10531                                              378
      del=a(k)-ak                                                           378
      rsq=rsq+del*(2.0*gk-del)                                              379
      y=y-del*x(:,k)                                                        379
      dlx=max(del**2,dlx)                                                   380
10531 continue                                                              381
10532 continue                                                              381
      if(dlx.lt.thr)goto 10522                                              382
      if(nlp .le. maxit)goto 10601                                          382
      jerr=-m                                                               382
      return                                                                382
10601 continue                                                              383
      goto 10521                                                            384
10522 continue                                                              384
      goto 10391                                                            385
10392 continue                                                              385
      return                                                                386
      end                                                                   387
      function penalty(ni,a,gam,alm)                                        388
      real a(ni)                                                            389
      ga=gam*alm                                                            389
      tga=2.0*ga                                                            389
      hga=0.5*ga                                                            389
      pen=0.0                                                               390
10610 do 10611 j=1,ni                                                       390
      if(a(j).eq.0.0)goto 10611                                             391
      if(abs(a(j)) .le. ga)goto 10631                                       391
      pen=pen+hga                                                           391
      goto 10641                                                            392
10631 continue                                                              392
      pen=pen+abs(a(j))-a(j)**2/tga                                         392
10641 continue                                                              393
10621 continue                                                              393
10611 continue                                                              394
10612 continue                                                              394
      penalty=pen                                                           395
      return                                                                396
      end                                                                   397
      subroutine chkvars(no,ni,x,ju)                                        398
      real x(no,ni)                                                         398
      integer ju(ni)                                                        399
10650 do 10651 j=1,ni                                                       399
      ju(j)=0                                                               399
      t=x(1,j)                                                              400
10660 do 10661 i=2,no                                                       400
      if(x(i,j).eq.t)goto 10661                                             400
      ju(j)=1                                                               400
      goto 10662                                                            400
10661 continue                                                              401
10662 continue                                                              401
10651 continue                                                              402
10652 continue                                                              402
      return                                                                403
      end                                                                   404
      subroutine uncomp(ni,ca,ia,nin,a)                                     405
      real ca(*),a(ni)                                                      405
      integer ia(*)                                                         406
      a=0.0                                                                 406
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                   407
      return                                                                408
      end                                                                   409
      subroutine modval(a0,ca,ia,nin,n,x,f)                                 410
      real ca(nin),x(n,*),f(n)                                              410
      integer ia(nin)                                                       411
      f=a0                                                                  411
      if(nin.le.0) return                                                   412
10670 do 10671 i=1,n                                                        412
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                       412
10671 continue                                                              413
10672 continue                                                              413
      return                                                                414
      end                                                                   415
      subroutine solns(ni,nx,ngam,lmu,a,ia,nin,b)                           416
      real a(nx,ngam,lmu),b(ni,ngam,lmu)                                    416
      integer ia(nx),nin(ngam,lmu)                                          417
10680 do 10681 lgam=1,ngam                                                  417
10690 do 10691 lam=1,lmu                                                    418
      call uncomp(ni,a(:,lgam,lam),ia,nin(lgam,lam),b(:,lgam,lam))          419
10691 continue                                                              419
10692 continue                                                              419
10681 continue                                                              420
10682 continue                                                              420
      return                                                                421
      end                                                                   423
