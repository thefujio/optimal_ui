      subroutine nexper(n,a,mtc)
C p57 Nijenhuis & Wilf Combinatorial Alg
C Next permutation of a set 1:n
      implicit integer(a-z)
      logical mtc
      integer a(*),h
      data nlast/0/
C
C For some reason on some machines I have found a 'save' statement was needed..
C
      save
      if(n.le.0) then
         mtc=.false.
         return
      endif
 10   if(n.eq.nlast) goto 20
 30   nlast=n
      m=1
      v=1
      nf=1
      do j=1,n
         nf=nf*j
         a(j)=j
      enddo
 40   mtc=m.ne.nf
      return
 20   if(.not.mtc) goto 30
      goto(70,80)v
 70   t=a(2)
      a(2)=a(1)
      a(1)=t
      v=2
      m=m+1
      goto 40
 80   h=3
      m1=m/2
 90   b=mod(m1,h)
      if(b.ne.0) goto 120
      m1=m1/h
      h=h+1
      goto 90
 120  m1=n
      h1=h-1
      do 160 j=1,h1
        m2=a(j)-a(h)
        if(m2.lt.0) m2=m2+n
        if(m2.ge.m1) goto 160
        m1=m2
        j1=j
 160    continue
        t=a(h)
        a(h)=a(j1)
        a(j1)=t
        v=1
        m=m+1
      return
      end
C
      logical mtc
      integer a(100),fact,fn
      read(5,*) n
      k=0
      fn=fact(n)
      write(6,*)' n=',n,'  fact(n)=',fn
 1    k=k+1
      call  nexper(n,a,mtc)
C      write(6,*) (a(i),i=1,n)
      if(k.eq.1.or.k.eq.fn)write(6,*) (a(i),i=1,n)
      if(mtc) goto 1
      stop
      end
C
      integer function fact(n)
      fact=1
      if(n.le.1) return
      do i=1,n
         fact=fact*i
      enddo
      return
      end
