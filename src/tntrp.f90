!
!	Third-order polynomial interpolation using Lagrange
!     polynomials.  FUNCTION is SUBROUTINE POLINT from
!     Numerial Recipes
!
real FUNCTION tntrp(XA,YA,X,N)
!
implicit NONE
!
integer:: i,m,n,ns
real   :: den,dif,dift,dy,hO,hp,w,x,y
  real, DIMENSION(4):: XA,YA,C,D
  NS=1
  DIF=ABS(X-XA(1))
  DO 11 I=1,N
  DIFT=ABS(X-XA(I))
  IF (DIFT.LT.DIF) THEN
    NS=I
    DIF=DIFT
  ENDIF
  C(I)=YA(I)
  D(I)=YA(I)
11 CONTINUE
  Y=YA(NS)
  NS=NS-1
  DO 13 M=1,N-1
  DO 12 I=1,N-M
  HO=XA(I)-X
  HP=XA(I+M)-X
  W=C(I+1)-D(I)
  DEN=HO-HP
  IF(DEN.EQ.0.) DEN=0.001
  DEN=W/DEN
  D(I)=HP*DEN
  C(I)=HO*DEN
12 CONTINUE
  IF (2*NS.LT.N-M)THEN
    DY=C(NS+1)
  ELSE
    DY=D(NS)
    NS=NS-1
  ENDIF
  Y=Y+DY
13 CONTINUE
  tntrp=y
  RETURN
  end
