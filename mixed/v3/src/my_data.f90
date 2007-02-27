
      MODULE MY_DATA
      INTEGER n,ncells,npers,nfirm,ncov,ncoef,ncoef2
      REAL*8 sig2e, sig2p, sig2f, lamp, lamf, yy

      INTEGER, ALLOCATABLE :: pers(:),firm(:),dfp(:),dff(:)

      REAL*8, ALLOCATABLE :: covin(:),y(:),d(:),cov(:,:),f(:),df(:), xx(:,:), &
        rq(:,:),xd(:,:),xf(:,:),theta(:),theta2(:),b(:),b2(:), &
        r(:),w(:),p(:),q(:),u1(:),din(:),fin(:)

      REAL*4, ALLOCATABLE :: covin4(:)

      END MODULE MY_DATA







