
      MODULE MY_DATA
      INTEGER*4 ncells,npers,nfirm,ncov,mpoint, g
      CHARACTER(LEN=5) :: VERSION

      INTEGER*4, ALLOCATABLE :: byp(:,:),byf(:,:), m(:),mtype(:), &
                              pg(:), fg(:), pindex(:), findex(:)

      LOGICAL*1, ALLOCATABLE :: ptraced(:), ftraced(:), ponstack(:), fonstack(:)

    
      END MODULE MY_DATA







