        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 19 17:12:19 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPECIMEN__genmod
          INTERFACE 
            SUBROUTINE SPECIMEN(ELEMENT,NODES,NDOFF,NELM,BC,NBC,LOADDOF,&
     &NLDOF,LOADDOF2,NLDOF2)
              INTEGER(KIND=4) :: NLDOF2
              INTEGER(KIND=4) :: NLDOF
              INTEGER(KIND=4) :: NBC
              INTEGER(KIND=4) :: NELM
              INTEGER(KIND=4) :: NDOFF
              INTEGER(KIND=4) :: ELEMENT(NELM,9)
              REAL(KIND=8) :: NODES(NDOFF,2)
              REAL(KIND=8) :: BC(NBC,2)
              INTEGER(KIND=4) :: LOADDOF(NLDOF)
              INTEGER(KIND=4) :: LOADDOF2(NLDOF2)
            END SUBROUTINE SPECIMEN
          END INTERFACE 
        END MODULE SPECIMEN__genmod
