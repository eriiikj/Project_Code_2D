module front_solver

! Last modified
! M. Ristinmaa 2021-07-03
!
!----------------------------------------------------------------------
!

   integer                    :: mnrele      ! number of elements
   integer                    :: ifloc=1     ! used in front solver
   private  mnrele, ifloc
   

!  nop     LIST OF NODES ELEMENT BY ELEMENT             
!  ien     ELEMENT LIST OBTAINED FROM FRONT PROFILER    
!           NELL SPANS OVER ALL ELEMENTS ALL GROUPS      
!           IEN(NELL,1) ELEMENT NUMBER                  
!           IEN(NELL,2) ELEMENT GROUP                   
!           IEN(NELL,3) POSITION IN NOP WHERE NODE      
!                       LIST STARTS                     
!
   integer,allocatable        :: nop(:)
   integer,allocatable        :: ien(:,:)
   private nop, ien
   
   integer        :: knuse
   PARAMETER(KNUSE=   450 ) ! size of front should set in optielem
   private  knuse

!
! enod when different elements are used
!
     type enodtype
       integer, allocatable  :: enod(:,:)
     end type
!
!     parameter(maxdom=10)
!     type(enodtype)    :: enod_dom(maxdom)


  interface optielem
     module procedure optielem_single
     module procedure optielem_multi
  end interface
  private optielem_single, optielem_multi


!----------------------------------------------------------------------
contains


! -- Front optimizing part ----


! NAME   :  OPTIELEM                                              C
!                                                                 C
! PURPOSE:  MINIMIZES FRONT USED IN DRSOLV                        C
!                                                                 C
! CALLS  :  SETUP, DIAM, RESEQ1, RESEQ2                           C
!                                                                 C
!           LEGRU     DEFINING ELEMENT GROUP IGR                  C
!                      LEGRU(IGR,1) FIRST ELEMENT NUMBER          C
!                      LEGRU(IGR,2) LAST ELEMENT NUMBER           C
!                      LEGRU(IGR,3) NUMBER OF NODES IN AN ELEMENT C
!                      LEGRU(IGR,4) NUMBER OF INTEGRATION POINTS  C
!                      LEGRU(IGR,5) NUMBER OF STATE VARIABLES/    C
!                                   INT. POINT                    C
!                      LEGRU(IGR,6) INTEGRATION TYPE              C
!                      LEGRU(IGR,7) STRESS TYPE                   C
!                      LEGRU(IGR,8) ELEMENT TYPE                  C
!                      LEGRU(IGR,10) PLACTICITY MODEL             C
!                      LEGRU(IGR,11) NR. DOF/NODE                 C
!                                                                 
! legru should contain 
!   - element type
!   - material model


subroutine optielem_single(enod)

   implicit none

   integer                     :: enod(:,:)
   type(enodtype)              :: enod_dom(1)
   
   integer                     :: legru(1,11)
   integer                     :: NDELE, KNELE
   
   allocate(enod_dom(1)%enod(size(enod,1),size(enod,2)))
   enod_dom(1)%enod=enod
   
   NDELE=size(enod,1)  ! number of nodes in element
   KNELE=size(enod,2)  ! number of elements

   LEGRU(1,1)=1
   LEGRU(1,2)=KNELE
   LEGRU(1,3)=NDELE
   LEGRU(1,4)=0 ! not used
   LEGRU(1,5)=0 ! not used
   LEGRU(1,6)=0 ! not used
   LEGRU(1,7)=0 ! not used
   LEGRU(1,8)=0 ! not used here
   LEGRU(1,9)=0 ! not used
   LEGRU(1,10)=0 ! not used
   LEGRU(1,11)=2 ! not used
  
   call optielem_multi(enod_dom, legru)

end subroutine optielem_single

subroutine optielem_multi(enod_dom, legru)

   implicit none

   integer                     :: legru(:,:)
   type(enodtype)              :: enod_dom(:)
   
   
   integer                     :: NDELE, KNELE, KNNOD, KELEN ,KNGRU
   integer                     :: NODES, NET, NVN, MINMAX
   integer                     :: I, NS, J, II
   
   integer                     :: MNGRU
   
   integer, allocatable        :: LELND(:), NEN(:), NDEG(:), LEV(:)
   integer, allocatable        :: NEWNN(:), NEWNUM(:), NSTART(:)
   integer, allocatable        :: NADJ(:)

   integer                     :: maxdeg=40

!
! notera att om har olika antal noder per element så fungerar
! så fungerar en matris (endo) inte lika bra en vektor (LELND)
! kanske använda en pekar struktur?

! currently only for one domain
      if (size(enod_dom).gt.1) stop


      NDELE=size(enod_dom(1)%enod,1)  ! number of nodes in element
      KNELE=size(enod_dom(1)%enod,2)  ! number of elements
      KNNOD=maxval(enod_dom(1)%enod)  ! max nod value number
      KELEN=NDELE*KNELE   ! number of entries in enod

!
! stored in module use in front solver
!
      mnrele=knele      
      
!
! Only one mesh domain
      KNGRU=size(legru,1)
  
      allocate(ien(KNELE,4))   ! 4: how many identifiers do we need
      allocate(nop(KELEN))
 
! Create a vector nop with nodes element by element
      ii=0
      do j=1,knele
        do i=1,NDELE
          ii=ii+1
          nop(ii)=enod_dom(1)%enod(i,j)
        end do
      end do
      
      allocate(NEN(KNELE))
      allocate(LEV(KNNOD))
      allocate(NDEG(KNNOD))
      allocate(NEWNN(KNNOD))
      allocate(NEWNUM(KNNOD))
      allocate(NSTART(KNNOD))
      allocate(NADJ(KNNOD*MAXDEG))

!
! CONTROL
!
!      KNPOS=MNPOS
!      WRITE(*,*)' CONTROL ',KELEN,KNGRU,KNPOS,KNELE,KNNOD,maxdeg
      NODES=KNNOD
      NET=KNELE
      NVN=1
      MINMAX=10**10
      
! maybe not needed
      MNGRU=1     

      CALL SETUP(nop,NADJ,NDEG,LEGRU,NODES,NET,MAXDEG)

      CALL DIAM(NDEG, NSTART, LEV, NADJ, NODES, MAXDEG, NS)

	WRITE(*,*)'LIST OF POSSIBLE STARTING NODES'
	DO 20 I=1,NS
 20          WRITE(*,*)'NODE--->',NSTART(I)

      CALL RESEQ1(NADJ, NDEG, NEWNN, NEWNUM, NSTART, NODES, &
                       MAXDEG, NS, MINMAX)

!
! minmax should be used to set fron size for front solver
! knuse
!                      
	WRITE(*,*)'AFTER OPTIMIZE----------'
	WRITE(*,*)'NEW FRONT SIZE IN NODES--->',MINMAX

      CALL RESEQ2(NEWNUM, nop, NEN, IEN, LEGRU, NET, NODES, &
                       NVN)

      

      deallocate(NEN)
      deallocate(LEV)
      deallocate(NDEG)
      deallocate(NEWNN)
      deallocate(NEWNUM)
      deallocate(NSTART)
      deallocate(NADJ)

!      STOP
      RETURN
end subroutine optielem_multi
      
      
SUBROUTINE SETUP(NPN, NADJ, NDEG, LEGRU, NODES, NET, MAXDEG)
!*********************************************************************
!     SUBPROGRAM SETUP - COMPUTE ADJACENCY LIST AND DEGREE FOR
!                        EACH NODE 
!                      - MODIFIED VERSION OF COLLINS ROUTINE
!                      - USE ONLY CORNER NODES IF ALL = .FALSE.
!*********************************************************************
      implicit none
      
      integer                      :: NPN(:), NADJ(:), NDEG(:), LEGRU(:,:)
      integer                      :: NODES, NET, MAXDEG
      
      integer                      :: KNGRU, IGR, NRELE
      integer                      :: IEF, IEL, NN, J, I, JNTI, JSUB, II
      integer                      :: JJT, MEM1, III
      
      NDEG=0

!
! LOOP OVER ALL ELEMENT GROUPS
!
      KNGRU=size(LEGRU,1)
      NRELE=0
      DO 70 IGR=1,KNGRU
        IEF=LEGRU(IGR,1)
        IEL=LEGRU(IGR,2)
        NN =LEGRU(IGR,3)
        
      DO 60 J=IEF, IEL
 
!  write(*,*)'setup ele ',j
  
        DO 50 I=1,NN
          JNTI = NPN(NRELE+I)
          JSUB = (JNTI-1)*MAXDEG

          DO 40 II=1, NN
            IF (II.EQ.I) GOTO 40
            JJT  = NPN(NRELE+II)              
            MEM1 = NDEG(JNTI)
            IF (MEM1.EQ.0) GOTO 30

            DO 20 III=1, MEM1
              IF (NADJ(JSUB+III).EQ.JJT) GOTO 40
   20       CONTINUE

   30       NDEG(JNTI) = NDEG(JNTI)+1
            NADJ(JSUB+NDEG(JNTI)) = JJT
            IF (NDEG(JNTI).GT.MAXDEG) THEN
               WRITE(*,*)'MODULE FRONT_SOLVER'
               WRITE(*,*)'OPTIFRONT UNABLE TO SOLVE GRAPH PROBLEM ', &
                       'MAXDEG MUST BE INCREASED TO',NADJ(JNTI), NDEG(JNTI)
               STOP
            ENDIF
   40     CONTINUE

   50   CONTINUE

        NRELE=NRELE+NN
   60 CONTINUE

   70 CONTINUE

      RETURN
END SUBROUTINE SETUP


      SUBROUTINE DIAM(NDEG, NSTART, LEV, NADJ, NODES, MAXDEG, NS)
!****************************************************************
!     SUBPROGRAM DIAM - COMPUTE SET OF PSUEDO-PERIPHERAL NODES
!****************************************************************
      DIMENSION NDEG(1), NADJ(1), LEV(1), NSTART(1)
      LOGICAL BETTER
!
!     BEGIN ITERATION
!     SELECT INITIAL ROOT NODE ARBITRARILY AND GENERATE
!     ITS LEVEL STRUCTURE

      IROOT = 1
      ITER  = 0
   10 ITER  = ITER+1
      CALL LEVEL(NDEG, LEV, IDEPTH, NADJ, IWIDTH, NODES, IROOT, &
                MAXDEG)
!
!     CREATE LIST OF NODES WHICH ARE AT MAXIMUM DISTANCE
!     FROM ROOT NODE

      LHW = 0
      DO 20 I=1, NODES
        IF (LEV(I).NE.IDEPTH) GOTO 20
        LHW = LHW+1
        NSTART(LHW) = I
   20 CONTINUE
!
!     STORE ROOT ON END OF LIST OF POSSIBLE STARTING NODES
!
      NS = LHW+1
      NSTART(NS) = IROOT

!     LOOP OVER NODES AT MAXIMUM DISTANCE FROM ROOT NODE
!     GENERATE LEVEL STRUCTURE FOR EACH NODE
!     SET SWITCH IF A LEVEL STRUCTURE OF GREATER DEPTH OCCURS

      BETTER = .FALSE.
      DO 30 I=1, LHW
        NEND = NSTART(I)
        CALL LEVEL(NDEG, LEV, NDEPTH, NADJ, NWIDTH, NODES, NEND, &
                  MAXDEG)
        IF (NDEPTH.LT.IDEPTH) GOTO 30
        IF ((NDEPTH.EQ.IDEPTH).AND.(NWIDTH.GE.IWIDTH)) GOTO 30
        IROOT  = NEND
        IDEPTH = NDEPTH
        IWIDTH = NWIDTH
        BETTER = .TRUE.
   30 CONTINUE
      IF (BETTER) GOTO 10

      RETURN
      END


      SUBROUTINE LEVEL(NDEG, LEV, LSD, NADJ, MLW, NODES, NROOT,  &
                      MAXDEG)
!*********************************************************************
!     SUBPROGRAM LEVEL - COMPUTE LEVEL STRUCTURE ROOTED AT NROOT
!*********************************************************************
      DIMENSION LEV(1), NDEG(1), NADJ(1)
!
!     INITIALISATION
!
      DO 10 I=1, NODES
   10 LEV(I) = 0
      LEV(NROOT) = 1
      KOUNT = 1
      MLW   = 1
!
!     ASSIGN LEVELS TO VERTICES
!
      DO 40 L=2, NODES
        LW = 0

        DO 30 I=1, NODES
          IF (LEV(I).GT.0) GOTO 30
          NCS  = NDEG(I)
          JSUB = (I-1)*MAXDEG

          DO 20 JJ=1, NCS
            NODE = NADJ(JSUB+JJ)
            IF (LEV(NODE).NE.L-1) GOTO 20
            LSD = L
            LW  = LW+1
            LEV(I) = L
            KOUNT  = KOUNT+1
            IF (KOUNT.EQ.NODES) GOTO 50
            GOTO 30
   20     CONTINUE

   30   CONTINUE
        IF (LW.GT.MLW) MLW = LW

   40 CONTINUE
   50 IF (LW.GT.MLW) MLW = LW

      RETURN

      END


      SUBROUTINE RESEQ1(NADJ, NDEG, NEWNN, NEWNUM, NSTART, NODES, &
                       MAXDEG, NS, MINMAX)
!********************************************************************
!     SUBPROGRAM - RESEQ1 - RESEQUNCE NODES FOR MINIMUM FRONTWIDTH
!********************************************************************
      DIMENSION NADJ(1), NDEG(1), NEWNN(1), NEWNUM(1), NSTART(1)
!
!     LOOP OVER SET OF STARTING NODES
!

      DO 100 II=1, NS
        I = NSTART(II)
        DO 10 J=1, NODES
          NEWNN(J) = 0
   10   CONTINUE          
        NIF      = NDEG(I)
        MAXFRT   = NIF
        NEWNN(I) = 1
!
!       NEGATE ALL NDEG ENTRIES FOR NODES WHICH
!       ARE ADJACENT TO STARTING NODE I
!
        NCN  = NDEG(I)
        JSUB = (I-1)*MAXDEG
        DO 20 J=1, NCN
          N = NADJ(JSUB+J)
          NDEG(N) = -NDEG(N)
   20   CONTINUE
        NDEG(I) = -NDEG(I)
!
!       LOOP OVER NODES TO BE RENUMBERED
!
        DO 60 K=2, NODES
          MINNEW = 10**10
          LMIN   = 10**10
!
!         LOOP OVER UNNUMBERED NODES
!         SKIP TO NEXT NODE IF OLD NODE IS ALREADY RENUMBERED
!         RESTRICT SEARCH TO ACTIVE NODES FOR KING SCHEME
!
          DO 40 J=1, NODES
            IF ((NEWNN(J).GT.0).OR.(NDEG(J).GT.0)) GOTO 40
            NEW  = 0
            MIN  = 10**10
            NCN  = IABS(NDEG(J))
            LSUB = (J-1)*MAXDEG
!
!           COMPUTE THE INCREMENT IN ACTIVE NODES FOR EACH NODE J
!           COMPUTE WHEN NODE WAS FIRST ACTIVATED BY CHECKING FOR
!           RENUMBERED NEIGHBOURS WITH LOWEST NUMBERS
!
            DO 30 L=1, NCN
              N = NADJ(LSUB+L)
              IF (NDEG(N).GT.0) NEW = NEW+1
              IF (NEWNN(N).EQ.0) GOTO 30
              IF (NEWNN(N).LT.MIN) MIN = NEWNN(N)
   30       CONTINUE
!
!           SELECT NODE WITH SMALLEST INCREMENT IN ACTIVE NODES
!           IN THE CASE OF A TIE, SELECT NODE WHICH HAS BEEN 
!           ACTIVE THE LONGEST
!
            IF (NDEG(J).LT.0) NEW = NEW-1
            IF (NEW.GT.MINNEW) GOTO 40
            IF ((NEW.EQ.MINNEW).AND.(MIN.GE.LMIN)) GOTO 40
            MINNEW = NEW
            LMIN   = MIN
            NEXT   = J
   40     CONTINUE
!
!         RENUMBER NODE AND COMPUTE NUMBER OF ACTIVE NODES  
!         ABANDON SCHEME IF NUMBER OF ACTIVE NODES EXCEEDS
!         PREVIOUS LOWEST MAXIMUM FRONTWIDTH
!
          NEWNN(NEXT) = K
          NIF = NIF+MINNEW
          IF (NIF.GT.MAXFRT) MAXFRT = NIF
          IF (MAXFRT.GE.MINMAX) GOTO 80
!
!         NEGATE ALL NDEG ENTRIES FOR NODES WHICH ARE
!         ADJACENT TO NODE JUST RENUMBERED
!
          IF (MINNEW.EQ.-1) GOTO 60
          NCN  = IABS(NDEG(NEXT))
          JSUB = (NEXT-1)*MAXDEG
          DO 50 J=1, NCN
            N = NADJ(JSUB+J)
            IF (NDEG(N).GT.0) NDEG(N) = -NDEG(N)
   50     CONTINUE

   60   CONTINUE
!
!       STORE NUMBERING SCHEME GENERATED
!       RESET NDEG TO POSITIVE VALUES
!
        DO 70 J=1, NODES
   70     NEWNUM(J) = NEWNN(J)
        MINMAX = MAXFRT
   80   DO 90 J=1, NODES
   90     NDEG(J) = IABS(NDEG(J))
  100 CONTINUE
      MINMAX = MINMAX+1

      RETURN
      END


      SUBROUTINE RESEQ2(NEWNUM, NPN, NEN, IEN, LEGRU, NET, &
                       NODES,NVN)
!**********************************************************************
!     SUBPROGRAM RESEQ2 - RESEQUENCE ELEMENT NUMBERS TO MINIMISE
!                         THE FRONTWIDTH  
!                       - REORDER THE ELEMENTS IN AN ASCENDING SEQUENCE
!                         OF THEIR LOWEST NUMBERED NODES
!**********************************************************************

      integer   :: newnum(:), npn(:), nen(:), ien(:,:), legru(:,:)
!      DIMENSION NEWNUM(1), NPN(1), NEN(1), IEN(MNELE,4), LEGRU(MNGRU,KNPOS)

      integer  :: kngru
      
      kngru=size(legru,1)

      DO 10 I=1, NET
        NEN(I) = 0
   10 CONTINUE        
      KOUNT = 0
!         
!     LOOP OVER EACH NODE NUMBER
!     LOOP ONLY OVER CORNER NODES IF ALL = .FALSE.
!
      DO 40 I=1, NODES

      NRELE=0
      INUMG=0
      DO 60 IGR=1,KNGRU
         IEF=LEGRU(IGR,1)
         IEL=LEGRU(IGR,2)
         NN =LEGRU(IGR,3)
         NRI=LEGRU(IGR,4)
!         
!       LOOP OVER EACH ELEMENT
!       SKIP TO NEXT ELEMENT IF ALREADY RENUMBERED
!
        DO 35 J=IEF, IEL
          IF (NEN(J).GT.0) GOTO 30
          I1 = NRELE
!
!         LOOP OVER EACH NODE IN ELEMENT
!
          DO 20 K=1, NN
            N = NPN(I1+K)
            N = NEWNUM(N)
            IF (N.NE.I) GOTO 20
            KOUNT      = KOUNT+1
            NEN(J)     = KOUNT
            IEN(KOUNT,1) = J
            IEN(KOUNT,2) = IGR
            IEN(KOUNT,3) = NRELE
            IEN(KOUNT,4) = INUMG
            IF (KOUNT.EQ.NET) THEN
              GOTO 50
            ENDIF
            GOTO 30
   20     CONTINUE
   30     NRELE=NRELE+NN
   35   CONTINUE
   60  CONTINUE
   40 CONTINUE

   50 RETURN
      END 

!---- Front solver part 


! NAME   :  DRSOLV                                                C
!                                       -1                        C
! PURPOSE:  Calculate  KX = F  =>  X = K  F                       C
!                                                                 C
! CALLS  :   BACKSUB                                              C
!                                                                 C
! INPUT  :   (*) INDICATES THAT THIS ENTRY IS NOT USED BY DRSOLV, C
!                BUT BY ELSTIF                                    C
!                                                                 C
!            RESID    THE FORCE VECTOR F                          C
!                                                                 C
!            NCOD     CONTROL OF BOUNDARY CONDITIONS              C
!                       1: X PRESCRIBED                           C
!                       0: F PRESCRIBED                           C
!                                                                 C
!            BDISP    PRESCIBED VALUES IN X                       C
!
! these both from enod                                                                 C
!            NE       NUMBER OF ELEMENTS                    C
!            NH       NUMBER OF NODES                       C
!                                                                 C
!            LEGRU    DEFINING ELEMENT GROUP IGR                  C
!                      LEGRU(IGR,1) FIRST ELEMENT NUMBER          C
!                      LEGRU(IGR,2) LAST ELEMENT NUMBER           C
!                      LEGRU(IGR,3) NUMBER OF NODES IN AN ELEMENT C
!                      LEGRU(IGR,4) NUMBER OF INTEGRATION POINTS  C
!                      LEGRU(IGR,5) NUMBER OF STATE VARIABLES/    C
!                                   INT. POINT                    C
!                      LEGRU(IGR,6) INTEGRATION TYPE              C
!                      LEGRU(IGR,7) STRESS TYPE                   C
!                      LEGRU(IGR,8) ELEMENT TYPE                  C
!                      LEGRU(IGR,10) PLACTICITY MODEL             C
!                      LEGRU(IGR,11) NR. DOF/NODE                 C
!                                                                 C
!                                                                 C
!            IERRO    ERROR INDICATOR 0 IF NO ERROR 1 IF ERROR    C
!                                                                 C
!            IFLOC    FORCE SETUP. THIS MODIFIES NOP. SHOULD BE   C
!                     1 FIRST CALL                                C
!                                                                 C
!     (*)    COORD    LIST OF NODE COORDINATES                    C
!                        COORD(DIM*(NODE-1)+1) = X                C
!                        COORD(DIM*(NODE-1)+2) = Y                C
!                        COORD(DIM*(NODE-1)+3) = Z  (DIM=3)       C
!                                                                 C
!     (*)    DISPI    ITERATIVE (TOTAL) VALUES OF DISPLACEMENTS   C
!                                                                 C
!                                                                 C
!                                                                 C
! OUTPUT :   SOLVE      THE SOLUTION VECTOR X                     C
!                                                                 C

! input enod???

!  subroutine drsolv(solve,resid,ncod,bdisp,coord,disp,enod,dofnod)

! legru should be stored in element module such that only element 
! is needed in the call for the element stiffness
                        
      SUBROUTINE DRSOLV(SOLVE,RESID,NCOD ,BDISP , &
                        LEGRU,IERRO, COORD, &
                        DISPI,enod,dofnod)

      IMPLICIT REAL*8(A-H,O-Z)

      double precision              :: solve(:), resid(:), bdisp(:)
      double precision              :: coord(:)  ! should be changed
      double precision              :: dispi(:)

      integer                       :: ncod(:), legru(:,:), enod(:,:), dofnod

      integer                       :: NE    ! number of elements
      integer                       :: NH    ! number of nodes
      
! must be changed
! note kndim=2 means 2-dimensional calculation
!

      PARAMETER( &
       KNDIM=     2, &   ! dofs per nod must be changed
       KNDEL=     8)     ! maximum size of stiffness matrix

! note if many domains, size of stiff etc could differ
! should be determined in an initiation routine
     
      DIMENSION &
        STIFF(KNDEL,KNDEL), &
        KDEST(KNDEL), &
        NK   (KNDEL)
!
! kndof total dof        
!     1  KMEMO(KNDOF,KNuse),
!     2  QQMEM(KNDOF,KNuse),

!
! front solver matrix size should be move do module head
! such that size can be allocated
!
      double precision           :: equat(knuse,knuse), qq(knuse), pvkol(knuse)
      integer                    :: khed(knuse), kpiv(knuse), jmod(knuse)
       
!
       KNGRU=size(legru,1)
       KNELE=mnele
       KELND=size(nop)
       
       NE=size(enod,2)
       NH=maxval(enod)
!
!      WRITE(*,*)' SUB  -> DRSOLV'
!
! initial values
!  NPIVO: 1 No pivo search
!         0 Pivo search active
!  ND1  : 2 Disk number two store equations
!  NH   :   Number of nodes 
!             
      NPIVO=1
      ND1=2
      NMAX=KNUSE
!
      IF (NPIVO.EQ.0) THEN
        NCRIT=NMAX-9
      ELSE
        NCRIT=0
      ENDIF
!
! Find last occurence of node in NOP
! NOTE! This is only done one time, if mesh is intact
! ifloc global in module
!                                             
      IF (IFLOC.EQ.1) THEN
        IFLOC=0
!
        DO 150 J=1,KELND
 150      NOP(J)=IABS(NOP(J))
!
!        I=0
! 12     CONTINUE
!          I=I+1
!          L=KELND+1
! 8        CONTINUE
!            L=L-1
!            IF (NOP(L).EQ.I) THEN
!              NOP(L)=-I
!              GOTO 13
!            ENDIF
!          IF (L.GT.1) GOTO 8
!          WRITE(*,*)'SOMETHING IS WRONG'
! 13     IF (I.LT.NH) GOTO 12
!        IFLOC=0
!      ENDIF
!
      DO 12 IN=1,NH
        DO 8 IE=KNELE,1,-1
          IGR=IEN(IE,2)
          NRN=IEN(IE,3)
          NN =LEGRU(IGR,3)
          DO 13 INN=NN,1,-1
            IF(IN.EQ.NOP(NRN+INN)) THEN
              NOP(NRN+INN)=-IN
              GOTO 12
            ENDIF
  13      CONTINUE
  8     CONTINUE
        WRITE(*,*)'SOMETHING IS WRONG WITH THE NODE NUMBER',IN
        STOP
  12   CONTINUE
      ENDIF
!
!
! zero front
!
       KROW = 0
       DO 16  I=1 , NMAX
          DO 16  J=1 , NMAX
             EQUAT(J,I) = 0.
   16  CONTINUE
!
       WRITE(ND1) 1
       REWIND ND1    
!
      IFAIL=0
      INUMG=0
      NELND=0
      IGR=0
      IME=0
!
! LOOP OVER ALL ELEMENTS
!    
      NE=KNELE
      NELL=0
 18   CONTINUE
        NELL=NELL+1
        NELO =IEN(NELL,1)
        IGR  =IEN(NELL,2)
        NELND=IEN(NELL,3)
        INUMG=IEN(NELL,4)
! 	WRITE(*,*)'ELEMENT--->',NELL,NELO
!        WRITE(*,*)' ELEMENT GRUPP-->',IGR
        IEFIR=LEGRU(IGR,1)
        IELAS=LEGRU(IGR,2)
        NODEL=LEGRU(IGR,3)
        NRINT=LEGRU(IGR,4)
        NRSTA=LEGRU(IGR,5)
        INTYP=LEGRU(IGR,6)
        LSTYP=LEGRU(IGR,7)
        LETYP=LEGRU(IGR,8)
        LPTYP=LEGRU(IGR,10)
        NDNOD=LEGRU(IGR,11)
!        WRITE(*,*)' MATERIAL TYPE LSTYPE-->',LSTYP
!        WRITE(*,*)' GAUSS TYPE --->',INTYP
!
!  NBN: NR NODES IN ELEMENT
!  NCN: NR DOF IN ELEMENT NOTE!!!!! MUST BE CHANGED
!
        NBN=NODEL
        NCN=NBN*NDNOD
!
!
! calculate element-stiffnes matrix
!
!       WRITE(*,*)'ELEMENT',NELL
!
!          IF (NELL.LT.10) WRITE(*,*)'RESID-INUMG-->',INUMG,IGR
!!       CALL ELSTIF(STIFF,COORD,NOP  ,VGASP,DISPN,DISPI,RESID,       &
!!                   NODEL,NRINT,INTYP,NELO ,IFILE,KFILE,IERRO,LPTYP, &
!!                   PROPS,NELND,INUMG,NRSTA,LBTYP, &
!!                   ITERA,YMATE,IGR,LKONT,LSTYP)
!
!        call dneohooke('tl',D,dg(:,ie))     
!        call c2dtl3_e(ke,coord(:,enod(:,ie)),ep,D,ed(:,ie),es([1,2,4],ie))
!
!        call dneohooke('tl',D,ie)     
!        call c2dtl3_e(ke,coord,enod,ed,ie)
!
! generic call, everything with which element and material used
! should have done in an initiation phase. Only element number is 
! needed then everything is known, i.e. element numer is connected
! to element and material
!
!
!       call elstif(ke,coord,enod,ed,ie)
!

       IF (IERRO.EQ.1) THEN 
          WRITE(*,*)' ERROR ELEMENT NR ',NELO,NELL
          WRITE(*,*)' NRELD',NELND,NODEL,INUMG,NRSTA
          GOTO 999
        ENDIF
!       DO 33 IDD1=1,8
! 33      WRITE(*,*)'STIF-->',(STIFF(IDD1,IDD2),IDD2=1,8)
 
!       WRITE(*,*)'ASSEMBLING ELEMENT',NELL
!
!  THIS SHOULD BE INTRODUCED
!    MDF(M): NR.DOFS FOR NODE M
!
         KC   = 0
         DO 22  J=1, NBN
           NN  = NOP(J+NELND)
           M   = IABS(NN)
!          K   = NOPP(M)
           K   = KNDIM*(M-1)+1
!          IDF = MDF(M)
!           IDF = KNDIM
           IDF = NDNOD
           DO 23  L=1 , IDF    
             KC = KC + 1
             II = K + L - 1
             IF (NN.LT.0) II = - II
             NK(KC) = II
   23      CONTINUE
   22    CONTINUE
!
!         NELND=NELND+NBN
!
! set up heading vectors
!
!       WRITE(*,*)' SETUP HEADING VECTOR'
       DO 52  LK=1 , NCN
         KDOF = NK(LK)
!
 36      IF (KROW.GT.0) THEN
          DO 42  K=1 , KROW
             KK = K
             IF (IABS(KDOF).EQ.IABS(KHED(K))) GOTO 48
   42     CONTINUE
         ENDIF  
         KROW       = KROW + 1
         KDEST(LK)  = KROW
         KHED(KROW) = KDOF
         GOTO 52
 48      KDEST(LK) = KK
         KHED(KK)  = KDOF
   52  CONTINUE
!
       IF (KROW.GT.NMAX) THEN
          NERROR = 2 
          WRITE(*,*) 'NERROR = ',NERROR, &
          ' THE DIFFERENCE NMAX-NCRIT IS NOT SUFFICINTLY LARGE', &
          ' TO PERHIT THE ASSEMBLY OF THE NEXT ELEMENT ---    ', &
          ' EITHER INCREASE NMAX OR LOWER NCRIT'
  417     FORMAT(A,I5,3A)
	  WRITE(*,*)NCRIT,NMAX
	  DO  I=1,KROW
            WRITE(*,*)'KROW-',I,KDEST(I),KHED(I)
          END DO	  
          STOP 
       ENDIF 
!
! assembly
!
!       WRITE(*,*)' --->ASSEMBLING'
       DO 57  L=1 , NCN   
          LL = KDEST(L)
          DO 56  K=1 , NCN
             KK        = KDEST(K) 
             EQUAT(KK,LL) = EQUAT(KK,LL) + STIFF(K,L)
   56     CONTINUE
   57  CONTINUE
!
!      WRITE(*,*)' AFTER ASSEMBLING'
! 17   IF (KROW.LT.NCRIT.AND.NELL.LT.KNELE.AND.NPIVO.EQ.0) GOTO 18
!       WRITE(*,*)' AFTER TEST ASSE'
!
! find out which matrix element are fully summed
!
 60    IR = 0
       KR = 0
       DO 68  K=1 , KROW
!          WRITE(*,*)'TEST-->',K
          KT = KHED(K)
          IF (KT.LT.0) THEN
             KR       = KR + 1
             KPIV(KR) = K
             KRO      = IABS(KT)
!             WRITE(*,*)' TEST2-->',KR,KRO
             IF (NCOD(KRO).EQ.1) THEN
!                WRITE(*,*)' INDEX--->',IR,KRO
                IR        = IR + 1
                JMOD(IR)  = K
                NCOD(KRO) = 2
             ENDIF      
          ENDIF   
   68  CONTINUE
!
! modify equations with applied boundary conditions
!
!       WRITE(*,*)' MODIFY FOR BOUNDARY COND'
       IF (IR.GT.0) THEN
         DO 70  IRR=1 , IR
           K  = JMOD(IRR)
           KH = IABS(KHED(K))
           DO 69  J=1 , KROW
             LH     =IABS(KHED(J))
             RESID(LH) =RESID(LH)-EQUAT(J,K)*BDISP(KH)
             EQUAT(J,K)=0.
 69        CONTINUE
           EQUAT(K,K)=-1.
 70      CONTINUE
       ENDIF
!                                                    
!       WRITE(*,*)' BEFORE IF'
       IF (KR.EQ.0.AND.NELL.LT.KNELE.AND.NPIVO.EQ.1) GOTO 18
!       WRITE(*,*)' EMELLAN'
!       IF (KR.EQ.0.AND.NPIVO.EQ.1.AND.IGR.LT.KNGRU) GOTO 9
!       WRITE(*,*)' AFTER IF'
!
       IF (KR.EQ.0) THEN
          NERROR = 3
          WRITE(*,*) 'NERROR = ',NERROR, &
          ' THERE ARE NO MORE ROWS FULLY SUMMED, THIS MAY BE',  &
          ' DUE TO ---- ',  &
          ' (1) INCORRECT CODING OF NOP OR NK ARRAYS ',  &
          ' (2) INCORRECT VALUE OF NCRIT. INCREASING NCRIT TO',  &
          '     TO PERMIT WHOLE FRONT TO BE ASSEMBLED'
  418     FORMAT(A,I5,5A)
          STOP  
       ENDIF                                          
!
! search for absolute pivot
!
       PIVOT = 0.
       DO 74  K=1 , KR
         KPIVR = KPIV(K)
         PIVA  = EQUAT(KPIVR,KPIVR)
         IF (ABS(PIVA).GT.ABS(PIVOT)) THEN 
            PIVOT  = PIVA
            KPIVRO = KPIVR
          ENDIF
 74     CONTINUE
!
! normalise pivotal row
!
       KRO = IABS(KHED(KPIVRO))
       IF (ABS(PIVOT).LT.1E-20) THEN
          WRITE(*,476)
  476     FORMAT(43H WARNING-MATRIX SINGULAR OR ILL CONDITIONED)
          WRITE(*,*)'PIVOT-->',PIVOT
          IF (NELL.LT.NE) GOTO 18
!          IF (IGR.LT.KNGRU) GOTO 9
          STOP   
       ENDIF
       DO 80  L=1 , KROW
          QQ(L) = EQUAT(KPIVRO,L)/PIVOT
   80  CONTINUE
       RHS           = RESID(KRO)/PIVOT
       RESID(KRO)    = RHS
       PVKOL(KPIVRO) = PIVOT
!
! eliminate then delete pivotal row and column
!
      IF (KPIVRO.EQ.1) GOTO 104
        KPIVR = KPIVRO - 1
        DO 100  K=1 , KPIVR
           KRW = IABS(KHED(K))
           FAC = EQUAT(K,KPIVRO)
           PVKOL(K) = FAC
           IF (KPIVRO.EQ.1.OR.FAC.EQ.0.) GOTO 88
             LPIVC = KPIVRO - 1
             DO 84  L=1 , LPIVC
               EQUAT(K,L) = EQUAT(K,L) - FAC*QQ(L)
 84          CONTINUE
 88        IF (KPIVRO.EQ.KROW) GOTO 96
             LPIVC = KPIVRO + 1
             DO 92  L=LPIVC , KROW
               EQUAT(K,L-1) = EQUAT(K,L) - FAC*QQ(L)
 92          CONTINUE
 96        RESID(KRW) = RESID(KRW) - FAC*RHS      
 100     CONTINUE
 104   IF (KPIVRO.EQ.KROW) GOTO 128
         KPIVR = KPIVRO + 1
         DO 124  K=KPIVR , KROW
           KRW = IABS(KHED(K))
           FAC = EQUAT(K,KPIVRO)
           PVKOL(K) = FAC
           IF (KPIVRO.EQ.1) GOTO 112
             LPIVC = KPIVRO - 1
             DO 108  L=1 , LPIVC
               EQUAT(K-1,L) = EQUAT(K,L) - FAC*QQ(L)
 108         CONTINUE
 112       IF (KPIVRO.EQ.KROW) GOTO 120
             LPIVC = KPIVRO + 1
             DO 116  L=LPIVC , KROW
               EQUAT(K-1,L-1) = EQUAT(K,L) - FAC*QQ(L)
 116         CONTINUE
 120       RESID(KRW) = RESID(KRW) - FAC*RHS      
 124     CONTINUE
 128   CONTINUE
!
! write pivotal equations on disc
!
       WRITE(ND1)  KRO,KROW,KPIVRO,(QQ(K),KHED(K),K=1,KROW)
!
!       IME=IME+1
!       KMEMO(IME,1)=KRO
!       KMEMO(IME,2)=KROW
!       KMEMO(IME,3)=KPIVRO
!       DO 2001 IMEM=1,KROW
!         KMEMO(IME,3+IMEM) =KHED(IMEM)
!         QQMEM(IME,IMEM)   =QQ(IMEM)
! 2001  CONTINUE
!
       DO 129  K=1 , KROW
          EQUAT(K,KROW) = 0.
          EQUAT(KROW,K) = 0.
  129  CONTINUE
       KR   =KR -1
!
! rearrange heading vectors
!           
       KROW = KROW - 1
       IF (.NOT.(KPIVRO.EQ.KROW+1)) THEN
          DO 140  K=KPIVRO , KROW
             KHED(K) = KHED(K+1)
  140     CONTINUE
       ENDIF
!
! determine wheter to assemble,eliminate or backsubstitute
!
!       WRITE(*,*)' DIFFERENT CH'
       IF (KROW.GT.NCRIT.AND.NPIVO.EQ.0) GOTO 60
       IF (KR.GT.0.AND.NELL.LT.NE.AND.NPIVO.EQ.1) GOTO 60
!       WRITE(*,*)'TEST 1',NELL,NE
       IF (NELL.LT.KNELE) GOTO 18
!       WRITE(*,*)'TEST 2'
!       IF (IGR.LT.KNGRU) GOTO 9
       IF (KROW.GT.1) GOTO 60
       KPIVRO = 1                                
       PIVOT  = EQUAT(1,1)
       KRO    = IABS(KHED(1))
       QQ(1)  = 1.
       IF (ABS(PIVOT).LT.1E-20) THEN
          WRITE(*,476)
          WRITE(*,*)'PIVOT ',PIVOT
       ENDIF
       RESID(KRO) = RESID(KRO)/PIVOT 
       IME=IME+1
!       KMEMO(IME,1)=KRO
!       KMEMO(IME,2)=KROW
!       KMEMO(IME,3)=KPIVRO
!       DO 2002 IMEM=1,1
!         KMEMO(IME,3+IMEM)=KHED(IMEM)
!         QQMEM(IME,IMEM)=QQ(IMEM)
! 2002  CONTINUE
!
       WRITE(ND1)  KRO,KROW,KPIVRO,QQ(1),KHED(1)
!
       CALL BACKSUB(SOLVE,NH,ND1,NCOD,BDISP,RESID, &
                   KMEMO,QQMEM,IME,NDNOD)
!                              
 999  CONTINUE

!      WRITE(*,*)' END--SUB  -> DRSOLV'
      RETURN
      END
!-----------------------------------------------------------------C
! NAME   :  BACKSUB                                               C
!                                                                 C
! PURPOSE:  Back-substitution for full pivoting                   C
!                                                                 C
! CALLS  :                                                        C
!                                                                 C
! INPUT  :                                                        C
!                                                                 C
! OUTPUT :                                                        C
!                                                                 C
!-----------------------------------------------------------------C
!
      SUBROUTINE BACKSUB(SOLVE,NH,ND1,NCOD,BDISP,RESID, &
                         KMEMO,QQMEM,IME,NDNOD)
!
      IMPLICIT REAL*8(A-H,O-Z)
      
!      
!      DIMENSION
!     1  KMEMO(KNDOF,KNuse),QQMEM(KNDOF,KNuse)
!

      integer                  :: ncod(:)
      double precision         :: bdisp(:), resid(:), solve(:)
      
      integer                  :: khed(knuse)
      double precision         :: qq(knuse)
      
!
!      WRITE(*,*)' SUB  -> BACKSUB'
!
!       BACK SUBSTITUTION
!
! NOTE THIS ASSUMES NDNOD DOFS PER NODE
!
! NDOF total number of dofs
!
       NDOF=NH*NDNOD         
       DO 32 IV=NDOF,1,-1
!
!       WRITE(*,*)' READING---->',IV
       BACKSPACE ND1
       READ(ND1)  KRO,KROW,KPIVRO,(QQ(K),KHED(K),K=1,KROW)
       BACKSPACE ND1
!
!
!       KRO   =KMEMO(IV,1)
!       KROW  =KMEMO(IV,2)
!       KPIVRO=KMEMO(IV,3)
!       DO 2001 IMEM=1,KROW
!         KHED(IMEM)=KMEMO(IV,3+IMEM)
!         QQ(IMEM)  =QQMEM(IV,IMEM)
! 2001  CONTINUE
!
          GASH       = 0.
          QQ(KPIVRO) = 0.
          SOLVE(KRO) = 0.
          DO 16  L=1 , KROW
             GASH = GASH - QQ(L)*SOLVE(IABS(KHED(L)))
   16     CONTINUE       
          IF (NCOD(KRO).EQ.2) THEN
            NCOD(KRO)=1
            RESID(KRO)=RESID(KRO)+GASH
            SOLVE(KRO)=BDISP(KRO)
          ELSE
            SOLVE(KRO) = RESID(KRO)+GASH
          ENDIF

          BDISP(KRO)=0.
 32   CONTINUE

       RETURN
       END                

     
end module front_solver
