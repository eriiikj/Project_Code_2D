  O1  h   k820309    ?          18.0        d¹łe                                                                                                          
       /home/er7128ja/Nextcloud/Projekt/Project_Code_2D/src/lib/somelib_ifort/mater_magnkho.f90 MATER_MAGNKHO              MAGMATER_KHOCALCH1 MAGMATER_KHOCALCH2 MAGMATER_KHOCALCH3 MAGMATER_KHOCALCH4 MAGMATER_TANGENTKHO1 MAGMATER_TANGENTKHO2 HAMODELODES PROPEM BSAVE HSAVE GSAVE KSAVE MATID NEW OLD                                                     
       INV3 INVERSE                      @                              
       SOLVEQ                                                        u #MAGMATER_KHOCALCH1    #MAGMATER_KHOCALCH3    #         @     @X                                                 #H    #B    #DT    #IEL              D @                                                  
               &                   &                                                     @                                                  
               &                   &                                                      @                                    
                  @                                           #         @     @X                                                 #HOUT 	   #BIN 
   #DT    #IEL    #IGP              D                                 	                   
 ,              &                                                     
 @                               
                   
 -             &                                                     
                                       
                
                                                       
                                                                                                   u #MAGMATER_TANGENTKHO1    #MAGMATER_TANGENTKHO2    #         @     @X                                                 #MUI    #SI    #DT    #IEL              D@                                                  
 	              &                   &                   &                                                     D @                                                  
 
              &                   &                   &                                                     
  @                                    
                D @                                           #         @     @X                                                #MUI    #SI    #DT    #IEL    #IGP              D@                                                  
               &                   &                                                     D                                                    
               &                   &                                                     
                                       
                                                                                                                                                                      	   u #EQ_SOLVE_NDOF    #EQ_SOLVE    #EQ_SOLVE_RS    #EQ_SOLVE_BC     #EQ_SOLVE_BCN %   #EQ_SOLVE_SPARSE_BC +   #EQ_SOLVE_SPARSE_BCNX 1   #EQ_SOLVE_SPARSE2 7   #EQ_SOLVE_SPARSE :   #         @     @                                               #K    #F                                                                 
 6              &                   &                                                                                                         
 7              &                                           #         @     @                                                #K    #F                                                                 
 <              &                   &                                                                                                        
 =              &                   &                                           #         @     @                                                 #K !   #F "   #BCDOF #   #BCVAL $                                             !                   
 ?              &                   &                                                                                      "                   
 @              &                                                                                     #                    B              &                                                                                      $                   
 A              &                                           #         @     @                            %                    #K &   #F '   #BCNOD (   #BCVAL )   #DOFNOD *                                             &                   
 H              &                   &                                                                                      '                   
 I              &                                                                                     (                    K              &                   &                                                                                      )                   
 J              &                                                                                      *            #         @     @                            +                    #K ,   #F .   #BCDOF /   #BCVAL 0                                              ,     ą               #SPARSE -                                             .                   
 j              &                                                                                      /                    l              &                                                                                      0                   
 k              &                                           #         @     @                            1                    #K 2   #F 3   #BCNOD 4   #BCVAL 5   #DOFNOD 6                                              2     ą               #SPARSE -                                             3                   
 o              &                                                                                     4                    q              &                   &                                                                                      5                   
 p              &                                                                                      6            #         @     @                            7                    #K 8   #F 9                                              8     ą               #SPARSE -                                             9                   
               &                   &                                           #         @     @                            :                    #K ;   #F <                                              ;     ą               #SPARSE -                                             <                   
 R              &                                                             @               D           -     'ą                    #A =   #IA >   #JA ?   #COLS @                                            =                              
            &                                                                                    >            H                             &                                                                                    ?                                         &                                                                                       @     Ų             #         @                                  A                    #RES B   #CA C                                              B                   
               &                   &                                                                                      C                   
               &                   &                                           #         @                                   D                    #A E   #AINVERSE F   #N G                                             E                    
       p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                                                    F                    
       p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                                                     G            #         @                                                      #K H   #F I   #NDOF J                                              H                   
 9              &                   &                                                                                      I                   
 :              &                                                                                      J                                                        K            #         @                                   L                    #PROPEMIN M   #MATIDIN N   #NIP O             @                               M                   
               &                   &                                                     @                               N                                  &                                                                                      O            #         @                                   P                     #         @                                  Q                    #CONV R   #H S   #G T   #K U   #DB V   #HOLD W   #GOLD X   #KOLD Y   #DT Z   #IEL [   #IGP \             D                                 R                      D @                               S                   
 L              &                                                     D @                               T                   
 M              &                                                     D                                 U                   
 N              &                                                     
                                  V                   
 O             &                                                     
                                  W                   
 P             &                                                     
                                  X                   
 Q             &                                                     
                                  Y                   
 R             &                                                     
                                  Z     
                
                                  [                     
                                 \           #         @                                   ]                    #FIELD ^   #HOUT _   #IEL `              @                               ^                     1           D@                               _                   
 f              &                   &                                                                                      `            #         @                                   a                    #PROP b   #IEL c             D                                 b                   
 g              &                                                                                      c                   o      fn#fn #     æ   b   uapp(MATER_MAGNKHO    Ī  M   J  MATRIX_UTIL      G   J  FEM_SYSTEM %   b  p       gen@MAGMATER_CALCKHO #   Ņ  g      MAGMATER_KHOCALCH1 %   9  ¤   a   MAGMATER_KHOCALCH1%H %   Ż  ¤   a   MAGMATER_KHOCALCH1%B &     @   a   MAGMATER_KHOCALCH1%DT '   Į  @   a   MAGMATER_KHOCALCH1%IEL #     u      MAGMATER_KHOCALCH3 (   v     a   MAGMATER_KHOCALCH3%HOUT '        a   MAGMATER_KHOCALCH3%BIN &     @   a   MAGMATER_KHOCALCH3%DT '   Ī  @   a   MAGMATER_KHOCALCH3%IEL '     @   a   MAGMATER_KHOCALCH3%IGP (   N  t       gen@MAGMATER_TANGENTKHO %   Ā  j      MAGMATER_TANGENTKHO1 )   ,  ¼   a   MAGMATER_TANGENTKHO1%MUI (   č  ¼   a   MAGMATER_TANGENTKHO1%SI (   ¤	  @   a   MAGMATER_TANGENTKHO1%DT )   ä	  @   a   MAGMATER_TANGENTKHO1%IEL %   $
  s      MAGMATER_TANGENTKHO2 )   
  ¤   a   MAGMATER_TANGENTKHO2%MUI (   ;  ¤   a   MAGMATER_TANGENTKHO2%SI (   ß  @   a   MAGMATER_TANGENTKHO2%DT )     @   a   MAGMATER_TANGENTKHO2%IEL )   _  @   a   MAGMATER_TANGENTKHO2%IGP &     ņ       gen@SOLVEQ+FEM_SYSTEM $     V      EQ_SOLVE+FEM_SYSTEM &   ē  ¤   a   EQ_SOLVE%K+FEM_SYSTEM &        a   EQ_SOLVE%F+FEM_SYSTEM '     V      EQ_SOLVE_RS+FEM_SYSTEM )   m  ¤   a   EQ_SOLVE_RS%K+FEM_SYSTEM )     ¤   a   EQ_SOLVE_RS%F+FEM_SYSTEM '   µ  l      EQ_SOLVE_BC+FEM_SYSTEM )   !  ¤   a   EQ_SOLVE_BC%K+FEM_SYSTEM )   Å     a   EQ_SOLVE_BC%F+FEM_SYSTEM -   Q     a   EQ_SOLVE_BC%BCDOF+FEM_SYSTEM -   Ż     a   EQ_SOLVE_BC%BCVAL+FEM_SYSTEM (   i  x      EQ_SOLVE_BCN+FEM_SYSTEM *   į  ¤   a   EQ_SOLVE_BCN%K+FEM_SYSTEM *        a   EQ_SOLVE_BCN%F+FEM_SYSTEM .     ¤   a   EQ_SOLVE_BCN%BCNOD+FEM_SYSTEM .   µ     a   EQ_SOLVE_BCN%BCVAL+FEM_SYSTEM /   A  @   a   EQ_SOLVE_BCN%DOFNOD+FEM_SYSTEM .     l      EQ_SOLVE_SPARSE_BC+FEM_SYSTEM 0   ķ  T   a   EQ_SOLVE_SPARSE_BC%K+FEM_SYSTEM 0   A     a   EQ_SOLVE_SPARSE_BC%F+FEM_SYSTEM 4   Ķ     a   EQ_SOLVE_SPARSE_BC%BCDOF+FEM_SYSTEM 4   Y     a   EQ_SOLVE_SPARSE_BC%BCVAL+FEM_SYSTEM 0   å  x      EQ_SOLVE_SPARSE_BCNX+FEM_SYSTEM 2   ]  T   a   EQ_SOLVE_SPARSE_BCNX%K+FEM_SYSTEM 2   ±     a   EQ_SOLVE_SPARSE_BCNX%F+FEM_SYSTEM 6   =  ¤   a   EQ_SOLVE_SPARSE_BCNX%BCNOD+FEM_SYSTEM 6   į     a   EQ_SOLVE_SPARSE_BCNX%BCVAL+FEM_SYSTEM 7   m  @   a   EQ_SOLVE_SPARSE_BCNX%DOFNOD+FEM_SYSTEM ,   ­  V      EQ_SOLVE_SPARSE2+FEM_SYSTEM .     T   a   EQ_SOLVE_SPARSE2%K+FEM_SYSTEM .   W  ¤   a   EQ_SOLVE_SPARSE2%F+FEM_SYSTEM +   ū  V      EQ_SOLVE_SPARSE+FEM_SYSTEM -   Q  T   a   EQ_SOLVE_SPARSE%K+FEM_SYSTEM -   „     a   EQ_SOLVE_SPARSE%F+FEM_SYSTEM #   1  q       SPARSE+SPARSE_UTIL %   ¢     a   SPARSE%A+SPARSE_UTIL &   6     a   SPARSE%IA+SPARSE_UTIL &   Ź     a   SPARSE%JA+SPARSE_UTIL (   ^   H   a   SPARSE%COLS+SPARSE_UTIL !   ¦   Y       INV3+MATRIX_UTIL %   ’   ¤   a   INV3%RES+MATRIX_UTIL $   £!  ¤   a   INV3%CA+MATRIX_UTIL $   G"  d       INVERSE+MATRIX_UTIL &   «"  ü   a   INVERSE%A+MATRIX_UTIL -   §#  ü   a   INVERSE%AINVERSE+MATRIX_UTIL &   £$  @   a   INVERSE%N+MATRIX_UTIL )   ć$  `       EQ_SOLVE_NDOF+FEM_SYSTEM +   C%  ¤   a   EQ_SOLVE_NDOF%K+FEM_SYSTEM +   ē%     a   EQ_SOLVE_NDOF%F+FEM_SYSTEM .   s&  @   a   EQ_SOLVE_NDOF%NDOF+FEM_SYSTEM !   ³&  @       DEBUGG_MATER_KHO !   ó&  l       MAGMATER_INITKHO *   _'  ¤   a   MAGMATER_INITKHO%PROPEMIN )   (     a   MAGMATER_INITKHO%MATIDIN %   (  @   a   MAGMATER_INITKHO%NIP #   Ļ(  H       MAGMATER_ACCEPTKHO &   )  §       MAGMATER_KHOCALCHSTEP +   ¾)  @   a   MAGMATER_KHOCALCHSTEP%CONV (   ž)     a   MAGMATER_KHOCALCHSTEP%H (   *     a   MAGMATER_KHOCALCHSTEP%G (   +     a   MAGMATER_KHOCALCHSTEP%K )   ¢+     a   MAGMATER_KHOCALCHSTEP%DB +   .,     a   MAGMATER_KHOCALCHSTEP%HOLD +   ŗ,     a   MAGMATER_KHOCALCHSTEP%GOLD +   F-     a   MAGMATER_KHOCALCHSTEP%KOLD )   Ņ-  @   a   MAGMATER_KHOCALCHSTEP%DT *   .  @   a   MAGMATER_KHOCALCHSTEP%IEL *   R.  @   a   MAGMATER_KHOCALCHSTEP%IGP !   .  f       MAGMATER_GETKHOH '   ų.  L   a   MAGMATER_GETKHOH%FIELD &   D/  ¤   a   MAGMATER_GETKHOH%HOUT %   č/  @   a   MAGMATER_GETKHOH%IEL %   (0  [       MAGMATER_KHOPARA_GET *   0     a   MAGMATER_KHOPARA_GET%PROP )   1  @   a   MAGMATER_KHOPARA_GET%IEL 