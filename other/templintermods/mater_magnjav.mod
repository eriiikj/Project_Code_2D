  9*  X   k820309    ?          18.0        8fłe                                                                                                          
       /home/er7128ja/Nextcloud/Projekt/Project_Code_2D/src/lib/somelib_gfortran/mater_magnjav.f90 MATER_MAGNJAV              MAGMATER_JAVCALCH1 MAGMATER_JAVCALCH2 MAGMATER_JAVCALCNR MATER_CALCREV MATER_CALCIRR MATER_CALCREVEULER MATER_JAVNRSTEP MAGMATER_JAVMATER_GETH1 MAGMATER_JAVMATER_GETH2 MAGMATER_JAVTANGENT PROPEM MATID MU0 BSAVE HSAVE NEW OLD                                                     
       INV3                      @                              
       SOLVEQ                                                        u #MAGMATER_JAVCALCH1    #MAGMATER_JAVCALCNR    #         @     @X                                                 #H    #B    #IEL              D @                                                  
 	              &                   &                                                     D@                                                  
 
              &                   &                                                     D @                                           #         @     @X                                                #H    #B 	   #IEL 
   #IGP              D                                                    
               &                                                     D@                               	                   
               &                                                     D @                               
                      D @                                                                                                  u #MAGMATER_JAVMATER_GETH1    #MAGMATER_JAVMATER_GETH2    #         @     @X                                                 #MUI    #SI    #IEL              D@                                                  
               &                   &                   &                                                     D @                                                  
               &                   &                   &                                                     D @                                           #         @     @X                                                 #MUI    #SI    #IEL    #IGP              D@                                                  
               &                   &                                                     D                                                    
               &                   &                                                                                                             @                                                                                             	   u #EQ_SOLVE_NDOF    #EQ_SOLVE    #EQ_SOLVE_RS    #EQ_SOLVE_BC    #EQ_SOLVE_BCN !   #EQ_SOLVE_SPARSE_BC '   #EQ_SOLVE_SPARSE_BCNX -   #EQ_SOLVE_SPARSE2 3   #EQ_SOLVE_SPARSE 6   #         @     @                                               #K    #F                                                                 
 6              &                   &                                                                                                         
 7              &                                           #         @     @                                                #K    #F                                                                 
 <              &                   &                                                                                                        
 =              &                   &                                           #         @     @                                                #K    #F    #BCDOF    #BCVAL                                                                  
 ?              &                   &                                                                                                         
 @              &                                                                                                         B              &                                                                                                          
 A              &                                           #         @     @                            !                    #K "   #F #   #BCNOD $   #BCVAL %   #DOFNOD &                                             "                   
 H              &                   &                                                                                      #                   
 I              &                                                                                     $                    K              &                   &                                                                                      %                   
 J              &                                                                                      &            #         @     @                            '                    #K (   #F *   #BCDOF +   #BCVAL ,                                              (     ą               #SPARSE )                                             *                   
 j              &                                                                                      +                    l              &                                                                                      ,                   
 k              &                                           #         @     @                            -                    #K .   #F /   #BCNOD 0   #BCVAL 1   #DOFNOD 2                                              .     ą               #SPARSE )                                             /                   
 o              &                                                                                     0                    q              &                   &                                                                                      1                   
 p              &                                                                                      2            #         @     @                            3                    #K 4   #F 5                                              4     ą               #SPARSE )                                             5                   
               &                   &                                           #         @     @                            6                    #K 7   #F 8                                              7     ą               #SPARSE )                                             8                   
 R              &                                                             @               D           )     'ą                    #A 9   #IA :   #JA ;   #COLS <                                            9                              
            &                                                                                    :            H                             &                                                                                    ;                                         &                                                                                       <     Ų             #         @                                  =                    #RES >   #CA ?                                              >                   
               &                   &                                                                                      ?                   
               &                   &                                           #         @                                                      #K @   #F A   #NDOF B                                              @                   
 9              &                   &                                                                                      A                   
 :              &                                                                                      B            #         @                                   C                    #PROPEMIN D   #MATIDIN E   #NIP F             @                               D                   
               &                   &                                                     @                               E                                  &                                                                                      F            #         @                                   G                     #         @                                  H                    #MUI I   #SI J   #IEL K   #IGP L             D@                               I                   
 j              &                   &                                                     D                                 J                   
 k              &                   &                                                                                      K                                                       L            #         @                                   M                    #FIELD N   #HOUT O   #IEL P              @                               N                     1           D@                               O                   
 ŗ              &                   &                                                                                      P            #         @                                   Q                    #PROP R   #IEL S             D                                 R                   
 »              &                                                                                      S                   r      fn#fn #     ń   b   uapp(MATER_MAGNJAV      E   J  MATRIX_UTIL    H  G   J  FEM_SYSTEM %     p       gen@MAGMATER_CALCJAV #   ’  _      MAGMATER_JAVCALCH1 %   ^  ¤   a   MAGMATER_JAVCALCH1%H %     ¤   a   MAGMATER_JAVCALCH1%B '   ¦  @   a   MAGMATER_JAVCALCH1%IEL #   ę  h      MAGMATER_JAVCALCNR %   N     a   MAGMATER_JAVCALCNR%H %   Ś     a   MAGMATER_JAVCALCNR%B '   f  @   a   MAGMATER_JAVCALCNR%IEL '   ¦  @   a   MAGMATER_JAVCALCNR%IGP (   ę  z       gen@MAGMATER_TANGENTJAV (   `  b      MAGMATER_JAVMATER_GETH1 ,   Ā  ¼   a   MAGMATER_JAVMATER_GETH1%MUI +   ~  ¼   a   MAGMATER_JAVMATER_GETH1%SI ,   :	  @   a   MAGMATER_JAVMATER_GETH1%IEL (   z	  k      MAGMATER_JAVMATER_GETH2 ,   å	  ¤   a   MAGMATER_JAVMATER_GETH2%MUI +   
  ¤   a   MAGMATER_JAVMATER_GETH2%SI ,   -  @   a   MAGMATER_JAVMATER_GETH2%IEL ,   m  @   a   MAGMATER_JAVMATER_GETH2%IGP &   ­  ņ       gen@SOLVEQ+FEM_SYSTEM $     V      EQ_SOLVE+FEM_SYSTEM &   õ  ¤   a   EQ_SOLVE%K+FEM_SYSTEM &        a   EQ_SOLVE%F+FEM_SYSTEM '   %  V      EQ_SOLVE_RS+FEM_SYSTEM )   {  ¤   a   EQ_SOLVE_RS%K+FEM_SYSTEM )     ¤   a   EQ_SOLVE_RS%F+FEM_SYSTEM '   Ć  l      EQ_SOLVE_BC+FEM_SYSTEM )   /  ¤   a   EQ_SOLVE_BC%K+FEM_SYSTEM )   Ó     a   EQ_SOLVE_BC%F+FEM_SYSTEM -   _     a   EQ_SOLVE_BC%BCDOF+FEM_SYSTEM -   ė     a   EQ_SOLVE_BC%BCVAL+FEM_SYSTEM (   w  x      EQ_SOLVE_BCN+FEM_SYSTEM *   ļ  ¤   a   EQ_SOLVE_BCN%K+FEM_SYSTEM *        a   EQ_SOLVE_BCN%F+FEM_SYSTEM .     ¤   a   EQ_SOLVE_BCN%BCNOD+FEM_SYSTEM .   Ć     a   EQ_SOLVE_BCN%BCVAL+FEM_SYSTEM /   O  @   a   EQ_SOLVE_BCN%DOFNOD+FEM_SYSTEM .     l      EQ_SOLVE_SPARSE_BC+FEM_SYSTEM 0   ū  T   a   EQ_SOLVE_SPARSE_BC%K+FEM_SYSTEM 0   O     a   EQ_SOLVE_SPARSE_BC%F+FEM_SYSTEM 4   Ū     a   EQ_SOLVE_SPARSE_BC%BCDOF+FEM_SYSTEM 4   g     a   EQ_SOLVE_SPARSE_BC%BCVAL+FEM_SYSTEM 0   ó  x      EQ_SOLVE_SPARSE_BCNX+FEM_SYSTEM 2   k  T   a   EQ_SOLVE_SPARSE_BCNX%K+FEM_SYSTEM 2   æ     a   EQ_SOLVE_SPARSE_BCNX%F+FEM_SYSTEM 6   K  ¤   a   EQ_SOLVE_SPARSE_BCNX%BCNOD+FEM_SYSTEM 6   ļ     a   EQ_SOLVE_SPARSE_BCNX%BCVAL+FEM_SYSTEM 7   {  @   a   EQ_SOLVE_SPARSE_BCNX%DOFNOD+FEM_SYSTEM ,   »  V      EQ_SOLVE_SPARSE2+FEM_SYSTEM .     T   a   EQ_SOLVE_SPARSE2%K+FEM_SYSTEM .   e  ¤   a   EQ_SOLVE_SPARSE2%F+FEM_SYSTEM +   	  V      EQ_SOLVE_SPARSE+FEM_SYSTEM -   _  T   a   EQ_SOLVE_SPARSE%K+FEM_SYSTEM -   ³     a   EQ_SOLVE_SPARSE%F+FEM_SYSTEM #   ?  q       SPARSE+SPARSE_UTIL %   °     a   SPARSE%A+SPARSE_UTIL &   D     a   SPARSE%IA+SPARSE_UTIL &   Ų     a   SPARSE%JA+SPARSE_UTIL (   l  H   a   SPARSE%COLS+SPARSE_UTIL !   “  Y       INV3+MATRIX_UTIL %      ¤   a   INV3%RES+MATRIX_UTIL $   ±   ¤   a   INV3%CA+MATRIX_UTIL )   U!  `       EQ_SOLVE_NDOF+FEM_SYSTEM +   µ!  ¤   a   EQ_SOLVE_NDOF%K+FEM_SYSTEM +   Y"     a   EQ_SOLVE_NDOF%F+FEM_SYSTEM .   å"  @   a   EQ_SOLVE_NDOF%NDOF+FEM_SYSTEM !   %#  l       MAGMATER_INITJAV *   #  ¤   a   MAGMATER_INITJAV%PROPEMIN )   5$     a   MAGMATER_INITJAV%MATIDIN %   Į$  @   a   MAGMATER_INITJAV%NIP #   %  H       MAGMATER_ACCEPTJAV ,   I%  k       MAGMATER_JAVMATER_TANGENTNR 0   “%  ¤   a   MAGMATER_JAVMATER_TANGENTNR%MUI /   X&  ¤   a   MAGMATER_JAVMATER_TANGENTNR%SI 0   ü&  @   a   MAGMATER_JAVMATER_TANGENTNR%IEL 0   <'  @   a   MAGMATER_JAVMATER_TANGENTNR%IGP !   |'  f       MAGMATER_GETJAVH '   ā'  L   a   MAGMATER_GETJAVH%FIELD &   .(  ¤   a   MAGMATER_GETJAVH%HOUT %   Ņ(  @   a   MAGMATER_GETJAVH%IEL %   )  [       MAGMATER_JAVPARA_GET *   m)     a   MAGMATER_JAVPARA_GET%PROP )   ł)  @   a   MAGMATER_JAVPARA_GET%IEL 