  6*  X   k820309    ?          18.0        İe                                                                                                          
       /home/er7128ja/Nextcloud/Projekt/Project_Code_2D/src/lib/somelib_ifort/mater_magnjav.f90 MATER_MAGNJAV              MAGMATER_JAVCALCH1 MAGMATER_JAVCALCH2 MAGMATER_JAVCALCNR MATER_CALCREV MATER_CALCIRR MATER_CALCREVEULER MATER_JAVNRSTEP MAGMATER_JAVMATER_GETH1 MAGMATER_JAVMATER_GETH2 MAGMATER_JAVTANGENT PROPEM MATID MU0 BSAVE HSAVE NEW OLD                                                     
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
 J              &                                                                                      &            #         @     @                            '                    #K (   #F *   #BCDOF +   #BCVAL ,                                              (     à               #SPARSE )                                             *                   
 j              &                                                                                      +                    l              &                                                                                      ,                   
 k              &                                           #         @     @                            -                    #K .   #F /   #BCNOD 0   #BCVAL 1   #DOFNOD 2                                              .     à               #SPARSE )                                             /                   
 o              &                                                                                     0                    q              &                   &                                                                                      1                   
 p              &                                                                                      2            #         @     @                            3                    #K 4   #F 5                                              4     à               #SPARSE )                                             5                   
               &                   &                                           #         @     @                            6                    #K 7   #F 8                                              7     à               #SPARSE )                                             8                   
 R              &                                                             @               D           )     'à                    #A 9   #IA :   #JA ;   #COLS <                                            9                              
            &                                                                                    :            H                             &                                                                                    ;                                         &                                                                                       <     Ø             #         @                                  =                    #RES >   #CA ?                                              >                   
               &                   &                                                                                      ?                   
               &                   &                                           #         @                                                      #K @   #F A   #NDOF B                                              @                   
 9              &                   &                                                                                      A                   
 :              &                                                                                      B            #         @                                   C                    #PROPEMIN D   #MATIDIN E   #NIP F             @                               D                   
               &                   &                                                     @                               E                                  &                                                                                      F            #         @                                   G                     #         @                                  H                    #MUI I   #SI J   #IEL K   #IGP L             D@                               I                   
 j              &                   &                                                     D                                 J                   
 k              &                   &                                                                                      K                                                       L            #         @                                   M                    #FIELD N   #HOUT O   #IEL P              @                               N                     1           D@                               O                   
 º              &                   &                                                                                      P            #         @                                   Q                    #PROP R   #IEL S             D                                 R                   
 »              &                                                                                      S                   o      fn#fn #     ñ   b   uapp(MATER_MAGNJAV       E   J  MATRIX_UTIL    E  G   J  FEM_SYSTEM %     p       gen@MAGMATER_CALCJAV #   ü  _      MAGMATER_JAVCALCH1 %   [  ¤   a   MAGMATER_JAVCALCH1%H %   ÿ  ¤   a   MAGMATER_JAVCALCH1%B '   £  @   a   MAGMATER_JAVCALCH1%IEL #   ã  h      MAGMATER_JAVCALCNR %   K     a   MAGMATER_JAVCALCNR%H %   ×     a   MAGMATER_JAVCALCNR%B '   c  @   a   MAGMATER_JAVCALCNR%IEL '   £  @   a   MAGMATER_JAVCALCNR%IGP (   ã  z       gen@MAGMATER_TANGENTJAV (   ]  b      MAGMATER_JAVMATER_GETH1 ,   ¿  ¼   a   MAGMATER_JAVMATER_GETH1%MUI +   {  ¼   a   MAGMATER_JAVMATER_GETH1%SI ,   7	  @   a   MAGMATER_JAVMATER_GETH1%IEL (   w	  k      MAGMATER_JAVMATER_GETH2 ,   â	  ¤   a   MAGMATER_JAVMATER_GETH2%MUI +   
  ¤   a   MAGMATER_JAVMATER_GETH2%SI ,   *  @   a   MAGMATER_JAVMATER_GETH2%IEL ,   j  @   a   MAGMATER_JAVMATER_GETH2%IGP &   ª  ò       gen@SOLVEQ+FEM_SYSTEM $     V      EQ_SOLVE+FEM_SYSTEM &   ò  ¤   a   EQ_SOLVE%K+FEM_SYSTEM &        a   EQ_SOLVE%F+FEM_SYSTEM '   "  V      EQ_SOLVE_RS+FEM_SYSTEM )   x  ¤   a   EQ_SOLVE_RS%K+FEM_SYSTEM )     ¤   a   EQ_SOLVE_RS%F+FEM_SYSTEM '   À  l      EQ_SOLVE_BC+FEM_SYSTEM )   ,  ¤   a   EQ_SOLVE_BC%K+FEM_SYSTEM )   Ğ     a   EQ_SOLVE_BC%F+FEM_SYSTEM -   \     a   EQ_SOLVE_BC%BCDOF+FEM_SYSTEM -   è     a   EQ_SOLVE_BC%BCVAL+FEM_SYSTEM (   t  x      EQ_SOLVE_BCN+FEM_SYSTEM *   ì  ¤   a   EQ_SOLVE_BCN%K+FEM_SYSTEM *        a   EQ_SOLVE_BCN%F+FEM_SYSTEM .     ¤   a   EQ_SOLVE_BCN%BCNOD+FEM_SYSTEM .   À     a   EQ_SOLVE_BCN%BCVAL+FEM_SYSTEM /   L  @   a   EQ_SOLVE_BCN%DOFNOD+FEM_SYSTEM .     l      EQ_SOLVE_SPARSE_BC+FEM_SYSTEM 0   ø  T   a   EQ_SOLVE_SPARSE_BC%K+FEM_SYSTEM 0   L     a   EQ_SOLVE_SPARSE_BC%F+FEM_SYSTEM 4   Ø     a   EQ_SOLVE_SPARSE_BC%BCDOF+FEM_SYSTEM 4   d     a   EQ_SOLVE_SPARSE_BC%BCVAL+FEM_SYSTEM 0   ğ  x      EQ_SOLVE_SPARSE_BCNX+FEM_SYSTEM 2   h  T   a   EQ_SOLVE_SPARSE_BCNX%K+FEM_SYSTEM 2   ¼     a   EQ_SOLVE_SPARSE_BCNX%F+FEM_SYSTEM 6   H  ¤   a   EQ_SOLVE_SPARSE_BCNX%BCNOD+FEM_SYSTEM 6   ì     a   EQ_SOLVE_SPARSE_BCNX%BCVAL+FEM_SYSTEM 7   x  @   a   EQ_SOLVE_SPARSE_BCNX%DOFNOD+FEM_SYSTEM ,   ¸  V      EQ_SOLVE_SPARSE2+FEM_SYSTEM .     T   a   EQ_SOLVE_SPARSE2%K+FEM_SYSTEM .   b  ¤   a   EQ_SOLVE_SPARSE2%F+FEM_SYSTEM +     V      EQ_SOLVE_SPARSE+FEM_SYSTEM -   \  T   a   EQ_SOLVE_SPARSE%K+FEM_SYSTEM -   °     a   EQ_SOLVE_SPARSE%F+FEM_SYSTEM #   <  q       SPARSE+SPARSE_UTIL %   ­     a   SPARSE%A+SPARSE_UTIL &   A     a   SPARSE%IA+SPARSE_UTIL &   Õ     a   SPARSE%JA+SPARSE_UTIL (   i  H   a   SPARSE%COLS+SPARSE_UTIL !   ±  Y       INV3+MATRIX_UTIL %   
   ¤   a   INV3%RES+MATRIX_UTIL $   ®   ¤   a   INV3%CA+MATRIX_UTIL )   R!  `       EQ_SOLVE_NDOF+FEM_SYSTEM +   ²!  ¤   a   EQ_SOLVE_NDOF%K+FEM_SYSTEM +   V"     a   EQ_SOLVE_NDOF%F+FEM_SYSTEM .   â"  @   a   EQ_SOLVE_NDOF%NDOF+FEM_SYSTEM !   "#  l       MAGMATER_INITJAV *   #  ¤   a   MAGMATER_INITJAV%PROPEMIN )   2$     a   MAGMATER_INITJAV%MATIDIN %   ¾$  @   a   MAGMATER_INITJAV%NIP #   ş$  H       MAGMATER_ACCEPTJAV ,   F%  k       MAGMATER_JAVMATER_TANGENTNR 0   ±%  ¤   a   MAGMATER_JAVMATER_TANGENTNR%MUI /   U&  ¤   a   MAGMATER_JAVMATER_TANGENTNR%SI 0   ù&  @   a   MAGMATER_JAVMATER_TANGENTNR%IEL 0   9'  @   a   MAGMATER_JAVMATER_TANGENTNR%IGP !   y'  f       MAGMATER_GETJAVH '   ß'  L   a   MAGMATER_GETJAVH%FIELD &   +(  ¤   a   MAGMATER_GETJAVH%HOUT %   Ï(  @   a   MAGMATER_GETJAVH%IEL %   )  [       MAGMATER_JAVPARA_GET *   j)     a   MAGMATER_JAVPARA_GET%PROP )   ö)  @   a   MAGMATER_JAVPARA_GET%IEL 