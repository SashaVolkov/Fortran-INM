	  H  ¦   k820309    q          13.1        Ë@V                                                                                                       
       Schemes.f90 SCHEMES                                                
                                                      
                                                      
                                                      
                     @                                'X            #LX    #LT    #VELOCITY    #DX 	   #DT 
   #GAMMA    #STEPSX    #STEPST    #XSIZE    #NS    #NF    #NP    #ID    #BSTEP    #FSTEP    #INIT                                                       
                                                     
                                                     
                                          	           
                                          
            
                                               (      
                                                0                                                      4                                                      8   	                                                  <   
                                                  @                                                      D                                                      H                                                     L                                                     P         1     À                                              #GRID_INIT    #     @                                                    
   #THIS    #LX    #LT    #STEPSX    #STEPST    #BSTEP    #FSTEP    #VELOCITY    #NP    #ID                                                X       #GRID          
                                      
        
                                      
        
                                               
                                               
                                               
                                               
                                      
        
                                               
                                                        @               @           !     ',            #NS "   #NF #   #D $   #INIT %   #DEINIT *   #EQ .                                          "                                                      #                                                   $                 
        &                       1     À                                %              #FUNC_INIT &   #     @                                   &                    #THIS '   #NS (   #NF )                                         '     ,       #FUNC !         
                                 (             
                                 )       1     À                                *              #FUNC_DEINIT +   #     @                                   +                   #FUNC_DEINIT%ALLOCATED ,   #THIS -                                           ,     ALLOCATED                                       -     ,       #FUNC !   1     À                               .              #FUNC_EQ /   #     @                                  /                    #THIS 0   #THAT 1                                         0     ,       #FUNC !                                         1     ,       #FUNC !                 @                           2     '            #IER 3   #STATUS 4   #FMAX 5   #FMIN 6   #FMED 7   #ABSF 8   #ABSFMAX_L 9   #PULSESTRETCHING :   #WIDTH_LAST ;   #PARAM <   #W_COUNTER =   #TEMP >   #N1_MASS ?   #N1_ALL_MASS_L @   #N2_MASS A   #N2_ALL_MASS_L B   #MESSAGE C   #TO_PRINT H   #SCHEMEPARAM Q   #BORNPARAM W   #INIT ^                                          3                                                      4                   p      p        p                                                     5           
                                          6            
                                          7     (      
                                          8     0      
                                          9     8      
                                          :     @      
                                          ;     H   	   
                                          <     P   
   
                                          =     X      
                                          >     `      
                                          ?     h      
                                          @     p      
                                          A     x      
                                          B           
   1     À                               C              #MET_MESSAGE D   #     @                                  D                   #GRID%FSTEP    #GRID%NF    #GRID%BSTEP    #GRID%NS    #THIS E   #TRANS_MASS F   #G G                                         E            #MET 2        
                                F            
      5 8 O p    U        5 8 O p    U          &  5 8 O p    U        5 8 O p    U         5 8 O p    U        5 8 O p    U               5 8 O p    U        5 8 O p    U         5 8 O p    U        5 8 O p    U        p                                                    G     X       #GRID    1     À                                H              #MET_PRINT I   #     @                                   I                   #MET_PRINT%REAL J   #THIS K   #F L   #G M   #T N   #NAME O   #TIMESET P                            @              J     REAL                                       K            #MET 2                                         L     ,       #FUNC !                                         M     X       #GRID          
                                 N             
                                 O     (                
                                 P       1     À                                Q              #MET_SCHEMEPARAM R   #     @                                   R                    #THIS S   #F T   #G U   #NAME V                                         S            #MET 2                                         T     ,       #FUNC !                                         U     X       #GRID          
                                 V               1     À                                W              #MET_FUNC_BORNPARAM X   #     @                                   X                   #MET_FUNC_BORNPARAM%MAXVAL Y   #MET_FUNC_BORNPARAM%MINVAL Z   #THIS [   #F \   #G ]                                           Y     MAXVAL                                         Z     MINVAL                                       [            #MET 2                                         \     ,       #FUNC !                                         ]     X       #GRID    1     À                                ^              #MET_INIT _   #     @                                   _                    #THIS `   #G a   #STATUS b   #IER c                                         `            #MET 2                                         a     X       #GRID          
                                 b               p      p        p                    
                                 c                                                d                                      5#     @                                 e    	               #V0 f   #IERROR g                                         f                                              g                                                 h                    
      D            1140850688              @               @           i     '            #CASENUMB j   #RCASE k   #EQVTYPE l   #K1 m   #K2 n   #K3 o   #K4 p   #ANGLE q   #CROSS w   #LAX ~   #LAXVEN    #RUNGEKUTTA    #FRUNGE    #INIT    #DEINIT ¢                                          j                                                      k                                                     l                                                   m                 
        &                                                            n        0         
        &                                                            o        T         
        &                                                            p        x         
        &                       1     À    $                            q              #SCH_SCHEME_ANGLE r   #     @                                   r                    #THIS s   #Y t   #YPREV v   #G u                                         s            #SCH i        
D                                t            
      5 8  p    r#GRID     u   U        5 8  p    r#GRID     u   U          &  5 8  p    r#GRID     u   U        5 8  p    r#GRID     u   U        5 8  p    r#GRID     u   U              5 8  p    r#GRID     u   U         5 8  p    r#GRID     u   U        5 8  p    r#GRID     u   U        p                   
                                v            
      5 8  p    r#GRID     u   U        5 8  p    r#GRID     u   U          &  5 8  p    r#GRID     u   U        5 8  p    r#GRID     u   U        5 8  p    r#GRID     u   U              5 8  p    r#GRID     u   U         5 8  p    r#GRID     u   U        5 8  p    r#GRID     u   U        p                                                   u     X       #GRID    1     À                               w         	     #SCH_SCHEME_CROSS x   #     @                                  x                    #THIS y   #Y z   #YPREV |   #YNEXT }   #G {                                         y            #SCH i        
                                z            
      5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U          &  5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U         5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U               5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U         5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U        p                   
                                |            
      5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U          &  5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U         5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U               5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U         5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U        p                   
D                                }            
 	     5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U          &  5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U         5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U               5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U         5 8  p    r#GRID     {   U        5 8  p    r#GRID     {   U        p                                                   {     X       #GRID    1     À                               ~         
     #SCH_SCHEME_LAX    #     @                                                      #THIS    #Y    #YPREV    #G                                                      #SCH i        
D                                            
 
     5 8  p    r#GRID        U        5 8  p    r#GRID        U          &  5 8  p    r#GRID        U        5 8  p    r#GRID        U         5 8  p    r#GRID        U        5 8  p    r#GRID        U               5 8  p    r#GRID        U        5 8  p    r#GRID        U         5 8  p    r#GRID        U        5 8  p    r#GRID        U        p                   
                                            
      5 8  p    r#GRID        U        5 8  p    r#GRID        U          &  5 8  p    r#GRID        U        5 8  p    r#GRID        U         5 8  p    r#GRID        U        5 8  p    r#GRID        U               5 8  p    r#GRID        U        5 8  p    r#GRID        U         5 8  p    r#GRID        U        5 8  p    r#GRID        U        p                                                        X       #GRID    1     À                                              #SCH_SCHEME_LAXVEN    #     @                                                      #SCH_SCHEME_LAXVEN%MOD    #THIS    #F    #FPREV    #FNEXT    #G    #M    #T                                                MOD       D @                                          #SCH i         D @                                   ,       #FUNC !         D @                                   ,       #FUNC !         D @                                   ,       #FUNC !         D @                                   X       #GRID          D @                                          #MET 2          @                                      1     À                                              #SCH_RUNGEKUTTA    #     @                                                       #THIS    #F    #FPREV    #G    #M          D @                                          #SCH i         D                                     ,       #FUNC !                                              ,       #FUNC !         D @                                   X       #GRID          D @                                          #MET 2   1     À                                             #SCH_FRUNGE    #     @                                                      #THIS    #K    #Y    #G                                                      #SCH i        
D                                            
      5 8  p    r#GRID        U        5 8  p    r#GRID        U          &  5 8  p    r#GRID        U        5 8  p    r#GRID        U         5 8  p    r#GRID        U        5 8  p    r#GRID        U               5 8  p    r#GRID        U        5 8  p    r#GRID        U         5 8  p    r#GRID        U        5 8  p    r#GRID        U        p                   
                                             
     5 8  p    r#GRID        U        5 8  p    r#GRID        U          &  5 8  p    r#GRID        U        5 8  p    r#GRID        U         5 8  p    r#GRID        U        5 8  p    r#GRID        U               5 8  p    r#GRID        U        5 8  p    r#GRID        U         5 8  p    r#GRID        U        5 8  p    r#GRID        U        p                                                        X       #GRID    1     À                                              #SCH_INIT    #     @                                                       #THIS    #CASENUMB    #EQVTYPE    #RCASE     #G ¡         D                                            #SCH i         
                                              
                                              
                                                                               ¡     X       #GRID    1     À                                ¢              #SCH_DEINIT £   #     @                                   £                   #SCH_DEINIT%ALLOCATED ¤   #THIS ¥                                          ¤     ALLOCATED       D                                ¥            #SCH i               fn#fn    ¸   <   J   MODNET    ô   <   J   MODFUNC    0  <   J   MPI    l  <   J   METHOD    ¨  à       GRID+MODNET      @   a   GRID%LX+MODNET    È  @   a   GRID%LT+MODNET %     @   a   GRID%VELOCITY+MODNET    H  @   a   GRID%DX+MODNET      @   a   GRID%DT+MODNET "   È  @   a   GRID%GAMMA+MODNET #     @   a   GRID%STEPSX+MODNET #   H  @   a   GRID%STEPST+MODNET "     @   a   GRID%XSIZE+MODNET    È  @   a   GRID%NS+MODNET      @   a   GRID%NF+MODNET    H  @   a   GRID%NP+MODNET      @   a   GRID%ID+MODNET "   È  @   a   GRID%BSTEP+MODNET "     @   a   GRID%FSTEP+MODNET !   H  O   a   GRID%INIT+MODNET !     ª       GRID_INIT+MODNET &   A  F   a   GRID_INIT%THIS+MODNET $     8   a   GRID_INIT%LX+MODNET $   ¿  8   a   GRID_INIT%LT+MODNET (   ÷  8   a   GRID_INIT%STEPSX+MODNET (   /  8   a   GRID_INIT%STEPST+MODNET '   g  8   a   GRID_INIT%BSTEP+MODNET '     8   a   GRID_INIT%FSTEP+MODNET *   ×  8   a   GRID_INIT%VELOCITY+MODNET $   	  8   a   GRID_INIT%NP+MODNET $   G	  8   a   GRID_INIT%ID+MODNET    	  y       FUNC+MODFUNC     ø	  @   a   FUNC%NS+MODFUNC     8
  @   a   FUNC%NF+MODFUNC    x
  l   a   FUNC%D+MODFUNC "   ä
  O   a   FUNC%INIT+MODFUNC "   3  ^       FUNC_INIT+MODFUNC '     F   a   FUNC_INIT%THIS+MODFUNC %   ×  8   a   FUNC_INIT%NS+MODFUNC %     8   a   FUNC_INIT%NF+MODFUNC $   G  Q   a   FUNC%DEINIT+MODFUNC $     i       FUNC_DEINIT+MODFUNC .     >      FUNC_DEINIT%ALLOCATED+MODFUNC )   ?  F   a   FUNC_DEINIT%THIS+MODFUNC       M   a   FUNC%EQ+MODFUNC     Ò  X       FUNC_EQ+MODFUNC %   *  F   a   FUNC_EQ%THIS+MODFUNC %   p  F   a   FUNC_EQ%THAT+MODFUNC    ¶  ^      MET+METHOD      @   a   MET%IER+METHOD "   T  x   a   MET%STATUS+METHOD     Ì  @   a   MET%FMAX+METHOD       @   a   MET%FMIN+METHOD     L  @   a   MET%FMED+METHOD       @   a   MET%ABSF+METHOD %   Ì  @   a   MET%ABSFMAX_L+METHOD +     @   a   MET%PULSESTRETCHING+METHOD &   L  @   a   MET%WIDTH_LAST+METHOD !     @   a   MET%PARAM+METHOD %   Ì  @   a   MET%W_COUNTER+METHOD       @   a   MET%TEMP+METHOD #   L  @   a   MET%N1_MASS+METHOD )     @   a   MET%N1_ALL_MASS_L+METHOD #   Ì  @   a   MET%N2_MASS+METHOD )     @   a   MET%N2_ALL_MASS_L+METHOD #   L  Q   a   MET%MESSAGE+METHOD #            MET_MESSAGE+METHOD (   <  E   a   MET_MESSAGE%THIS+METHOD .     Î  a   MET_MESSAGE%TRANS_MASS+METHOD %   O  F   a   MET_MESSAGE%G+METHOD $     O   a   MET%TO_PRINT+METHOD !   ä         MET_PRINT+METHOD &   r  9      MET_PRINT%REAL+METHOD &   «  E   a   MET_PRINT%THIS+METHOD #   ð  F   a   MET_PRINT%F+METHOD #   6  F   a   MET_PRINT%G+METHOD #   |  8   a   MET_PRINT%T+METHOD &   ´  @   a   MET_PRINT%NAME+METHOD )   ô  8   a   MET_PRINT%TIMESET+METHOD '   ,  U   a   MET%SCHEMEPARAM+METHOD '     f       MET_SCHEMEPARAM+METHOD ,   ç  E   a   MET_SCHEMEPARAM%THIS+METHOD )   ,  F   a   MET_SCHEMEPARAM%F+METHOD )   r  F   a   MET_SCHEMEPARAM%G+METHOD ,   ¸  @   a   MET_SCHEMEPARAM%NAME+METHOD %   ø  X   a   MET%BORNPARAM+METHOD *   P         MET_FUNC_BORNPARAM+METHOD 1   ê  ;      MET_FUNC_BORNPARAM%MAXVAL+METHOD 1   %  ;      MET_FUNC_BORNPARAM%MINVAL+METHOD /   `  E   a   MET_FUNC_BORNPARAM%THIS+METHOD ,   ¥  F   a   MET_FUNC_BORNPARAM%F+METHOD ,   ë  F   a   MET_FUNC_BORNPARAM%G+METHOD     1  N   a   MET%INIT+METHOD       j       MET_INIT+METHOD %   é  E   a   MET_INIT%THIS+METHOD "   .  F   a   MET_INIT%G+METHOD '   t  t   a   MET_INIT%STATUS+METHOD $   è  8   a   MET_INIT%IER+METHOD .       ]       MPI_STATUS_SIZE+MPI_CONSTANTS %   }   X       MPI_BARRIER+MPI_BASE (   Õ   8   a   MPI_BARRIER%V0+MPI_BASE ,   !  8   a   MPI_BARRIER%IERROR+MPI_BASE -   E!  f       MPI_COMM_WORLD+MPI_CONSTANTS    «!  ç       SCH    "  @   a   SCH%CASENUMB    Ò"  @   a   SCH%RCASE    #  @   a   SCH%EQVTYPE    R#  l   a   SCH%K1    ¾#  l   a   SCH%K2    *$  l   a   SCH%K3    $  l   a   SCH%K4    %  V   a   SCH%ANGLE !   X%  g       SCH_SCHEME_ANGLE &   ¿%  E   a   SCH_SCHEME_ANGLE%THIS #   &    a   SCH_SCHEME_ANGLE%Y '    (    a   SCH_SCHEME_ANGLE%YPREV #   <*  F   a   SCH_SCHEME_ANGLE%G    *  V   a   SCH%CROSS !   Ø*  r       SCH_SCHEME_CROSS &   J+  E   a   SCH_SCHEME_CROSS%THIS #   +    a   SCH_SCHEME_CROSS%Y '   .    a   SCH_SCHEME_CROSS%YPREV '   §0    a   SCH_SCHEME_CROSS%YNEXT #   33  F   a   SCH_SCHEME_CROSS%G    y3  T   a   SCH%LAX    Í3  g       SCH_SCHEME_LAX $   44  E   a   SCH_SCHEME_LAX%THIS !   y4    a   SCH_SCHEME_LAX%Y %   7    a   SCH_SCHEME_LAX%YPREV !   9  F   a   SCH_SCHEME_LAX%G    ×9  W   a   SCH%LAXVEN "   .:         SCH_SCHEME_LAXVEN &   É:  8      SCH_SCHEME_LAXVEN%MOD '   ;  E   a   SCH_SCHEME_LAXVEN%THIS $   F;  F   a   SCH_SCHEME_LAXVEN%F (   ;  F   a   SCH_SCHEME_LAXVEN%FPREV (   Ò;  F   a   SCH_SCHEME_LAXVEN%FNEXT $   <  F   a   SCH_SCHEME_LAXVEN%G $   ^<  E   a   SCH_SCHEME_LAXVEN%M $   £<  8   a   SCH_SCHEME_LAXVEN%T    Û<  T   a   SCH%RUNGEKUTTA    /=  n       SCH_RUNGEKUTTA $   =  E   a   SCH_RUNGEKUTTA%THIS !   â=  F   a   SCH_RUNGEKUTTA%F %   (>  F   a   SCH_RUNGEKUTTA%FPREV !   n>  F   a   SCH_RUNGEKUTTA%G !   ´>  E   a   SCH_RUNGEKUTTA%M    ù>  P   a   SCH%FRUNGE    I?  c       SCH_FRUNGE     ¬?  E   a   SCH_FRUNGE%THIS    ñ?    a   SCH_FRUNGE%K    }B    a   SCH_FRUNGE%Y    	E  F   a   SCH_FRUNGE%G    OE  N   a   SCH%INIT    E  {       SCH_INIT    F  E   a   SCH_INIT%THIS "   ]F  8   a   SCH_INIT%CASENUMB !   F  8   a   SCH_INIT%EQVTYPE    ÍF  8   a   SCH_INIT%RCASE    G  F   a   SCH_INIT%G    KG  P   a   SCH%DEINIT    G  h       SCH_DEINIT %   H  >      SCH_DEINIT%ALLOCATED     AH  E   a   SCH_DEINIT%THIS 