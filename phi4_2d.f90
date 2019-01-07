!----------------------------------------------------------------------------------------------------------------
!
!  Action:   pure phi_4
!    
!----------------------------------------------------------------------------------------------------------------     
	PROGRAM MAIN 
        USE MUMBERS      
                                         
	IMPLICIT NONE
	   include 'mpif.h'      
!----------------------------Static/Helper_Variables-----------------------------------------------------------------------
        integer, parameter             :: d=2          		                                  ! dimension
        integer                        :: L                                                     ! Linear system size
	double precision               :: a,b,t                 		 	    
        double precision               :: integral,uplim,lowlim                       
        double precision               :: sitefactors(0:200),RF2(0:200),RF4(0:200),RF6(0:200),RF8(0:200),RF12(0:200)                  
        double precision               :: step,scan_step,s_step,p_step,t_step,i_t,i_p,i_s               
        integer                        :: charge,j_1,prnt_number         
        character*256                  :: results_name,cnf_name,stats_name,fname,error_file,fef       
!------------------------Configuration_Variables-----------------------------------------------------------------------
        integer, allocatable           :: site_charges(:,:),bonds_x(:,:),bonds_y(:,:)                              ! arrays        
        integer                        :: x_m1,y_m1,x_m2,y_m2,m1_charge,m2_charge,shifter,wiggler(0:1)                        ! coordinates ...
        integer                        :: m_charge,y_j,x_j
!------------------------Statistical_Variables--------------------------------------------------------------------------
         double precision, allocatable  ::  Gr(:,:), Gr_5(:,:),Gr_3(:,:),Gr_7(:,:)
         integer                        :: y_dist,x_dist
         double precision               :: Z_p,Z,R2,d_dt,Q_d_dt,Q_R2,Q_dR2_dt,Q_1,Q_2,Q_3,Q_4,Q_5,Q_6,Q_7,Q_8,Q_9,G3_l,G5_l,G7_l
         double precision               :: amax, tmax, amin, tmin, Q_10,Q_11,Q_12
         integer                        :: stats_mode,config_mode,maxcharge
         integer, parameter             :: prt_m=100000                                          
         double precision               :: Q_1_prntouts(prt_m),Q_2_prntouts(prt_m),Q_3_prntouts(prt_m),Q_8_prntouts(prt_m),Q_9_prntouts(prt_m)
         double precision               :: Q_4_prntouts(prt_m),Q_5_prntouts(prt_m),Q_6_prntouts(prt_m),Q_7_prntouts(prt_m),Q_10_prntouts(prt_m)
         double precision               :: Q_11_prntouts(prt_m),Q_12_prntouts(prt_m)
         
!------------------INITIALIZE MPI------------------------------------------------------------
        integer :: rank, comm_size, ierr, outproc, status(MPI_STATUS_SIZE)

        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, comm_size, ierr )        
!------------INPUTS from PAR file-----------------------------------------------------------------
         CALL READ_INPUT_FILE         
!-------------------Filenames for MPI----------------------------------------------------------         
         if (L>99) then
           if (rank>9) then
             write(results_name,'(A,i3,A,i2,A)') "output/",L,"_rank",rank,".dat" 
           else
             write(results_name,'(A,i3,A,i1,A)') "output/",L,"_rank",rank,".dat"
           end if
         else if (L > 9) then
            if (rank>9) then
         	write(results_name,'(A,i2,A,i2,A)') "output/",L,"_rank",rank,".dat" 
            else
            	write(results_name,'(A,i2,A,i1,A)') "output/",L,"_rank",rank,".dat"
            end if
         else
          if (rank>9) then
             write(results_name,'(A,i1,A,i2,A)') "output/",L,"_rank",rank,".dat"  
          else
             write(results_name,'(A,i1,A,i1,A)') "output/",L,"_rank",rank,".dat"
          end if
         end if                                                 
 
        if(L>99) then                               !  formatting this will break if numprocs > 99
            if(rank>9) then                 
                write(cnf_name,1343) rank,'_L', L
                write(stats_name,3343) rank,'_L', L
               write(fname,4343) rank,'_L', L
            else                
                write(cnf_name,1344) rank,'_L', L
                write(stats_name,3344) rank,'_L', L
                write(fname,4344) rank,'_L', L
            end if
        else
            if(rank>9) then                 
               write(cnf_name,1346) rank,'_L', L
                write(stats_name,3346) rank,'_L', L
                write(fname,4346) rank,'_L', L
            else                
                write(cnf_name,1347) rank,'_L', L
                write(stats_name,3347) rank,'_L', L
                write(fname,4347) rank,'_L', L
            end if                
        endif
        
        1343    format('config/cnf_',i2,A,i3)
        3343    format('config/stats_',i2,A,i3)
        4343    format('config/rand_state_',i2,A,i3)                 
        1344    format('config/cnf_',i1,A,i3)
        3344    format('config/stats_',i1,A,i3)
        4344    format('config/rand_state_',i1,A,i3)                
        1346    format('config/cnf_',i2,A,i2)
        3346    format('config/stats_',i2,A,i2)
        4346    format('config/rand_state_',i2,A,i2)        
        1347    format('config/cnf_',i1,A,i2)
        3347    format('config/stats_',i1,A,i2)
        4347    format('config/rand_state_',i1,A,i2)               
!--------------Scan Parameter (for different MPI_processes)--------------------------------------------
        t = t    + scan_step*rank                        
!------------calculate site factor integrals-------------------------------------------------------------------------------            
	
        DO charge = 0, 110                                                
        	uplim = 5.d0; lowlim = -uplim        	                          	                 
                CALL simpson(f,lowlim,uplim,integral,1000)
                sitefactors(charge) = integral               
        END DO
         
        DO  j_1 = 0, 90,2   
             RF2(j_1) = sitefactors(j_1+2)/sitefactors(j_1)
             RF4(j_1) = sitefactors(j_1+4)/sitefactors(j_1)
             RF6(j_1) = sitefactors(j_1+6)/sitefactors(j_1)
             RF8(j_1) = sitefactors(j_1+8)/sitefactors(j_1)
             RF12(j_1) = sitefactors(j_1+12)/sitefactors(j_1)
        END DO
                                                                                                    
      if (config_mode == 1) then      
         CALL READ_CONFIG         
      else      
         CALL INIT_CNF         
      end if
      
      if (stats_mode == 1) then      
         CALL READ_STATS         
      else          
         CALL INIT_STAT    
      end if
       
      if (config_mode == 0 .AND. stats_mode == 0) then
	 DO  	                                ! annealling
	     if ( rndm() < 0.5d0 ) then		    
			shifter    = 0                                 
			wiggler(0) = y_m1
			wiggler(1) = x_m1
			
			m_charge   = m1_charge
	     else		    
			shifter = 1                                  
			wiggler(0) = y_m2
			wiggler(1) = x_m2
			
			m_charge   = m2_charge	

	     end if	    	         	        
	            CALL SHIFT         	  				
		        i_t = i_t + 1.d0		               			                   
			    if (i_t >= t_step) then		    
			        exit			        
			    end if	                                                    		
          END DO		     			    		    				                       		    	             
       end if    
!------------------------------main loop---------------------------------------------------------------------                                                                                                                                              
	    DO  
	    
	     if ( rndm() < 0.5d0 ) then		    
			shifter    = 0                                 
			wiggler(0) = y_m1
			wiggler(1) = x_m1
			
			m_charge   = m1_charge

	     else		    
			shifter = 1                                  
			wiggler(0) = y_m2
			wiggler(1) = x_m2
			
			m_charge   = m2_charge	

	     end if
	    	                 
	            CALL SHIFT         	  
                    CALL MEASURE
	            step = step + 1.d0	          
	            i_p = i_p + 1.d0
	            i_s = i_s + 1.d0	               
                                   
              if (i_p >= p_step) then                            !=> write the results, including errors                   
                  i_p = 0.d0                                    
                  CALL PRNT_RESULTS                 
              end if
              
              if (i_s  >= s_step)  then                          !=> save the configuration and statistics            
                  i_s = 0.d0                
                  CALL SAVE_CONFIG 
                  CALL SAVE_STATS                     
              end if  
                         		
        END DO
               
        call MPI_FINALIZE(ierr)
                       
    CONTAINS
!---------------------------------------------------------------------------------------------
        double precision FUNCTION f(x)                 
           double precision, INTENT(IN) :: x

             f = x**(charge)*exp(-a*x**2 - (b/2.d0)*x**4)
                   
        END FUNCTION f
!------------------------NUMERICAL-INTEGRATION------------------------------------------------
 	SUBROUTINE SIMPSON(f,t1,t2,integral,n)
        
    double precision           :: integral, t1, t2	
	double precision           :: f 
	double precision           :: h, y ,s                            
	integer n, i
	! if n is odd we add +1 to make it even
	IF((n/2)*2.NE.n) n=n+1
	! loop over n (number of intervals)
	s = 0.0
	h = (t2-t1)/dfloat(n)
	DO i=2, n-2, 2
		y   = t1+dfloat(i)*h
		s = s + 2.0*f(y) + 4.0*f(y+h)
	END DO
	integral = (s + f(t1) + f(t2) + 4.0*f(t1+h))*h/3.0	
	
	END SUBROUTINE SIMPSON
!---------------------------------------------------------------------------------------------	
	    SUBROUTINE READ_INPUT_FILE 
	    
        OPEN(1, FILE='par')
        
        READ(1,*) L        
	    READ(1,*) a
	    READ(1,*) b
        READ(1,*) t           
        READ(1,*) scan_step
        READ(1,*) t_step                       ! thermolization
        READ(1,*) p_step                       ! print
        READ(1,*) s_step                       ! save
        READ(1,*) config_mode
        READ(1,*) stats_mode
        
        ALLOCATE(bonds_x(0:L-1,0:L-1))
        ALLOCATE(bonds_y(0:L-1,0:L-1))
        
        ALLOCATE(site_charges(0:L-1,0:L-1))

        ALLOCATE(Gr(0:L-1,0:L-1))  
        ALLOCATE(Gr_5(0:L-1,0:L-1))
        ALLOCATE(Gr_3(0:L-1,0:L-1))
        ALLOCATE(Gr_7(0:L-1,0:L-1))
        
        t_step = t_step*d*dble(L)*dble(L)*dble(L)
        p_step = p_step*d*dble(L)*dble(L)*dble(L)
        s_step = s_step*d*dble(L)*dble(L)*dble(L)
            
        END SUBROUTINE READ_INPUT_FILE

!---------------------------------------------------------------------------------
! algorithm: randomly choose to shift either ira or masha
 
	     SUBROUTINE SHIFT
         
        integer :: direction,draw_or_erase,target_x,target_y,target_z,bond_number,h3,h4        		                                                                
        double precision :: acc_r
        
         direction = RN(2*d)                               ! select a direction at random
         
         if (direction == 3) then
             direction = -1
         
         else if (direction == 4) then
                 direction = -2
         end if
         
         draw_or_erase = RN(2)
         
         if (draw_or_erase == 2) then
             draw_or_erase = -1
         end if
          
        if (direction == 1) then                         ! PBC and bond number lookup
            target_y = wiggler(0) + 1
            target_x = wiggler(1) 

            if (target_y == L) then
                target_y = 0
            end if			
            bond_number = bonds_y(wiggler(0),wiggler(1))
            
        else if (direction == 2) then
            target_y = wiggler(0)
            target_x = wiggler(1) + 1
         
            if (target_x == L) then
                target_x = 0
            end if
            bond_number = bonds_x(wiggler(0),wiggler(1))
            
        
        else if (direction == -1) then
            target_y = wiggler(0) - 1
            target_x = wiggler(1)

            if (target_y == -1) then
                target_y = L-1
            end if
            bond_number = bonds_y(target_y,wiggler(1))
        
        else if (direction == -2) then
            target_y = wiggler(0)
            target_x = wiggler(1) -1

            if (target_x == -1) then
                 target_x = L-1           
            end if 
            bond_number = bonds_x(wiggler(0),target_x)
        

        
        end if
        
!---------ACCEPTANCE RATIOS--------------------------------------------------------------------------------
        
        if (draw_or_erase == 1) then        
                acc_r = t/(1.d0 + bond_number)*RF2(site_charges(target_y,target_x))       
        else if (draw_or_erase == -1) then        
                acc_r = (bond_number*1.d0/t)*1.d0/RF2((site_charges(wiggler(0),wiggler(1))-2))                     
        end if
                
!------------DECISION-TO-UPDATE/UPDATE------------------------------------------------------------------

        if (rndm() < acc_r) then
            
            h3 = site_charges(wiggler(0),wiggler(1))
            h4 = site_charges(target_y,target_x)
            d_dt = d_dt - h3 - h4
                               
            if (direction == 1) then
                bonds_y(wiggler(0),wiggler(1)) = bonds_y(wiggler(0),wiggler(1)) + draw_or_erase
            else if (direction == 2) then
                bonds_x(wiggler(0),wiggler(1)) = bonds_x(wiggler(0),wiggler(1)) + draw_or_erase
            
            else if (direction == -1) then
                bonds_y(target_y,wiggler(1))   = bonds_y(target_y,wiggler(1)) + draw_or_erase
            else if (direction == -2) then
                bonds_x(wiggler(0),target_x)   = bonds_x(wiggler(0),target_x) + draw_or_erase
            
            end if 
            
            if (draw_or_erase == 1) then               
                site_charges(target_y,target_x) = site_charges(target_y,target_x) + 2                            
            else if (draw_or_erase == -1) then            
                site_charges(wiggler(0),wiggler(1)) = site_charges(wiggler(0),wiggler(1)) - 2                          
            end if
            
                       
            if (shifter == 0) then
	        y_m1 = target_y
                x_m1 = target_x
                
            else if (shifter == 1) then
	        y_m2 = target_y
		x_m2 = target_x
		
	    end if
		    
		h3 = site_charges(wiggler(0),wiggler(1))
                h4 = site_charges(target_y,target_x)    
		d_dt = d_dt + h3 + h4	
		    
          end if
        
        END SUBROUTINE SHIFT
!---------------------------------------------------------------------------------------------------------------------
        SUBROUTINE INIT_STAT       
        ! initialize statistical quantities
        
        Z   = 0.d0 
        Z_p = 0.d0
        Gr   = 0.d0
        Gr_3 = 0.d0
        G3_l = 0.d0
        Gr_5 = 0.d0
        G5_l = 0.d0
        G7_l = 0.d0
        Gr_7 = 0.d0
        R2  = 0.d0
        d_dt = sum(site_charges)   
        maxcharge = maxval(site_charges)
        Q_R2  = 0.d0
        Q_d_dt= 0.d0
        Q_dR2_dt= 0.d0

        Q_1 = 0.d0
        Q_2 = 0.d0
        Q_3 = 0.d0
        Q_4 = 0.d0
        Q_5 = 0.d0
        Q_6 = 0.d0
        Q_7 = 0.d0
        Q_8 = 0.d0
        Q_9 = 0.d0
       Q_10 = 0.d0
       Q_11 = 0.d0
       
        step = 0.d0
        i_p  = 0.d0
        i_s  = 0.d0
        i_t  = 0.d0
        
        prnt_number = 0

        Q_1_prntouts = 0.d0
        Q_2_prntouts = 0.d0
        Q_3_prntouts = 0.d0  
        Q_4_prntouts = 0.d0 
        Q_5_prntouts = 0.d0 
        Q_6_prntouts = 0.d0
        Q_7_prntouts = 0.d0
        Q_8_prntouts = 0.d0
        Q_10_prntouts = 0.d0
        Q_11_prntouts = 0.d0
              
        END SUBROUTINE INIT_STAT
!-------------------------------------------------------------------------
!		Initializing a new configuration
!-------------------------------------------------------------------------
	    SUBROUTINE INIT_CNF
	    
        integer :: ij,kl,jj
        
!       Initializes configuration, sets all the configuration variables to 0  
                                          
         bonds_x = 0
         bonds_y = 0  
  
         site_charges = 0
        
         m1_charge = 1
         m2_charge = 1       
                  
         x_m1 = RN(L)                                               
         x_m1 = x_m1 -1
         y_m1 = RN(L)
         y_m1 = y_m1 -1

                  
         x_m2 = x_m1
         y_m2 = y_m1

         
         site_charges(y_m1,x_m1) =  2
                                                 
!       Initializing random number generator  
                                                                           
         ij = 1802 + 18*rank
         kl = 9373 + 17*rank
         CALL SRANMAR(ij,kl)     
        
        END SUBROUTINE INIT_CNF
!-------------------------------------------------------------------------------------------------------------
        SUBROUTINE MEASURE
                                                        
           y_dist = abs(y_m1 - y_m2)      
           x_dist = abs(x_m1 - x_m2)    
     
                      
           Gr(y_dist,x_dist) = Gr(y_dist,x_dist) + 1.d0  
                                              
           if ( x_m1 == x_m2 .AND. y_m1 == y_m2 ) then                          	   
           	   !Z = Z + 1.d0/RF2(site_charges(y_m1,x_m1,z_m1)-2)                                             !estimator for the partition function
           	   Gr_7(y_dist,x_dist) = Gr_7(y_dist,x_dist) + RF12(site_charges(y_m1,x_m1))
           	   Gr_5(y_dist,x_dist) = Gr_5(y_dist,x_dist) + RF8(site_charges(y_m1,x_m1))
           	   Gr_3(y_dist,x_dist) = Gr_3(y_dist,x_dist) + RF4(site_charges(y_m1,x_m1))
           else 	     
          	   Gr_5(y_dist,x_dist) = Gr_5(y_dist,x_dist) + RF4(site_charges(y_m1,x_m1))*RF4(site_charges(y_m2,x_m2))
          	   Gr_3(y_dist,x_dist) = Gr_3(y_dist,x_dist) + RF2(site_charges(y_m1,x_m1))*RF2(site_charges(y_m2,x_m2))
          	   Gr_7(y_dist,x_dist) = Gr_7(y_dist,x_dist) + RF6(site_charges(y_m1,x_m1))*RF6(site_charges(y_m2,x_m2))
           end if
           
           if (y_dist > L/2) then                                       !this code block must be placed so as not to interfere with full space GGn's
               	   y_dist = L - y_dist
               	   R2 = ((y_dist)*(y_dist))*(12.d0/(L*L))
           else
           	   R2 = ((y_dist)*(y_dist))*(12.d0/(L*L)) 
           end if
           
           Z_p = Z_p + 1.d0
           Q_R2     =     Q_R2 + R2                              
           Q_d_dt   =     Q_d_dt  + d_dt           
           Q_dR2_dt =     Q_dR2_dt + d_dt*R2    
                                                                         
	    END SUBROUTINE MEASURE	    
!-------------------OUTPUT------------------------------------------
        SUBROUTINE PRNT_RESULTS
                
                integer          :: i,j,k
		double precision :: result_set(0:6),result_set2(0:5),error_set_5(0:7),error_set_6(0:7),error_set_7(0:9)                                                            
		double precision :: error_set_1(0:7),error_set_2(0:7),error_set_3(0:7),error_set_4(0:7),error_set_9(0:7) 
		double precision :: error_set_8(0:7),error_set_10(0:7),error_set_11(0:7),error_set_12(0:7)
		
           prnt_number = prnt_number + 1
		
	   if (maxval(site_charges) > maxcharge) then
           	   maxcharge = maxval(site_charges)
           end if
           
           Q_5 = 0.d0
           Q_2 = 0.d0
           Q_10 = 0.d0
           
           G3_l = 0.d0
           G5_l = 0.d0
           G7_l = 0.d0
           
                
                DO i=0,L-1
                DO j=0,L-1
                
                    Q_5 = Q_5 + Gr_5(i,j)
                    Q_2 = Q_2 + Gr_3(i,j)
                    Q_10 = Q_10 + Gr_7(i,j)
                    if (Gr(i,j) >= 0.1) then
                    	    G3_l = G3_l + (Gr_3(i,j)/Gr(i,j))
                    	    G5_l = G5_l + (Gr_5(i,j)/Gr(i,j))
                    	    G7_l = G7_l + (Gr_7(i,j)/Gr(i,j))
                    end if
               
                END DO
                END DO
                
            Q_5 = Q_5/(dble(L)*dble(L)*Gr_5(0,0))                   ! estimator for <GG_5>
            Q_2 = Q_2/(dble(L)*dble(L)*Gr_3(0,0))                   ! estimator for <GG_3>          
            Q_1 = (Q_dR2_dt)/Z_p - (Q_d_dt/Z_p)*(Q_R2/Z_p)                    ! estimator for <dGyr/dt>           
            Q_3 = Q_R2/Z_p                                                    ! estimator for <Gyr>          
            Q_4 = Z_p/(Gr(0,0)*dble(L)*dble(L))                               ! estimator for <GG> 
            Q_6 = Q_5/Q_4                                                      ! estimator for <GG5/GG> global
            Q_7 = Q_2/Q_4                                                       ! estimator for  <GG3/GG> global
            Q_8 = G3_l*(Gr(0,0)/Gr_3(0,0))*(1/(dble(L)*dble(L)))        ! estimator for locG3
            Q_9 = G5_l*(Gr(0,0)/Gr_5(0,0))*(1/(dble(L)*dble(L)))       ! estimator for locG5
            Q_10 = Q_10/(dble(L)*dble(L)*Gr_7(0,0))                    ! estimator for <GG7>
            Q_11 = G7_l*(Gr(0,0)/Gr_7(0,0))*(1/(dble(L)*dble(L)))      ! estimator for locG7
            Q_12 = Q_10/Q_4                                                        ! estimator for <GG7/GG> global
      
        		
        Q_1_prntouts(prnt_number) = Q_1                                    
        call STAT(Q_1_prntouts,prnt_number,amax,tmax,amin,tmin)                          
        error_set_1(0) = prnt_number
        error_set_1(1) = Q_1_prntouts(prnt_number)
        error_set_1(2) = amax-amin                                        
        error_set_1(3) = amax                                             
        error_set_1(4) = tmax
        error_set_1(5) = amin
        error_set_1(6) = tmin
       	error_set_1(7) = t	
       	
        Q_2_prntouts(prnt_number) = Q_2                                    ! These assignments need to be ordered sequentially
        call STAT(Q_2_prntouts,prnt_number,amax,tmax,amin,tmin)            ! because amax/tmax etc. are global        
        error_set_2(0) = prnt_number
        error_set_2(1) = Q_2_prntouts(prnt_number)                        
        error_set_2(2) = amax-amin
        error_set_2(3) = amax
        error_set_2(4) = tmax
        error_set_2(5) = amin
        error_set_2(6) = tmin
       	error_set_2(7) = t	
       	
        Q_3_prntouts(prnt_number) = Q_3                                       
        call STAT(Q_3_prntouts,prnt_number,amax,tmax,amin,tmin)                     
        error_set_3(0) = prnt_number
        error_set_3(1) = Q_3_prntouts(prnt_number)                          
        error_set_3(2) = amax-amin
        error_set_3(3) = amax
        error_set_3(4) = tmax
        error_set_3(5) = amin
        error_set_3(6) = tmin
       	error_set_3(7) = t
       	
       	Q_4_prntouts(prnt_number) = Q_4                                       
        call STAT(Q_4_prntouts,prnt_number,amax,tmax,amin,tmin)                   
        error_set_4(0) = prnt_number
        error_set_4(1) = Q_4_prntouts(prnt_number)                          
        error_set_4(2) = amax-amin
        error_set_4(3) = amax
        error_set_4(4) = tmax
        error_set_4(5) = amin
        error_set_4(6) = tmin
       	error_set_4(7) = t
       	
        Q_5_prntouts(prnt_number) = Q_5                                       
        call STAT(Q_5_prntouts,prnt_number,amax,tmax,amin,tmin)                      
        error_set_5(0) = prnt_number
        error_set_5(1) = Q_5_prntouts(prnt_number)                          
        error_set_5(2) = amax-amin
        error_set_5(3) = amax
        error_set_5(4) = tmax
        error_set_5(5) = amin
        error_set_5(6) = tmin
       	error_set_5(7) = t
       	
       	Q_7_prntouts(prnt_number) = Q_7                                       
        call STAT(Q_7_prntouts,prnt_number,amax,tmax,amin,tmin)                      
        error_set_7(0) = prnt_number
        error_set_7(1) = Q_7_prntouts(prnt_number)                          
        error_set_7(2) = amax-amin
        error_set_7(3) = amax
        error_set_7(4) = tmax
        error_set_7(5) = amin
        error_set_7(6) = tmin
       	error_set_7(7) = t
       	
       	Q_6_prntouts(prnt_number) = Q_6                                       
        call STAT(Q_6_prntouts,prnt_number,amax,tmax,amin,tmin)                      
        error_set_6(0) = prnt_number
        error_set_6(1) = Q_6_prntouts(prnt_number)                          
        error_set_6(2) = amax-amin
        error_set_6(3) = amax
        error_set_6(4) = tmax
        error_set_6(5) = amin
        error_set_6(6) = tmin
       	error_set_6(7) = t

       	Q_8_prntouts(prnt_number) = Q_8                                       
        call STAT(Q_8_prntouts,prnt_number,amax,tmax,amin,tmin)                      
        error_set_8(0) = prnt_number
        error_set_8(1) = Q_8_prntouts(prnt_number)                          
        error_set_8(2) = amax-amin
        error_set_8(3) = amax
        error_set_8(4) = tmax
        error_set_8(5) = amin
        error_set_8(6) = tmin
       	error_set_8(7) = t
       	
       	Q_9_prntouts(prnt_number) = Q_9                                       
        call STAT(Q_9_prntouts,prnt_number,amax,tmax,amin,tmin)                      
        error_set_9(0) = prnt_number
        error_set_9(1) = Q_9_prntouts(prnt_number)                          
        error_set_9(2) = amax-amin
        error_set_9(3) = amax
        error_set_9(4) = tmax
        error_set_9(5) = amin
        error_set_9(6) = tmin
       	error_set_9(7) = t
       	
       	Q_10_prntouts(prnt_number) = Q_10                                       
        call STAT(Q_10_prntouts,prnt_number,amax,tmax,amin,tmin)                      
        error_set_10(0) = prnt_number
        error_set_10(1) = Q_10_prntouts(prnt_number)                          
        error_set_10(2) = amax-amin
        error_set_10(3) = amax
        error_set_10(4) = tmax
        error_set_10(5) = amin
        error_set_10(6) = tmin
       	error_set_10(7) = t
       	
       	Q_11_prntouts(prnt_number) = Q_11                                       
        call STAT(Q_11_prntouts,prnt_number,amax,tmax,amin,tmin)                      
        error_set_11(0) = prnt_number
        error_set_11(1) = Q_11_prntouts(prnt_number)                          
        error_set_11(2) = amax-amin
        error_set_11(3) = amax
        error_set_11(4) = tmax
        error_set_11(5) = amin
        error_set_11(6) = tmin
       	error_set_11(7) = t
       	
       	Q_12_prntouts(prnt_number) = Q_12                                       
        call STAT(Q_12_prntouts,prnt_number,amax,tmax,amin,tmin)                      
        error_set_12(0) = prnt_number
        error_set_12(1) = Q_12_prntouts(prnt_number)                          
        error_set_12(2) = amax-amin
        error_set_12(3) = amax
        error_set_12(4) = tmax
        error_set_12(5) = amin
        error_set_12(6) = tmin
       	error_set_12(7) = t
       
               open(53+rank,file = results_name)
       		
		result_set(0) = t
		result_set(1) = Q_3
		result_set(2) = Q_1
		result_set(3) = Q_4       
		result_set(4) = Q_2
		result_set(5) = Q_5
		result_set(6) = Q_10
		
		result_set2(0) = Q_7
		result_set2(1) = Q_6		
		result_set2(2) = Q_12
		result_set2(3) = Q_8
		result_set2(4) = Q_9
		result_set2(5) = Q_11
          
             write(53+rank,*) ' '
             write(53+rank,'(A,i1,A,i3,A,F6.4,A,F6.4,A,ES10.3E2)')'rank = ',rank, ' size L =',L ,' a = ', a,'  b = ',b,'  mc run steps', step
             write(53+rank,*) ' '
             write(53+rank,*) ' t    <Gyr/G_0>  <dGyr/dt>  <GG>  <GG3>  <GG5>  <GG7>  '
             write(53+rank,*) ' '             
             write(53+rank,'(F6.4,F8.5,F9.2,F9.5,F8.5,F8.5,F8.5)') result_set 
             write(53+rank,*) ' '  
             write(53+rank,*) ' <GG3/GG>  <GG5/GG>  <GG7/GG>  locG3  locG5  locG7     '
             write(53+rank,*) ' '             
             write(53+rank,'(F8.5,F8.5,F8.5,F10.5,F10.5,F10.5)') result_set2 
             write(53+rank,*) ' '
             write(53+rank,'(A,i12)') 'number of printouts: ', prnt_number
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F12.5)') ' <Gyr/Gyr_0> = ', error_set_3(1), '  max-min = ', error_set_3(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_3(3),' at -> ',error_set_3(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_3(5),' at -> ',error_set_3(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F12.5)') ' <dGyr/dt> = ', error_set_1(1), '  max-min = ', error_set_1(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_1(3),' at -> ',error_set_1(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_1(5),' at -> ',error_set_1(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F11.7)') ' <GG1> = ', error_set_4(1), '  max-min = ', error_set_4(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_4(3),' at -> ',error_set_4(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_4(5),' at -> ',error_set_4(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F11.7)') ' <GG3> = ', error_set_2(1), '  max-min = ', error_set_2(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_2(3),' at -> ',error_set_2(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_2(5),' at -> ',error_set_2(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F11.7)') ' <GG5> = ', error_set_5(1), '  max-min = ', error_set_5(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_5(3),' at -> ',error_set_5(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_5(5),' at -> ',error_set_5(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F11.7)') ' <GG7> = ', error_set_10(1), '  max-min = ', error_set_10(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_10(3),' at -> ',error_set_10(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_10(5),' at -> ',error_set_10(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F11.7)') ' <GG3>/<GG> = ', error_set_7(1), '  max-min = ', error_set_7(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_7(3),' at -> ',error_set_7(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_7(5),' at -> ',error_set_7(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F11.7)') ' <GG5>/<GG> = ', error_set_6(1), '  max-min = ', error_set_6(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_6(3),' at -> ',error_set_6(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_6(5),' at -> ',error_set_6(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F11.7)') ' <GG7>/<GG> = ', error_set_12(1), '  max-min = ', error_set_12(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_12(3),' at -> ',error_set_12(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_12(5),' at -> ',error_set_12(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F11.7)') ' locG3 = ', error_set_8(1), '  max-min = ', error_set_8(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_8(3),' at -> ',error_set_8(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_8(5),' at -> ',error_set_8(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F11.7)') ' locG5 = ', error_set_9(1), '  max-min = ', error_set_9(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_9(3),' at -> ',error_set_9(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_9(5),' at -> ',error_set_9(6)
             write(53+rank,*) ' '
             write(53+rank,'(A,F12.5,A,F11.7)') ' locG7 = ', error_set_11(1), '  max-min = ', error_set_11(2)
             write(53+rank,'(F12.5,A,F12.5)') error_set_11(3),' at -> ',error_set_11(4)
             write(53+rank,'(F12.5,A,F12.5)') error_set_11(5),' at -> ',error_set_11(6)
             write(53+rank,*) ' '

             	close(53+rank) 
          
	    END SUBROUTINE PRNT_RESULTS
!-----------------stat-------------------------------------------------------------------------------------------
      SUBROUTINE STAT(prnt_storage,n,amax,tmax,amin,tmin) 
!     Analyzing 3/4 print-out
	      
      double precision :: prnt_storage, amax, tmax, amin, tmin, aa
      integer n, i
      dimension prnt_storage(n)

      amax=-1.d200; amin=1.d200

      DO i=n/4+1, n;  aa=prnt_storage(i)
         if (aa > amax) then
            amax=aa; tmax=i
         end if
         if (aa < amin) then
            amin=aa; tmin=i
         end if
      END DO

      tmax=tmax/n; tmin=tmin/n
      END SUBROUTINE STAT	    
!------------SAVE CONFIGURATION--------------------------------------------------------------------
       
       SUBROUTINE SAVE_CONFIG
                  
       open(14, file=cnf_name)  
                               
       write(14,*) y_m1,x_m1    
       write(14,*) y_m2,x_m2
       write(14,*) m1_charge,m2_charge     
       write(14,*) bonds_y       
       write(14,*) bonds_x 
  
       write(14,*) site_charges
       
       close(14)                                                  
                                                     	                   						   
	   call save_rndm(fname)                     
	                                              ! (this uses file reference number 23)      
       END SUBROUTINE SAVE_CONFIG       
!--------------------Save Stats -----------------------------------------------------------------------------  
     SUBROUTINE SAVE_STATS
     
     open(15, file=stats_name)
     
     write(15,*) step  
     write(15,*) i_s
     write(15,*) i_p 

     write(15,*) Z
     write(15,*) Z_p
     write(15,*) Gr  
     write(15,*) Gr_3
     write(15,*) Gr_5
     write(15,*) Gr_7
     write(15,*) d_dt
     write(15,*) Q_R2
     write(15,*) Q_d_dt
     write(15,*) Q_dR2_dt
     write(15,*) G3_l
     write(15,*) G5_l
     write(15,*) G7_l
     write(15,*) Q_1
     write(15,*) Q_2
     write(15,*) Q_3
     write(15,*) Q_4
     write(15,*) Q_5
     write(15,*) Q_6
     write(15,*) Q_7
     write(15,*) Q_8
     write(15,*) Q_9
     write(15,*) Q_10
     write(15,*) Q_11
     write(15,*) Q_12
     write(15,*) Q_1_prntouts
     write(15,*) Q_2_prntouts
     write(15,*) Q_3_prntouts
     write(15,*) Q_4_prntouts
     write(15,*) Q_5_prntouts
     write(15,*) Q_6_prntouts
     write(15,*) Q_7_prntouts
     write(15,*) Q_8_prntouts
     write(15,*) Q_9_prntouts
     write(15,*) Q_10_prntouts
     write(15,*) Q_11_prntouts
     write(15,*) Q_12_prntouts
     write(15,*) prnt_number
     write(15,*) maxcharge
     
     close(15)
         
     END SUBROUTINE SAVE_STATS
!----------------------------------------------------------------------------------------------------------------------
       SUBROUTINE READ_CONFIG
       
       open(14, file=cnf_name)  
                               
       read(14,*) y_m1,x_m1   
       read(14,*) y_m2,x_m2
       read(14,*) m1_charge,m2_charge      
       read(14,*) bonds_y       
       read(14,*) bonds_x  

       read(14,*) site_charges
       
       close(14)
       
       call read_rndm(fname)
                    
       END SUBROUTINE READ_CONFIG
!----------------------------------------------------------------------------------------------------------------------
       SUBROUTINE READ_STATS
       
     open(15, file=stats_name)
     
     read(15,*) step  
     read(15,*) i_s
     read(15,*) i_p  
     
     read(15,*) Z
     read(15,*) Z_p
     read(15,*) Gr 
     read(15,*) Gr_3
     read(15,*) Gr_5
     read(15,*) Gr_7
     read(15,*) d_dt
     
     read(15,*) Q_R2 
     read(15,*) Q_d_dt
     read(15,*) Q_dR2_dt
     read(15,*) G3_l
     read(15,*) G5_l
     read(15,*) G7_l
     read(15,*) Q_1
     read(15,*) Q_2
     read(15,*) Q_3
     read(15,*) Q_4
     read(15,*) Q_5
     read(15,*) Q_6
     read(15,*) Q_7 
     read(15,*) Q_8
     read(15,*) Q_9
     read(15,*) Q_10
     read(15,*) Q_11
     read(15,*) Q_12
     read(15,*) Q_1_prntouts
     read(15,*) Q_2_prntouts
     read(15,*) Q_3_prntouts
     read(15,*) Q_4_prntouts
     read(15,*) Q_5_prntouts
     read(15,*) Q_6_prntouts
     read(15,*) Q_7_prntouts
     read(15,*) Q_8_prntouts
     read(15,*) Q_9_prntouts
     read(15,*) Q_10_prntouts
     read(15,*) Q_11_prntouts
     read(15,*) Q_12_prntouts
     read(15,*) prnt_number
     read(15,*) maxcharge
     
     close(15)
              
       END SUBROUTINE READ_STATS
!----------------------------------------------------------------------------------------------------------------------                   	
      END PROGRAM MAIN
      