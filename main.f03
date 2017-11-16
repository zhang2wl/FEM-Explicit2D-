!**********************************************************************
!  Main program
!
!  Purpose: Explicit FEM calculation 
!  Function:   1. Reads in nodes, element, BC info
!              2. Initialize matrices, U,V,A,R_ext,R_int
!              3. Loop over time
!                 3.1 Loop over element
!                    3.1.1 Get internal force vector
!                 3.2 Calculate the new acceleration
!                 3.3 Calculate the new velocity
!                 3.4 Calculate the new displacement
!              4. Write out data
!  Author:Wenlong Zhang
!  Last edit: 11/16/2017
!     
      program main
         use explicit !defines U,V,A. assembly R_int, calculates A
         use extract  !calculate stresses
         use BC       !applies BC condition
         use vtk_io   !module to write output
         implicit none
         real(kind=8) :: time1,time2,time3,time4,time5,time6, time7, &
                         time8,time9,time_start, time_end !variables use
                         ! for time measurement
         real(kind=8),parameter::zero=0.d0
         integer::i,num_Gauss_points,ii,j 
         real(kind=8),dimension(8,3)::temp1 !temp1=trans(B)*mat_mtx
         call cpu_time(time_start)
         call cpu_time(time1)
!         
!-- Read element, node, BC, force, crack data from file
!
         call determine_size(node_nums, element_nums, BC_nums, &
                           force_nums, crack_nums)
         call build_connectivity(node_nums,node_data,element_nums, &
         element_data, BC_nums, BC_data,force_nums,force_data, &
         crack_nums, crack_data)
         call change_shape(BC_data,BC_nums, force_data,force_nums,&
                  new_BC_data,new_force_data)
         call get_critical_dt(element_nums, element_data,node_data,&
                  velocity,dt_cri)
         call create_mass_matrix(element_data,node_data,node_nums,&
               element_nums,rho,mass)
         NT = get_total_steps(T,dt_cri)
         print*, 'Total simulation time is', T
         print*, 'Critical time step is ', dt_cri
         print*, 'There are total', NT, ' steps'
         Gauss_4 = Gauss_integration_points()
         
        ! print '("Time for building connectivity:", f10.5)',time2-time1
         num_Gauss_points=4
         mat_mtx = material_matrix(E,miu)
         call cpu_time(time2)
         call initialize(U,V,A,U_dt,V_dt,A_dt,R_int,node_nums)
         call initialize_A_temp(A_temp,node_nums)
         call initialize_external_force(R_ext,node_nums)
         call initialize_stress_strain(sig,eps,node_nums)
         call cpu_time(time3)

        ! print '("Time cost for Initialization:", f10.5)',time2-time1
         do it = 1,NT 
            call cpu_time(time4)
            R_int = zero 
            call get_external_force(new_force_data,force_nums,node_nums,&
                     R_ext,it) 
            if (it==1) then
               V_dt%values = V%values + 0.5*CFL*dt_cri*A%values
               U_dt%values = U%values + dt_cri*CFL*V_dt%values
               A_dt%values = 0.d0
            else
               V%values = V_dt%values
               U%values = U_dt%values
               A%values = A_dt%values

               V_dt%values = V%values + dt_cri*CFL*A%values
               U_dt%values = U%values + dt_cri*CFL*V_dt%values
               call apply_BC(new_BC_data,BC_nums,new_force_data,force_nums,&
                   U_dt,it,NT)
            endif

         call cpu_time(time6)
      !--- loop over element ---
         do i=1,element_nums
             call element_info(i,element_data,node_data,en,node)
             call get_deformed_node_coord(node,en,U_dt)
             U_el(1:7:2,1) = U_dt%values(en(1,1:4),1)
             U_el(2:8:2,1) = U_dt%values(en(1,1:4),2)
             K_el(:,:) = zero 
          !--- loop over integration point ---   
             do intp=1,num_Gauss_points
               gauss(1,1:2) = Gauss_4(intp,1:2)
               call shape_function(gauss,N)
               call Jacobian(gauss,node,jaco,det_jacobian)
               N_ab = get_N_ab(gauss)
               call get_B(jaco,N_ab,B,inv_J)
               B_trans = transpose(B)      
               call &
               dgemm('n','n',8,3,3,1.d0,B_trans,8,mat_mtx,3,0.d0,temp1,8)
               call &
               dgemm('n','n',8,8,3,1.d0,temp1,8,B,3,0.d0,K_el_int,8)
               K_el_int = K_el_int*det_jacobian
               K_el = K_el + K_el_int
               call get_strain(eps,U_el,B,en,intp)
               call get_stress(sig,eps,mat_mtx,en,intp)
             enddo
         !--- calculate internal stress ---
             call dgemv('n',8,8,1.d0,K_el,8,U_el,1,0.d0,F_el,1)
             call cpu_time(time8)

             call assemble_internal_force(F_el,R_int,node_nums,en)
             call cpu_time(time9)

        enddo
        call cpu_time(time7)
         !print '(" Time for all element internal force calculation:", & 
         !         f10.5)',time3-time2
         call calculate_acceleration(R_ext,R_int,Mass,node_nums,A_temp)
         call extrapolate_stress_to_node(sig,node_nums)
      !--- calculate acceleration    
         A_dt%values(:,1) = A_temp(1:2*node_nums-1:2)
         A_dt%values(:,2) = A_temp(2:2*node_nums:2)
         call cpu_time(time5)
         print '("Timestep", I8, ",   cpu time=", f10.5)',it,time5-time4
         enddo
      !---- Write output to file -----
         title = 'disp_stress.vtk'
         open(unit=1,file=title, status='replace')
         call vtk_sig_write(1,title,node_nums,element_nums,&
               element_order, node_data,element_data,U_dt,sig)
         close(1)
      !-------------------------------
         call cpu_time(time_end)
         print '("Read in data time: ", f10.5)',time2-time1
         print '("Initialize time:   ", f10.5)',time3-time2
         print '("For each cycle, loop over element cost:" &
         ,f10.5)',time7-time6
         print '("For each cycle, calculate A cost:", &
         f10.5)', time5-time7 
         print '("Within each element, assembly cost:", &
         f10.5)', time9-time8
         print '("Total time:        ", f10.5)',time_end-time_start

      end program
