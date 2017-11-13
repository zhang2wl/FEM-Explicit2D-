      program main
         use explicit
         use extract
         use BC
         implicit none
         real(kind=8) :: time_start,time_end
         real(kind=8),parameter::zero=0.d0
         integer::i,num_Gauss_points,ii,j
         real(kind=8),dimension(8,3)::temp1
         call cpu_time(time_start)
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
!      print*, BC_nums,force_nums
         call  get_critical_dt(element_nums, element_data,node_data,&
                  velocity,dt_cri)
         call create_mass_matrix(element_data,node_data,node_nums,&
               element_nums,rho,mass)
         NT = get_total_steps(T,dt_cri)
         print*, 'Total simulation time is', T
         print*, 'Critical time step is ', dt_cri
         print*, 'There are total', NT, ' steps'
         Gauss_4 = Gauss_integration_points()
         
         num_Gauss_points=4
         mat_mtx = material_matrix(E,miu)

         call initialize(U,V,A,U_dt,V_dt,A_dt,R_int,node_nums)
         call initialize_A_temp(A_temp,node_nums)
         call initialize_external_force(R_ext,node_nums)
         call initialize_stress_strain(sig,eps,node_nums)
         do it = 1,NT 
         !   print*, 'this is cycle',it
            R_int = zero 
           ! print*, new_force_data,force_nums
            call get_external_force(new_force_data,force_nums,node_nums,&
                     R_ext,it) 
            !print*, R_ext
            !        print*, 'Get the external force'
            if (it==1) then
               V_dt%values = V%values + 0.5*CFL*dt_cri*A%values
               U_dt%values = U%values + dt_cri*CFL*V_dt%values
               !A_dt%values = (R_ext-R_int)./Mass
               A_dt%values = 0.d0
               !print*, A_dt%values
            else
               V%values = V_dt%values
               U%values = U_dt%values
               A%values = A_dt%values

               V_dt%values = V%values + dt_cri*CFL*A%values
               U_dt%values = U%values + dt_cri*CFL*V_dt%values
         !      print*, 'updated u values'
               call apply_BC(new_BC_data,BC_nums,new_force_data,force_nums,&
                   U_dt,it,NT)
            endif

         do i=1,element_nums
             call element_info(i,element_data,node_data,en,node)
             call get_deformed_node_coord(node,en,U_dt)
             U_el(1:7:2,1) = U_dt%values(en(1,1:4),1)
             U_el(2:8:2,1) = U_dt%values(en(1,1:4),2)
             K_el(:,:) = zero 
             
             do intp=1,num_Gauss_points
               gauss(1,1:2) = Gauss_4(intp,1:2)
               !write(*,'(2F10.4)') gauss
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
!               if(i==1)then
!                do ii=1,8
!                   write(*,'(8F10.4)') (K_el_int(ii,j),j=1,8)
!                enddo
!                do ii=1,8
!                   write(*,'(3F10.4)') (temp1(ii,j),j=1,3)
!                enddo
!               endif
               call get_strain(eps,U_el,B,en,intp)
               call get_stress(sig,eps,mat_mtx,en,intp)
             enddo
             call dgemv('n',8,8,1.d0,K_el,8,U_el,1,0.d0,F_el,1)
             if (i==2) then
!                do ii=1,8
!                   write(*,'(F8.2)') F_el
!               enddo
!                print*,'Element 1 B matrix at last integration point: '
!                do ii=1,3
!                   write(*,'(8F10.4)') (B(ii,j),j=1,8)
!                enddo
!                do ii=1,2*node_nums
!                   write(*,*) R_int(ii)
!                enddo
                !do ii=1,2
                !  write(*,'(4F10.4)') (N_ab(ii,j),j=1,4)
                !enddo
                !do ii=1,4
                !   write(*,'(2F10.4)') (node(ii,j),j=1,2)
                !enddo
                              ! print*, det_jacobian
!                print*,'Element 1 K_el at the last integration point: '
                !do ii=1,8
                !   write(*,'(8F10.4)') (K_el_int(ii,j),j=1,8)
                !enddo

!               write(*,*) det_jacobian
!               write(*,*) B_trans
!               write(*,*) U_el
             endif
             call assemble_internal_force(F_el,R_int,node_nums,en)

        enddo
         !print*,'size of R_int is',size(R_int)
         !do ii=1,node_nums
         !   print*,R_ext(ii)
         !enddo
         call calculate_acceleration(R_ext,R_int,Mass,node_nums,A_temp)
         !do ii=1,2*node_nums
         !!    write(*,*) A_temp(ii)
         !enddo
         call extrapolate_stress_to_node(sig,node_nums)
          
         A_dt%values(:,1) = A_temp(1:2*node_nums-1:2)
         A_dt%values(:,2) = A_temp(2:2*node_nums:2)
         ! print*, R_ext(1:5)
         write(*,100) it, U_dt%values(6,:)
         100 format ('Time step ', I4 ,' .Disp is ', 2F12.6)
         enddo
         call cpu_time(time_end)

         print*, 'Running time is: ', time_end-time_start
      end program
