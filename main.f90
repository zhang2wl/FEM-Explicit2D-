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
    use mpi
    use connectivity
    use explicit !defines U,V,A. assembly R_int, calculates A
    use extract  !calculate stresses
    use BC       !applies BC condition
    use vtk_io   !module to write output
    use partition 
    implicit none
    real(kind=8) :: time1,time2,time3,time4,time5,time6, time7, &
                    time8,time9,time_start, time_end !variables use
                    ! for time measurement
    real(kind=8),parameter::zero=0.d0
    integer::i,num_Gauss_points ,ierr,j,k,gid,request,send_count, start_id,&
            status(MPI_STATUS_SIZE),recv_request
    real(kind=8),dimension(8,3)::temp1 !temp1=trans(B)*mat_mtx
    type(element_type)::element_data_local 
    type(node_type)::node_data_local
    type(disp_type)::U_local, U_dt_local,V_local,V_dt_local,A_local,&
         A_dt_local
   real(kind=8),allocatable::R_int_local(:), buffer(:), temp(:)
   integer,parameter:: big_number = 10000
   integer:: recv_buffer(big_number), recv_id_buffer(big_number), &
            local_id(big_number),send_r_int(big_number)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Initiate MPI
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
         call MPI_init(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs,ierr)
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
!
!-- Get the shared node between processors and local elements
!
call get_shared_node(element_data,element_nums,node_nums, &
        shared_node_local,shared_node,shared_node_size, &
        send_nums,num_local, epart,npart, numprocs,myid)
call get_local_element_id(epart, element_nums, numprocs,&
          local_elements,local_nums)

!call get_local_node_id(npart,node_nums, numprocs,local_nodes,&
!         local_node_nums, shared_node,shared_node_size)
!if (myid .eq. 0 ) then
!   print*, 'shared node size is:', shared_node_size
!   print*, local_node_nums
!endif
!        print*, 'Get shared node'
!      print*, 'alread get local node id'
!        print*, send_nums(1,1:2)
!!
!!-- Initialize arrays
!!
Gauss_4 = Gauss_integration_points()
num_Gauss_points=4
mat_mtx = material_matrix(E,miu)
call initialize(U,V,A,U_dt,V_dt,A_dt,R_int,node_nums)
call initialize_A_temp(A_temp,node_nums)
call initialize_external_force(R_ext,node_nums)
call initialize_stress_strain(sig,eps,node_nums)
print*,'initialized UVA'
!
!-- Create local arrays
!
!allocate(element_data_global%element_id(local_nums(myid+1)))
!allocate(element_data_global%element_id(local_nums(myid+1)))
allocate(element_data_local%element_id(local_nums(myid+1)))
allocate(element_data_local%vertex(local_nums(myid+1),4))
element_data_local%element_id(1:local_nums(myid+1)) = (/(i,i=1,local_nums(myid+1))/)
element_data_local%vertex(1:local_nums(myid+1),:) = &
   element_data%vertex(local_elements(:,myid+1),1:4)
print*, 'initialized local element'
!element_data_global%element_id = element_data_local%element_id
!element_data_global%vertex = element_data_global%vertex

!allocate(node_data_local%node_id(local_node_nums(myid+1)))
!allocate(node_data_local%coord(local_node_nums(myid+1),22))
!node_data_local%node_id(1:local_node_nums(myid+1)) = &
!      (/(i,i=1,local_node_nums(myid+1))/) ! This is global id
!node_data_local%coord(1:local_nums(myid+1),1:2) = &
!   node_data%coord(local_nodes(:,myid+1),1:2)
!call expand_shared_node(shared_node_local,num_local,node_data_local,&
!         new_shared_node_local)

! Change element_data_local's vertex values to local node id
!do i=1,local_node_nums(myid+1)
!   gid = node_data_local%node_id(i)
!   do j=1,local_nums(myid+1)
!      do k=1,4
!         if (element_data_local%vertex(j,k)==gid) then
!            element_data_local%vertex(j,k) = i
!         endif
!      enddo
!   enddo
!enddo
! Localize U_dt,U,V_dt,V,A_dt,A
!allocate(U_local%values(local_node_nums(myid+1),2))
!allocate(U_dt_local%values(local_node_nums(myid+1),2))
!allocate(V_local%values(local_node_nums(myid+1),2))
!allocate(V_dt_local%values(local_node_nums(myid+1),2))
!allocate(A_local%values(local_node_nums(myid+1),2))
!allocate(A_dt_local%values(local_node_nums(myid+1),2))
!!
!allocate(U_local%node_id(local_node_nums(myid+1)))
!allocate(U_dt_local%node_id(local_node_nums(myid+1)))
!allocate(V_local%node_id(local_node_nums(myid+1)))
!allocate(V_dt_local%node_id(local_node_nums(myid+1)))
!allocate(A_local%node_id(local_node_nums(myid+1)))
!allocate(A_dt_local%node_id(local_node_nums(myid+1)))

!allocate(R_int_local(2*local_node_nums(myid+1)))
!!
!!-- Loop over time
!!
allocate(buffer(2*node_nums))
do it = 1,NT 
   R_int= zero
   call get_external_force(new_force_data,force_nums,node_nums, R_ext,it) 
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
      !print*, U_dt%values(2202:2203,2)
   endif
   call cpu_time(time6)
   do i=1,local_nums(myid+1)
      ! This output en will be local node id
       call element_info(i,element_data_local,node_data, en,node)

       call get_deformed_node_coord(node,en,U_dt)
       U_el(1:7:2,1) = U_dt%values(en(1,1:4),1)
       U_el(2:8:2,1) = U_dt%values(en(1,1:4),2)
       K_el(:,:) = zero 
   !=============================================================
   !== Then following calculation are totally localized:
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
   !==============================================================
         call cpu_time(time9)
  enddo
  !-- Send and receive shared node info -- 
  !
  !  print*, size(R_int)
  if (myid .ne. 0) then
     call mpi_send(R_int, 2*node_nums,MPI_double, 0,0,mpi_comm_world,&
            ierr)
     !call mpi_wait(request, status, ierr)
  else
     do i=1,(numprocs-1)
        call mpi_recv(buffer, 2*node_nums, MPI_double, MPI_any_source,0,&
                      mpi_comm_world, status, ierr)
      !  call mpi_wait(recv_request, status,ierr)
        R_int = R_int + buffer
     enddo
  endif
!---------------------------------------
  call cpu_time(time7)
   !print '(" Time for all element internal force calculation:", & 
   !         f10.5)',time3-time2
   if (myid .eq. 0) then
      call calculate_acceleration(R_ext,R_int,Mass,node_nums,A_temp)
      !call extrapolate_stress_to_node(sig,node_nums)
      print*, sig%values(2202,2)
   !--- calculate acceleration    
      A_dt%values(:,1) = A_temp(1:2*node_nums-1:2)
      A_dt%values(:,2) = A_temp(2:2*node_nums:2)
      call cpu_time(time5)
      print '("Timestep", I8, ",   cpu time=", f10.5)',it,time7-time6

      call MPI_Bcast(A_temp,2*node_nums,MPI_double, 0,MPI_comm_world,ierr)
   !print*, A_temp(4400:4404) 
   else
   ! manager has all the acceleration data
   ! now send  the updated accelration data to other processors 

      call MPI_Bcast(A_temp,2*node_nums,MPI_double, 0,MPI_comm_world,ierr)
      A_dt%values(:,1) = A_temp(1:2*node_nums-1:2)
      A_dt%values(:,2) = A_temp(2:2*node_nums:2)
   endif
enddo
 
!  
!---- Write output to file -----
if (myid .eq.  0) then
   title = 'disp_stress.vtk'
   open(unit=10001,file=title,Access='sequential', status='replace',&
      position='append',iostat = ierr)
   if (ierr /= 0 ) then
      print*, 'Error in opening in APPEND mode existing file; code error=',&
      ierr
      STOP '---EXECUTION ABORTED---'
   endif
   call vtk_essential_write(10001,title,node_nums,element_nums,&
         element_order, node_data,element_data)
   call vtk_scalar_write(10001, node_nums, npart, 'npart')
   allocate(temp(node_nums))
   temp  = U%values(:,2)
   call vtk_scalar_write_double(10001,node_nums, temp,'disp-v')
!   call vtk_vector_write(10001, node_nums,2, temp, 'disp')
   close(10001)
   call cpu_time(time_end)
!  print '("Read in data time: ", f10.5)',time2-time1
!   print '("Initialize time:   ", f10.5)',time3-time2
!   print '("For each cycle, loop over element cost:" &
!   ,f10.5)',time7-time6
!   print '("For each cycle, calculate A cost:", &
!   f10.5)', time5-time7 
!   print '("Within each element, assembly cost:", &
!   f10.5)', time9-time8

   print '("Total time:        ", f10.5)',time_end-time_start
   endif      
call MPI_Finalize(ierr)
end program main
