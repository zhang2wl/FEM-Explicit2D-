      module explicit
      use parameters
      use connectivity
      use element
      implicit none
         type disp_type
            integer,dimension(:),allocatable::node_id
            real(kind=8),dimension(:,:),allocatable::values
         end type disp_type

         type(disp_type)::U,V,A,U_dt,V_dt,A_dt
         real(kind=8),parameter::T  =1, &
                                 CFL=0.8
         real(kind=8)::dt
         real(kind=8)::dt_cri
         integer::NT,it
         real(kind=8),dimension(:),allocatable::R_int
         real(kind=8),dimension(:),allocatable::R_ext
         real(kind=8),dimension(:),allocatable::A_temp
      contains

         subroutine get_critical_dt(element_nums,element_data, & 
            node_data,velocity,dt_cri) 
            type(element_type),intent(in) ::element_data
            type(node_type),intent(in)    ::node_data
            real(kind=8),intent(in)       ::velocity
            real(kind=8),intent(out)      ::dt_cri
            real(kind=8),dimension(4,2)   ::node
            real(kind=8)      ::Area,length,dt
            integer,dimension(1,4)        ::en
            integer::i
            integer,intent(in)::element_nums
            dt_cri = 100.0
            ! print*, velocity, element_nums
            do i=1,element_nums
               call element_info(i,element_data,node_data,en,node)
               Area = element_area(node)
               !print*,'element',i,'has node: ',Area
               length = sqrt(Area)
               dt = length/velocity
               if (dt < dt_cri) then
                  dt_cri = dt
               !   print*,i, dt
               endif
            enddo
         end subroutine get_critical_dt
         
         function get_total_steps(T,dt_cri) result(NT)
            real(kind=8),intent(in)::T,dt_cri
            real(kind=8):: dt
            integer::NT
            dt = dt_cri*CFL
            print*,dt
            NT = floor(T/dt)
         end function get_total_steps
         
         subroutine initialize(U,V,A,U_dt,V_dt,A_dt,R_int,node_nums)
            integer,intent(in)::node_nums
            type(disp_type),intent(out)::U,V,A
            type(disp_type),intent(out)::U_dt,V_dt,A_dt
            real(kind=8),dimension(:),allocatable,intent(out)::R_int

            allocate(R_int(node_nums*2))
            R_int(:) = 0.d0
            allocate(U%node_id(node_nums))
            allocate(U%values(node_nums,2))
            U%node_id(:) = 0
            U%values(:,:) = 0.d0

            allocate(V%node_id(node_nums))
            allocate(v%values(node_nums,2))
            V%node_id(:) = 0
            V%values(:,:) = 0.d0

            allocate(A%node_id(node_nums))
            allocate(A%values(node_nums,2))
            A%node_id(:) = 0
            A%values(:,:) = 0.d0

            allocate(U_dt%node_id(node_nums))
            allocate(U_dt%values(node_nums,2))
            U_dt%node_id(:) = 0
            U_dt%values(:,:) = 0.d0

            allocate(V_dt%node_id(node_nums))
            allocate(v_dt%values(node_nums,2))
            V_dt%node_id(:) = 0
            V_dt%values(:,:) = 0.d0

            allocate(A_dt%node_id(node_nums))
            allocate(A_dt%values(node_nums,2))
            A_dt%node_id(:) = 0
            A_dt%values(:,:) = 0.d0
         end subroutine initialize
         
         subroutine initialize_A_temp(A_temp,node_nums)
            real(kind=8),dimension(:),allocatable,intent(out)::A_temp
            integer,intent(in)::node_nums
            allocate(A_temp(2*node_nums))
            A_temp(:) = 0.d0
         end subroutine initialize_A_temp

         subroutine initialize_external_force(R_ext,node_nums)
            real(kind=8),dimension(:),allocatable,intent(out)::R_ext
            integer,intent(in)::node_nums
            allocate(R_ext(2*node_nums))
         end subroutine

         subroutine assemble_internal_force(F_el,R_int,node_nums,en)
            real(kind=8),dimension(8,1),intent(in)::F_el
            integer,dimension(1,4),intent(in)::en
            integer,intent(in)::node_nums
            integer,dimension(4):: indices
            real(kind=8),dimension(node_nums*2),intent(inout):: R_int

            indices = en(1,1:4)
            R_int(2*indices-1) = R_int(2*indices-1) + F_el(1:7:2,1)
            R_int(2*indices) = R_int(2*indices) + F_el(2:8:2,1)
         end subroutine assemble_internal_force
  
         subroutine get_external_force(force_data,force_nums,node_nums,&
                           R_ext,it)
            integer,intent(in)::node_nums,force_nums,it
            integer,dimension(force_nums),intent(in)::force_data
            real(kind=8),dimension(2*node_nums),intent(out)::R_ext
            real(kind=8),parameter::force_dt=0.d0,zero=0.d0
          !  print*, force_data 
            R_ext(:) = zero
            R_ext(2*force_data)=force_dt*it
         end subroutine get_external_force


         subroutine &
         calculate_acceleration(R_ext,R_int,Mass,node_nums,A_temp)
            integer,intent(in)::node_nums
            real(kind=8),dimension(node_nums*2),intent(in)::R_ext
            real(kind=8),dimension(node_nums*2),intent(in)::R_int
            real(kind=8),dimension(node_nums*2),intent(in)::Mass
            real(kind=8),dimension(node_nums*2),intent(out)::A_temp
            integer::i

!         do ii=1,node_nums
!            print*,R_ext(ii)            
!         enddo
            A_temp= R_ext-R_int
            !A_temp = 0.d0
!            print*,'--------------------------------'
!         do ii=1,node_nums
!            print*,R_ext(ii)            
!         enddo
            do i=1,2*node_nums
               A_temp(i) = A_temp(i)/Mass(i)
            enddo
         end subroutine calculate_acceleration
      end module explicit


