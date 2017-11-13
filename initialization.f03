      module initialization
      use connectivity,only: node_nums
      implicit none
         type disp_type
            integer,dimension(:),allocatable::node_id
            real(kind=8),dimension(:,:),allocatable::values
         end type disp_type

         type(disp_type)::U,V,A,U_dt,V_dt,A_dt
      contains
 
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
      end module initialization
