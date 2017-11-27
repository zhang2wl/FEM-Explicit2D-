      module BC
      use connectivity
      use explicit
      implicit none  
      
      contains
         subroutine apply_BC(BC_data,BC_nums,force_data,force_nums, &
                             U_dt,it,NT)
            type(disp_type),intent(inout)::U_dt
            integer,intent(in)::BC_nums,force_nums
            integer,dimension(BC_nums),intent(in)::BC_data
            integer,dimension(force_nums),intent(in)::force_data
            integer,intent(in)::it,NT
            real(kind=8),parameter:: zero=0.d0
            real(kind=8),parameter:: disp = 0.01d0
            
            U_dt%values(BC_data,1:2) = zero

            U_dt%values(force_data,1) =  disp/10/NT*it
            U_dt%values(force_data,2) =  disp/NT*it
!            print*, disp/NT*it
!            print*, U_dt%values(force_data,2)
         end subroutine
      end module BC
