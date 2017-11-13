      module extract
      use element,only: U_el,B 
      use connectivity, only: node_nums,en
      use parameters,only: mat_mtx
      implicit none
      type stress_type
         integer,dimension(:),allocatable:: node_id, add_time
         real(kind=8),dimension(:,:),allocatable::values
      end type stress_type
      type(stress_type)::sig,eps
      integer::intp

      contains
         subroutine initialize_stress_strain(sig,eps,node_nums)
            type(stress_type),intent(out)::sig,eps
            integer,intent(in):: node_nums
            integer::i

            allocate(sig%node_id(node_nums))
            allocate(sig%add_time(node_nums))
            allocate(sig%values(node_nums,3))
            do i=1,node_nums
               sig%node_id = i
            enddo
            sig%add_time = 0
            sig%values(:,:) = 0.d0

            allocate(eps%node_id(node_nums))
            allocate(eps%add_time(node_nums))            
            allocate(eps%values(node_nums,3))
            do i=1,node_nums
               eps%node_id = i
            enddo
            eps%values(:,:) = 0.d0
            eps%add_time = 0
         end subroutine initialize_stress_strain
            
         subroutine get_strain(eps,U_el,B,en,intp)
            real(kind=8),dimension(8,1),intent(in)::U_el
            real(kind=8),dimension(3,8),intent(in)::B
            integer,intent(in)::intp
            integer,dimension(1,4),intent(in)::en            
            type(stress_type),intent(inout)::eps
            real(kind=8),dimension(3)::eps_intp

            call dgemv('n',3,8,1.d0,B,3,U_el,1,0.d0,eps_intp,1)

            eps%values(en(1,intp),1) = eps_intp(1)
            eps%values(en(1,intp),2) = eps_intp(1)
            eps%values(en(1,intp),3) = eps_intp(1)

         end subroutine get_strain

         subroutine get_stress(sig,eps,mat_mtx,en,intp)
            type(stress_type),intent(inout)::sig
            type(stress_type),intent(in)::eps
            real(kind=8),dimension(3,3),intent(in)::mat_mtx
            integer,dimension(1,4),intent(in)::en
            integer,intent(in)::intp
            real(kind=8),dimension(3)::eps_buffer,sig_intp

            eps_buffer = eps%values(en(1,intp),:) 
            call &
            dgemv('n',3,3,1.d0,mat_mtx,3,eps_buffer,1,0.d0,sig_intp,1)
            sig%values(en(1,intp),:) = sig%values(en(1,intp),:) + &
                                       sig_intp(:)
            sig%add_time(en(1,intp)) = sig%add_time(en(1,intp)) + 1
         end subroutine get_stress
         
         subroutine extrapolate_stress_to_node(sig,node_nums)
            type(stress_type), intent(inout)::sig
            integer,intent(in)::node_nums
            integer::i
            do i=1,node_nums
               sig%values(i,:) = sig%values(i,:)/sig%add_time(i)
            enddo
         end subroutine extrapolate_stress_to_node


      end module extract
