      module element
      use connectivity
      implicit none
         real(kind=8),dimension(1,2)::gauss
         real(kind=8),dimension(4,2)::node,Gauss_4
         real(kind=8),dimension(2,2)::jaco,inv_J
         real(kind=8),dimension(1,4)::N
         real(kind=8),dimension(2,4)::N_ab,N_xy
         real(kind=8),dimension(3,8)::B
         real(kind=8),dimension(8,3)::B_trans
         real(kind=8)::det_jacobian,area
         integer,dimension(1,4)::en
         real(kind=8),dimension(:),allocatable::Mass
      contains
         function Gauss_integration_points() result(Gauss_4)
            real(kind=8),dimension(4,2):: Gauss_4      
         !   integer::i,j

            Gauss_4 = reshape((/-1,1,1,-1,-1,-1,1,1/),shape(Gauss_4))
         !   do i=1,size(Gauss_4,1)
         !      print*, (Gauss_4(i,j),j=1,size(Gauss_4,2))
         !   enddo
            Gauss_4 = 1/sqrt(3.0)*Gauss_4
         end function Gauss_integration_points
   
         subroutine shape_function(gauss,N)
            real(kind=8),dimension(1,2),intent(in)::gauss
            real(kind=8),dimension(1,4),intent(out)::N
            real(kind=8)::a,b
            real(kind=8),parameter::one=1.0d0
            a = gauss(1,1)
            b = gauss(1,2)
            N(1,1) = 0.25*(one-a)*(one-b)
            N(1,2) = 0.25*(one+a)*(one-b)
            N(1,3) = 0.25*(one+a)*(one+b)
            N(1,4) = 0.25*(one-a)*(one+b)
         end subroutine shape_function

         subroutine element_info(el_id,element_data,node_data,en,node)
            integer,intent(in)::el_id
            type(node_type),intent(in)::node_data
            type(element_type),intent(in)::element_data
            integer,dimension(1,4),intent(out)::en
            real(kind=8),dimension(4,2),intent(out)::node

            en(1,1:4)  =element_data%vertex(el_id,1:4)
            node = node_data%coord(en(1,1:4),1:2)
         end subroutine element_info

         function element_area(node) result(area)
            real(kind=8),dimension(4,2),intent(in)::node
            real(kind=8)::area
            real(kind=8),dimension(1,2)::gauss
            real(kind=8),dimension(4,2)::gauss_4
            real(kind=8)::det_jacobian
            integer::intp
            gauss_4 = Gauss_integration_points()
            
            area = 0.0
            do intp = 1,4
               gauss(1,1:2) = gauss_4(intp,1:2)
               call Jacobian(gauss,node,jaco,det_jacobian)
               !print*, det_jacobian
               area = area + det_jacobian
            enddo
         end function element_area
         
         subroutine create_mass_matrix(element_data,node_data, &
               node_nums, element_nums,rho,Mass)
            integer::el_id
            integer,intent(in)::element_nums,node_nums
            type(node_type),intent(in)::node_data
            type(element_type),intent(in)::element_data
            integer,dimension(1,4) ::en
            integer,dimension(4)::en2
            real(kind=8),dimension(4,2)::node
            real(kind=8)::area, mass_per_node
            real(kind=8),intent(in)::rho 
            real(kind=8),dimension(:),allocatable,intent(out)::Mass
            
            allocate(Mass(2*node_nums))
            Mass(:) = 0.d0
            print*,size(Mass)
            do el_id=1,element_nums
               call element_info(el_id,element_data,node_data,en,node)
               area = element_area(node)               
               mass_per_node = area*rho/4.d0
               en2 = en(1,1:4)
               Mass(2*en2-1) = Mass(2*en2-1) + mass_per_node
               Mass(2*en2) = Mass(2*en2) + mass_per_node
            enddo
         end subroutine create_mass_matrix

         subroutine Jacobian(gauss,node,jaco,det_jacobian)
            real(kind=8),dimension(1,2),intent(in)::gauss
            real(kind=8),dimension(4,2),intent(in)::node
            real(kind=8),dimension(2,2),intent(out)::jaco
            real(kind=8),intent(out)::det_jacobian
            real(kind=8),dimension(2,4)::dphi
            real(kind=8)::a,b
            real(kind=8),parameter::one=1.0d0
           ! integer::i,j,k
            a = gauss(1,1)
            b = gauss(1,2)
            dphi(1,1:4) = 0.25*(/b-one,one-b,one+b,-one-b/)
            dphi(2,1:4) = 0.25*(/a-one,-one-a,one+a,one-a/)
            jaco=0.0
            call dgemm('n','n',2,2,4,1.d0,dphi,2,node,4,0.d0,jaco,2)
!            do j=1,2
!               do i=1,2
!                  do k=1,4
!                  jaco(i,j) = jaco(i,j) + dphi(i,k)*node(k,j) 
!                  enddo
!               enddo
!            enddo
       !call dgemm('n','n',m,m,m,alpha,amat,m,bmat,m,beta,cmatblas,m)
            det_jacobian = jaco(1,1)*jaco(2,2)-jaco(1,2)*jaco(2,1)
          !  print *, det_jacobian
         end subroutine Jacobian
         
         function get_N_ab(gauss) result(dphi)
            real(kind=8),dimension(1,2),intent(in)::gauss
            real(kind=8),dimension(2,4)::dphi
            real(kind=8)::a,b
            real(kind=8),parameter::one=1.0d0
            a = gauss(1,1)
            b = gauss(1,2)
            dphi(1,1:4) = 0.25*(/b-one,one-b,one+b,-one-b/)
            dphi(2,1:4) = 0.25*(/a-one,-one-a,one+a,one-a/) 
         end function get_N_ab
         
         subroutine get_B(jacobian,N_ab,B,inv_J)
            real(kind=8),dimension(2,2),intent(in)::jacobian
            real(kind=8),dimension(2,4),intent(in)::N_ab
            real(kind=8),dimension(2,2),intent(out)::inv_J
            real(kind=8),dimension(2,4)::N_xy
            real(kind=8),dimension(3,8),intent(out)::B
            real(kind=8)::det_jacobian
            integer::i,j,k
            N_xy(:,:) = 0.0
            !-- First need to inverse Jacobian
            det_jacobian = jacobian(1,1)*jacobian(2,2)-&
            jacobian(1,2)*jacobian(2,1)
            
           ! do i=1,2
           !    write(*,'(2F16.10)') (jacobian(i,j),j=1,2)
           ! enddo
            inv_J(1,1) =  jacobian(2,2)/det_jacobian
            inv_J(2,2) =  jacobian(1,1)/det_jacobian
            inv_J(1,2) = -jacobian(1,2)/det_jacobian
            inv_J(2,1) = -jacobian(2,1)/det_jacobian
            !-- Then multiply inv_J with N_ab
            do j=1,4
               do i=1,2
                  do k=1,2
                     N_xy(i,j) = N_xy(i,j) + inv_J(i,k)*N_ab(k,j)
                  enddo
               enddo
            enddo
            B(1,1:7:2) = N_xy(1,1:4)
            B(2,2:8:2) = N_xy(2,1:4)
            B(3,1:7:2) = N_xy(2,1:4)
            B(3,2:8:2) = N_xy(1,1:4)
         end subroutine get_B
         
         
      end module element
