      module step1
      implicit none
      
      type my_type
         integer::id
         integer,dimension(:,:),allocatable::my_data
      end type my_type
         integer,dimension(:),allocatable::AA,BB,CC

         type(my_type)::A
         integer:: n !size of A
      contains

         subroutine initialize(A,n)
            implicit none
            type(my_type),intent(out)::A
            integer,intent(in)::n
            integer,parameter::zero=0
            allocate(A%my_data(n,n))
            A%id = zero
            A%my_data(:,:) = zero         
         end subroutine initialize

         subroutine initialize_AA(AA,BB,CC,n)
            integer,dimension(:),allocatable::AA,BB,CC
            integer,intent(in)::n
            allocate(AA(n))
            allocate(BB(n))
            allocate(CC(n))
         end subroutine initialize_AA

         subroutine substract(AA,BB,CC,n)
            integer,intent(in)::n
            integer,dimension(n)::AA,BB,CC
            BB = 1
            CC = 2
            AA = BB-CC
            print*,BB
         end subroutine substract
         
      end module step1


      module step2
      use step1
      implicit none

      contains
         subroutine change_A(A,n)
            implicit none
            type(my_type),intent(inout)::A
            integer,intent(in)::n
            integer,parameter::one=1
            
            A%my_data(1,1) = 1                    
         end subroutine change_A
      end module step2
