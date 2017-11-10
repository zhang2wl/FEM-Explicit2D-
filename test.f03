      program test
      use step1
      use step2
         implicit none
         integer,dimension(2,2)::B
         integer::ii
         n = 2
         call initialize(A,n)
         call change_A(A,n) 

         B(:,:) = 3
         A%my_data(:,2) = B(:,2) 
         print*, A%my_data
         print*, (1==1)
         print*, (1/=1)
         call try_exit()
         call initialize_AA(AA,BB,CC,n)
         call substract(AA,BB,CC,n)
      contains
      subroutine try_exit()
         do ii = 1,10
            print*,ii
            if (ii==5) exit
         enddo
         print*,'I am outside the loop'
      end subroutine try_exit
      end program test
      
          
      
