      program main
         implicit none
         integer::ierror,i,j
         real,dimension(121,2)::A
         character(len=80)::B = "abc"
         open(1,file="node_heat.txt",status="old",&
         action="read",iostat=ierror)
         do i=1,121
            read(1,*) (A(i,j),j=1,2)
            print*, A(i,1),A(i,2)
         enddo
         
         B = "abcde"
         print*, B         
      end program main
