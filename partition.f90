! This module is to find the shared nodes at partition interface.
! Partition is the sub-domain divided using Metis for Parallel computing
! Need to use the connectivity module to provide the element connecitivity info

! Need to read in the npart and epart file. Build the partition info
! I already have the number of elements and number of nodes num_nodes
! Data needed: 
!  element_nums               
!  node_nums
!  epart file
!  npart file
!  element_data

module partition
use connectivity,only:element_type,element_data,element_nums, node_nums,&
               node_type
                     
implicit none
integer,allocatable::epart(:), npart(:), local_elements(:,:),local_nums(:),&
local_nodes(:,:),local_node_nums(:), new_shared_node_local(:,:)
integer,dimension(:,:), allocatable:: shared_node, shared_node_local,send_nums
integer     :: numprocs, shared_node_size,myid,num_local
contains

! First need to read in npart and epart file
subroutine read_epart_data(file_id,file_name,array_size, epart)
   integer,intent(in)               :: file_id, array_size
   character(len=*),intent(in)      :: file_name
   integer                          :: i, check, ierr
   integer,allocatable,intent(out)  :: epart(:)
   allocate(epart(array_size))
   open(file_id,file=file_name, status="old",action="read",iostat=ierr)
!   print*, "able to open file"
   do i=1,array_size
      read(file_id, *, iostat=check) epart(i)
      if (check>0) then
         print*, 'Error reading data from partition file'
      endif
   enddo
!   print*, "Finish reading"
   close(file_id)
end subroutine read_epart_data      

! Then loop over element to build the shared node array
subroutine get_shared_node(element_data, element_nums, &
           node_nums, shared_node_local, shared_node,shared_node_size, &
           send_nums,num_local,epart,npart,numprocs,myid)
   type(element_type),intent(in)    :: element_data
   integer,intent(in)               :: element_nums, node_nums,numprocs,myid
   integer                          :: i,nid, num_shn,j,k,proc_temp, num_procs
   integer                          :: shared_node_temp(node_nums,3)
   integer,allocatable              :: shared_node(:,:) , &
                                       send_nums_temp(:,:)
   integer,allocatable,intent(out)  :: epart(:),npart(:),send_nums(:,:)
   integer:: shared_node_size, num_local
   integer,allocatable:: shared_node_size_each_proc(:)
   integer,allocatable:: shared_node_local(:,:)
   call read_epart_data(1,'element.mesh.epart',element_nums,epart)
   call read_epart_data(2,'element.mesh.npart',node_nums,npart)
!   print*, "Able to read file"
   num_shn = 1
   shared_node_temp(:,:) = 0

   do i=1,size(element_data%element_id)
      do j = 1,size(element_data%vertex,2)
         nid = element_data%vertex(i,j)
         if (npart(nid) .ne. epart(i)) then ! This node is shared node
            shared_node_temp(num_shn,1) = num_shn
            shared_node_temp(num_shn,2) = npart(nid)
            shared_node_temp(num_shn,3) = epart(i)
            num_shn = num_shn + 1
         endif
      enddo
   enddo
   shared_node_size = 0
   do i=1,node_nums
      if (shared_node_temp(i,1) .eq. 0) then
         ! It means the none zero numbers in shared_node_temp is i-1
         shared_node_size = i-1
         exit         
      endif  
   enddo
   allocate(shared_node(shared_node_size,3))
   do i=1,shared_node_size
      shared_node(i,1:3) = shared_node_temp(i,1:3)
   enddo
   !
   !-- Now I need to sort shared_node so it has the following format
   !     node_id     proc #1     proc #2
   !                 0           
   !                 0           
   !                 .
   !                 1
   !                 1
   !                 .
   !                 .
   !---------------------------------------------------------------  
   call quicksortmatrix(shared_node,1,shared_node_size,2, shared_node_size,3)
   ! Need to find on each processor, how many nodes are there
   allocate(shared_node_size_each_proc(numprocs))
   shared_node_size_each_proc = 0
   do i=1,shared_node_size
      shared_node_size_each_proc(shared_node(i,2)+1) = &
         shared_node_size_each_proc(shared_node(i,2)+1) + 1
   enddo

   i=1   
   do 
      if (shared_node(i,2) .ne. myid) then
         i=i+1
      else
         exit
      endif
   enddo
   num_local = shared_node_size_each_proc(myid+1)
   ! Shared_node_local defines the shared node only in that processor
   !  node id,  shread_proc_id (in ascending order)
   allocate(shared_node_local(num_local,2))
   shared_node_local(1:num_local,1) = shared_node(i:(i+num_local),1)
   shared_node_local(1:num_local,2) = shared_node(i:(i+num_local),3)

   call quicksortmatrix(shared_node_local,1,num_local,2,num_local,2) 
!   print*, 'processor:', myid
!   do i=1,num_local
!      print*, shared_node_local(i,:)
!   enddo
   !--- Need to get how many nodes to send to each processor
   ! Fist find out how many adjacent processors
   proc_temp = shared_node_local(1,2)
   num_procs = 1
   do i=1,num_local
      if (shared_node_local(i,2) == proc_temp) then
      ! do nothing   
      else
         proc_temp = shared_node_local(i,2) ! update proc_temp
         num_procs = num_procs + 1
      endif
   enddo

   ! Then allocate send_nums
   !                 (proc_id, total # of nodes)
   allocate(send_nums_temp(num_procs,2))
   send_nums_temp(:,:) = 0
   ! At last find out how many nodes to send to each processor
   !
   proc_temp = shared_node_local(1,2)
   k = 1
   do i=1,num_local
      if (shared_node_local(i,2) == proc_temp) then
         send_nums_temp(k,2) = send_nums_temp(k,2) + 1
         send_nums_temp(k,1) = proc_temp
      else !now it is a different processor
         proc_temp = shared_node_local(i,2) ! update proc_temp
         k = k + 1
         send_nums_temp(k,2) = send_nums_temp(k,2) + 1
      endif
   enddo
   ! need to remove the zero lines of send_nums
   allocate(send_nums(k,3))
   send_nums(:,1:2) = send_nums_temp(1:k,1:2)
   ! The third row calculates the sum of previous rows of column 2, then plus 1
   send_nums(:,3) = 0 !Initialize the third column
   send_nums(1,3) = 1 !The first row starts from 1
   do i=2,k
      j = i
      do
         send_nums(i,3) = send_nums(i,3) + send_nums(j-1,2)
         j = j-1
         if (j == 1) then
            exit
         endif
      enddo   
   enddo
!   print*, 'processor', myid
!   do i=1,k
!      print*, send_nums(i,:)
!   enddo
end subroutine get_shared_node 

subroutine expand_shared_node(shared_node_local, num_local, &
      node_data_local,new_shared_node_local)
   integer,intent(in)::num_local
   integer,intent(in)::shared_node_local(num_local,2)
   integer::i,j,gid
   type(node_type),intent(in)::node_data_local
   integer,allocatable,intent(out)::new_shared_node_local(:,:)
   allocate(new_shared_node_local(num_local,3))

   do i=1,num_local
   do j=1,size(node_data_local%node_id)
      gid = node_data_local%node_id(j)
      if (shared_node_local(i,1) == gid) then
         new_shared_node_local(i,1) = j
      endif
   enddo
   enddo
   new_shared_node_local(:,2:3) = shared_node_local(:,:)
end subroutine expand_shared_node

recursive subroutine quicksortmatrix(a, first, last, col,row_max,col_max)
  implicit none
  integer,intent(inout):: a(row_max,col_max)
  integer x(1,col_max), t(1,col_max)
  integer first, last
  integer i, j
  integer col ! Based on this column the quick sort is carried out
  integer row_max, col_max ! number of rows and columns in the input matrix a
  x(1,:) = a( (first+last) / 2 , :)
  i = first
  j = last

  do
     do while (a(i,col) < x(1,col))
        i=i+1
     end do
     do while (x(1,col) < a(j,col))
        j=j-1
     end do
     if (i >= j) exit
     t(1,:) = a(i,:);  a(i,:) = a(j,:);  a(j,:) = t(1,:)
     i=i+1
     j=j-1
  end do

  if (first < i-1) call quicksortmatrix(a, first, i-1,col,row_max,col_max)
  if (j+1 < last)  call quicksortmatrix(a, j+1, last,col,row_max,col_max)
end subroutine quicksortmatrix

subroutine get_local_element_id(epart,element_nums,&
            numprocs,local_elements,local_nums)
   ! This subroutine is to generate a matrix whose rows stores the element
   ! number and the columns represent the id of processors. So to find 
   ! which element that processor has, I can go to myid-th column and read
   ! the data
   ! output: 
   !  local_nums: the numbers of elements each processor has
   !  local_elements: the element global id for each processor
   integer,intent(in)   :: epart(element_nums), element_nums,numprocs
   integer,allocatable,intent(out) :: local_elements(:,:),local_nums(:)
   integer                       :: i
   allocate(local_elements(element_nums,numprocs))
   allocate(local_nums(numprocs))
   local_nums(:) = 1
   local_elements(:,:) = 0
   do i=1, element_nums
      local_elements(local_nums(epart(i)+1),epart(i)+1) = i
      local_nums(epart(i)+1) = local_nums(epart(i)+1) + 1
!      print*, local_nums(epart(i)+1)
   enddo
!   print*, local_nums

end subroutine get_local_element_id 

subroutine get_local_node_id(npart,node_nums, numprocs,local_nodes,&
      local_node_nums, shared_node,shared_node_size)
   ! This subroutine is to generate a matrix whose rows stores the element
   ! number and the columns represent the id of processors. So to find 
   ! which element that processor has, I can go to myid-th column and read
   ! the data
   ! output: 
   !  local_node_nums: the numbers of elements each processor has
   !  local_nodes: the node global id for each processor
   integer,intent(in)   :: npart(node_nums), node_nums,numprocs
   integer,allocatable,intent(out) :: local_nodes(:,:),local_node_nums(:)
   integer                       :: i,node1,proc,j,check
   integer, intent(in)::shared_node(shared_node_size,3), shared_node_size
   allocate(local_nodes(node_nums,numprocs))
   allocate(local_node_nums(numprocs))
   local_node_nums(:) = 1
   local_nodes(:,:) = 0
   do i=1, node_nums
      local_nodes(local_node_nums(npart(i)+1),npart(i)+1) = i
      local_node_nums(npart(i)+1) = local_node_nums(npart(i)+1) + 1
   enddo
   ! Need to consider the shared_node as well
   do i=1,shared_node_size
      node1 = shared_node(i,1)
      proc = shared_node(i,2)
      check = 0
      do j=1,local_node_nums(proc+1)! loop over the number of nodes that
        if (node1 .eq. local_nodes(j,proc+1)) then
           check=1
           exit
        endif 
      enddo
      
      if (check .eq. 0) then
         ! if this shared_node is not included yet
         local_nodes(local_node_nums(proc+1)+1,proc+1) = node1
         local_node_nums(proc+1) = local_node_nums(proc+1)+1
      endif
   enddo
end subroutine get_local_node_id 


end module partition
