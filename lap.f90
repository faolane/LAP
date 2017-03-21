! Linear assignment problem (LAP) solver
! -----------------------------------------------
! The first method (method = 1) is a brute force
! solution that computes all the permutations (n!)
! using Heap's algorithm - time complexity O(n!)
! The second method (method = 2) uses the Munkres
! algorithm (also called the Hungarian algo) to
! solve the problem - time complexity O(n^3)
! -----------------------------------------------
! F. Brieuc - March 2017

module var
   ! module with global variables
   implicit none

   ! norm for matrix representation C(j,i) : j = columns, i = rows
   integer, dimension(:,:), allocatable :: C    ! cost matrix (max sum)
   integer, dimension(:,:), allocatable :: CC   ! cost matrix (min sum)
   integer :: n      ! dimension of C - assumed to be a (nxn) square matrix
   integer :: method ! 1: brute force (heap's algo) - 2: Munkre's algorithm
   integer :: mode   ! 0: minimize sum, 1: maximize sum

   ! following variables are use only for Munkres algo
   integer, dimension(:,:), allocatable :: M    ! mask matrix
   integer, dimension(:), allocatable :: rowCover, colCover !cover row and cols
   integer :: pathRow0 = 0, pathCol0 = 0  ! starting point for path finding part
end module var

program lap
   use var
   implicit none

   ! solution
   integer :: sumSol ! maximal sum
   integer, dimension(:), allocatable   :: jSol ! assocaiated "column" indices j
                                                ! for each row i
   integer :: sum, i

   ! read/generate C and initialize
   call readInputFile('lap.in')

   allocate(jSol(n))

   sumSol = 0

   select case(method)
   case (1)
      ! solve lap by computing all permutations using Heap's algorithm
      call heap(jSol)
   case(2)
      ! Munkre's algorithm
      call munkres(jSol)
   case default
      stop('error : wrong value for parameter method (1 or 2)')
   end select

   ! print the solution
   sumSol = 0
   do i = 1, n
      sumSol = sumSol + C(jSol(i),i)
   enddo

   if (mode == 1) then
      write(6, '(a29,i5)') 'The maximum possible sum is: ', sumSol
   else if (mode == 0) then
      write(6, '(a29,i5)') 'The minimum possible sum is: ', sumSol
   endif
   write(6, '(a38)') 'obtained using the following elements:'
   write(6, '(a29)') 'row i, column j, C(i,j), sum '
   sum = 0
   do i = 1, n
      sum = sum + C(jSol(i),i)
      write(6,'(i4, a4, i4, a5, i4, a1, i6)') &
      & i, '    ' ,jSol(i), '     ', C(jSol(i),i), ' ', sum
   enddo

   deallocate(jSol)
   deallocate(C)
   deallocate(CC)

end program lap


subroutine initSeed()
   ! initialize random number generator seed
   ! using the computer clock
   ! should work for any seed size
   implicit none

   integer, dimension(:), allocatable :: time
   integer :: seedSize, i, j

   ! get the size of the seed
   call random_seed(size=seedSize)
   if (seedSize > 8) then
      allocate(time(seedSize))
      call date_and_time(values = time(1:8))
      j = 2
      do i = 9, seedSize
         time(i) = time(j) + time(j - 1)
         j = j + 1
         if (j > 8) j = 2
      enddo
      call random_seed(put = time)
      deallocate(time)
   else
      allocate(time(8))
      call date_and_time(values = time)
      call random_seed(put = time(1:seedSize))
      deallocate(time)
   endif

end subroutine initSeed

subroutine readInputFile(filename)
   use var
   implicit none

   character(len = *), intent(in) :: filename
   integer :: gen  ! gen = 0 -> read C, gen = 1 random C, gen = -1 worst case C
   integer :: i, j, io, tmp
   real :: r

   ! initialize seed of random number generator
   call initSeed()

   open(10, file=trim(filename))

   read(10, *) n, gen, method, mode
   if (n <= 0) then
      write(6,*) 'Dimension of matrix C has to be greater then 0!'
      stop('error while reading input file')
   endif
   if (gen /= 0 .and. gen /= 1 .and. gen /= -1) then
      stop('error : wrong value for parameter gen (0, 1 or -1)')
   endif

   allocate(C(n,n))
   allocate(CC(n,n))

   write(6, '(a33)') '*** Linear assignment problem ***'
   write(6, '(a16, i3)') 'dimension : n = ', n
   write(6, '(a12)') 'Cost matrix:'

   if (gen == 0) then
      do i = 1, n
         read(10,*,iostat=io) C(:,i)
         if (io /= 0) stop('error while reading input file')
      enddo
   else if (gen == 1) then
      do i = 1, n
         do j = 1, n
            call random_number(r)
            C(j,i) = 1 + int(9 * r)
         enddo
      enddo
   else if (gen == -1) then
      do i = 1, n
         do j = 1, n
            C(j,i) = j * i
         enddo
      enddo
   endif
   call printMatrix(C,n)

   close(10)

   if (mode == 1) then ! maximize sum
      if (method == 1) write(6,'(a32)') 'Maximizing sum using Heap''s algo'
      if (method == 2) write(6,'(a34)') 'Maximizing sum using Munkres'' algo'
      write(6,*) ''
      tmp = huge(i)
      do i = 1, n
         do j = 1, n
            CC(j,i) = -C(j,i)
            if (CC(j,i) < tmp) tmp = CC(j,i)
         enddo
      enddo
      CC = CC - tmp
   else if (mode == 0) then ! minimize sum
      if (method == 1) write(6,'(a32)') 'Minimizing sum using Heap''s algo'
      if (method == 2) write(6,'(a34)') 'Minimizing sum using Munkres'' algo'
      write(6,*) ''
      do i = 1, n
         do j = 1, n
            CC(j,i) = C(j,i)
         enddo
      enddo
   else
      stop('error : wrong value for parameter mode (0 or 1)')
   endif

end subroutine readInputFile

subroutine printMatrix(M,d)
   implicit none

   integer, intent(in) :: d
   integer, dimension(d,d), intent(in) :: M
   integer :: i, j

   if (d == 1) then
      write(6,'(i4)') M(1,1)
   else
      do i = 1, d
         do j = 1, d - 1
            write(6, '(i4)', advance="no") M(j,i)
         enddo
         write(6, '(i4)') M(d,i)
      enddo
   endif

   write(6,*) ''

end subroutine printMatrix

subroutine heap(jSol)
   ! Compute all the possible (n!) permutations for jSol(i)
   ! using Heap's algorithm and find the smallest sum
   ! Heap, B. R., The Computer Journal 6, 3 (1963)
   ! see : https://en.wikipedia.org/wiki/Heap's_algorithm
   use var
   implicit none

   integer, dimension(n), intent(out) :: jSol ! solution indices

   integer, dimension(:), allocatable :: p  ! for permutations
   integer, dimension(:), allocatable :: j  ! choosen indice of row i
   integer :: sum, sum_tmp
   integer :: i, tmp

   allocate(j(n))
   allocate(p(n))

   ! *** Compute all the possible (n!) permutations for jj(i) ***
   ! Heap's algorithm
   sum = 0
   do i = 1, n
      jSol(i) = i
      j(i) = i
      sum = sum + CC(i,i)
      p(i) = 1
   enddo
   sum_tmp = sum

   i = 1;
   do while (i <= n)
      if (p(i) < i) then
         if (mod(i,2) == 0) then
            ! i odd - swap i <-> p(i)
            sum_tmp = sum_tmp - CC(j(p(i)),p(i)) - CC(j(i),i)
            tmp = j(p(i))
            j(p(i)) = j(i)
            j(i) = tmp
            sum_tmp = sum_tmp + CC(j(p(i)),p(i)) + CC(j(i),i)
         else
            ! i even - swap i <-> 1
            sum_tmp = sum_tmp - CC(j(1),1) - CC(j(i),i)
            tmp = j(1)
            j(1) = j(i)
            j(i) = tmp
            sum_tmp = sum_tmp + CC(j(1),1) + CC(j(i),i)
         end if
         if (sum_tmp < sum) then
            sum = sum_tmp
            jSol(:) = j(:)
         endif
         p(i) = p(i) + 1
         i = 1
      else
         p(i) = 1
         i = i + 1
      endif
   enddo

   deallocate(j)
   deallocate(p)

end subroutine heap

subroutine munkres(jSol)
   ! Implementation of the Munkres algorithm (also referred to as the Hungarian
   ! algorithm). J. Munkres, Journal of the SIAM 5, 1 (1957)
   ! The following implementation is based on
   ! http://csclab.murraystate.edu/%7Ebob.pilgrim/445/munkres.html
   use var
   implicit none

   integer, dimension(n), intent(out) :: jSol ! solution indices

   integer :: step, i, j, tmp
   logical :: done

   done = .false.
   step = 1
   tmp = 0

   allocate(M(n,n))      ! mask matrix - contains starred zeros
   allocate(rowCover(n)) ! to keep track of covered rows
   allocate(colCover(n)) ! to keep track of covered columns

   do i = 1, n
      M(:,i) = 0
      rowCover(i) = 0
      colCover(i) = 0
   enddo

   do while(.not. done)
      select case(step)
      case(1)
         call step1(step)
      case(2)
         call step2(step)
      case(3)
         call step3(step)
      case(4)
         call step4(step)
      case(5)
         call step5(step)
      case(6)
         call step6(step)
      case default ! done
         do i = 1, n
            do j = 1, n
               if (M(j,i) == 1) jSol(i) = j
            enddo
         enddo
         done = .true.
      end select
   enddo

   deallocate(M)
   deallocate(rowCover)
   deallocate(colCover)

end subroutine munkres

subroutine step1(step)
   ! row reduction : for each row find the smallest value and substract it from
   ! all elements of that row. Go to step 2.
   use var
   implicit none

   integer, intent(out) :: step

   integer :: minVal, i, j

   do i = 1, n
      minVal = CC(1,i)
      do j = 1, n
         if (CC(j,i) < minVal) minVal = CC(j,i)
      enddo
      CC(:,i) = CC(:,i) - minVal
   enddo

   step = 2

end subroutine step1

subroutine step2(step)
   ! Search for zeros.
   ! Find a zero (Z) in the matrix. If no zeros has been previously starred in
   ! its row and column then star Z. Go to step 3.
   use var
   implicit none

   integer, intent(out) :: step

   integer :: i, j

   do i = 1, n
      do j = 1, n
         if (CC(j,i) == 0 .and. rowCover(i) == 0 .and. colCover(j) == 0) then
            M(j,i) = 1
            rowCover(i) = 1
            colCover(j) = 1
         endif
      enddo
   enddo
   ! uncovers
   do i = 1, n
      rowCover(i) = 0
      colCover(i) = 0
   enddo

   step = 3

end subroutine step2

subroutine step3(step)
   ! cover each column containing a starred zero. If n column are covered
   ! the starred zero describe an optimal assignment and we are done otherwise
   ! go to step 4.
   use var
   implicit none

   integer, intent(out) :: step

   integer :: colCount, i, j

   colCount = 0
   do i = 1, n
      do j = 1, n
         ! if starred and column is uncovered
         if (M(j,i) == 1 .and. colCover(j) == 0) then
            colCover(j) = 1
            colCount = colCount + 1
         endif
      enddo
   enddo

   if (colCount == n) then
      step = 0
   else
      step = 4
   endif

end subroutine step3

subroutine step4(step)
   ! Find a uncovered zero and prime it. If there is no starred zero in the row
   ! go to step 5. Otherwise, cover the row and uncover the column containing
   ! the starred zero. Continue until no uncovered zeros is left. Go to step 6.
   use var
   implicit none

   integer, intent(out) :: step

   logical :: done, starInRow
   integer :: i, j, row, col

   done = .false.

   do while (.not. done)
      ! find an uncovered zero
      row = 0; col = 0
      starInRow = .false.
      loop1: do i = 1, n
         loop2: do j = 1, n
            if (CC(j,i) == 0 .and. rowCover(i) == 0 .and. colCover(j) == 0) then
               row = i
               col = j
               exit loop1
            endif
         enddo loop2
      enddo loop1
      if (row == 0) then !no zero uncoverred left
         done = .true.
         step = 6
      else
         M(col,row) = 2 !primed zero
         ! search if there is a starred zero in the same row
         do j = 1, n
            if (M(j,row) == 1) then
               starInRow = .true.
               col = j
            endif
         enddo
         if (starInRow) then ! if there is a starred zero in line
            rowCover(row) = 1
            colCover(col) = 0
         else ! if no starred zero in line
            done = .true.
            step = 5
            pathRow0 = row
            pathCol0 = col
         endif
      endif
   enddo

end subroutine step4

subroutine step5(step)
   ! Augmenting path algorithm: construct a serie of alternating primed and
   ! starred zeros as follows. Let Z0 be the uncoverd primed zero found in
   ! step 4. Let Z1 be the starred zero in the column of Z0 (if any).
   ! Let Z2 be the primed zero in the row of Z1 (there will always be one).
   ! Continue until the series terminates at a primed zero that has no starred
   ! zero in its column. Then unstar each starred zeros of the series, star
   ! each primed zeros of the series, erase all primes and uncover every line
   ! and columns. Return to step 3.
   use var
   implicit none

   integer, intent(out) :: step

   logical :: done
   integer :: i, j
   integer :: row, col
   integer :: pathCount
   integer, dimension(2*n+1,2) :: path

   pathCount = 1

   path(pathCount,1) = pathRow0
   path(pathCount,2) = pathCol0

   done = .false.

   do while (.not. done)
      ! search for a starred zero in column
      row = 0
      col = path(pathCount,2)
      do i = 1, n
         if (M(col,i) == 1) row = i
      enddo
      if (row /= 0) then ! update path
         pathCount = pathCount + 1
         path(pathCount,1) = row
         path(pathCount,2) = path(pathCount-1,2)
      else
         done = .true.
      endif
      if (.not. done) then
         ! search for a prime zero in row
         do j = 1, n
            if (M(j,row) == 2) col = j
         enddo
         ! update path
         pathCount = pathCount + 1
         path(pathCount,1) = path(pathCount-1,1)
         path(pathCount,2) = col
      endif
   enddo

   ! augment path
   do i = 1, pathCount
      if(M(path(i,2),path(i,1)) == 1) then
         M(path(i,2),path(i,1)) = 0
      else
         M(path(i,2),path(i,1)) = 1
      endif
   enddo

   ! clear covers and erase primes
   do i = 1, n
      rowCover(i) = 0
      colCover(i) = 0
      do j = 1, n
         if (M(j,i) == 2) M(j,i) = 0
      enddo
   enddo

   step = 3

end subroutine step5

subroutine step6(step)
   ! Search for the smallest uncovered value and add it to covered rows
   ! and substract it from uncovered columns. Return to step 4.
   use var
   implicit none

   integer, intent(out) :: step

   integer :: i, j, minVal

   minVal = huge(i)

   do i = 1, n
      do j = 1, n
         if (rowCover(i) == 0 .and. colCover(j) == 0 .and. CC(j,i) < minVal) then
            minVal = CC(j,i)
         endif
      enddo
   enddo
   do i = 1, n
      do j = 1, n
         if (rowCover(i) == 1) CC(j,i) = CC(j,i) + minVal
         if (colCover(j) == 0) CC(j,i) = CC(j,i) - minVal
      enddo
   enddo

   step = 4

end subroutine step6
