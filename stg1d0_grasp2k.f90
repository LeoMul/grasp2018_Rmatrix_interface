program read 
    !stg2d0_grasp2k.f90 
    !gfortran stg1d0_grasp2k.f90 -o stg1d0_grasp2k.x
    !author: Leo Patrick Mulholland, QUB PhD Student, ARC
    
    !reads in:
    !rwfn.out, rmcdhf.sum from grasp2018 runs
    !into the format employed by the DARC R-matrix codes (Norrington, Berrington, Ballance)
    
    !currently fixed dimensions:
    !max number of orbitals
    !max num of radial points
    !right now they are unreasonably large numbers but should be updated for 
    ! a fully dynamic allocation.

    !Is this a particularly elegant code? Not really. Does it work? Yes.
    !This code serves really as a proof of concept that it is quite easy
    !to employ the large atomic models used by other atomic physicists
    !to produce R-matrix data.

    !In principle, their large models can be used to systematically produce good orbitals
    !With large amounts of CI in the eigenvectors.
    !If calculation size becomes a problem  (R-matrix is O(N^3) afterall), the important 
    !CSFS can be isolaed and used in conjunction with the fully optimised orbitals.

    !A point to consider is the implementation of pseudo-correlation-orbitals in the main
    !version of grasp. There are orbitals with large principal quantum numberincluded to 
    !help convergence of lower ones. The problem is these orbitals tend to have the incorrect
    !number of nodes, as this restriction is not enforced on them. Only the so-called spectroscopic
    !orbitals have this constraint. I have not tested this in R-matrix - but imagine results involving
    !excitation to those orbitals would probably be a little unphysical in terms of Collision strength
    !and so forth. 

    implicit none

    integer :: i,j,II,stat
    real*8 :: x
    character*6 :: header

    integer,parameter :: NUM_ORBITALS_FIXED   = 50
    integer,parameter :: NUM_RAD_POINTS_FIXED = 500

    integer :: NP(NUM_ORBITALS_FIXED),NAK(NUM_ORBITALS_FIXED),MFJ(NUM_ORBITALS_FIXED)
    real*8 :: E(NUM_ORBITALS_FIXED)
    real*8,allocatable :: large_component(:,:),small_component(:,:),radial_grid(:)
    real*8 ::PZ

    real*8,allocatable :: generalised_occupation(:)

    integer ::num_orbitals, num_points 

    logical :: end_of_file

    write(*,100)
    print*,'Entering routine read_orbitals'
    call read_orbitals()    
    
    write(*,100)
    print*,'Entering routine read_occupation_numbers'

    call read_occupation_numbers(num_orbitals)
    
    write(*,100)
    print*,'Entering routine write_out'
    call write_out()
    write(*,100)
    print*,'STOPPING NORMALLY'

    100 format('------------------------------------')

    contains

    subroutine write_out()

        !writes TARGET.INP
        !this file contains:
        !generalised occupation numbers
        !the largest radial grid used to store the orbitals
        !the relativistic orbitals, 
                !n, kappa
                !large component 
        !then 
                !n, kappa
                !small component 
        
        !where the components have been extended by zeros to fit into the complete grid described above.


        open(2,file='TARGET.INP')
        print*,'Writing to TARGET.INP'
        write(2,3040) num_orbitals,num_points
        WRITE (2,3050) (generalised_occupation(i),i=1,num_orbitals)
        WRITE (2,3050) (radial_grid(I),I=1,num_points)

        do i = 1,num_orbitals
            WRITE (2,3060) NP(i),NAK(i)
            WRITE (2,3050) (large_component(I,j),j=1,num_points)
            WRITE (2,3060) NP(i),NAK(i)
            WRITE (2,3050) (small_component(I,j),j=1,num_points)
        end do 

        close(2)

        3040 FORMAT (1X,2I7)
        3050 FORMAT (1X,1P,4E16.8)
        3060 FORMAT (1X,I4,2X,I4)
    end subroutine

    subroutine read_orbitals()
        logical :: rwfn_existence

        inquire(file="rwfn.out",exist=rwfn_existence)


        if (.not. rwfn_existence) then 
            stop 'rwfn.out not found - did you save it to something else? are you in the correct directory?'
        end if 

        open(1,file = "rwfn.out",form='unformatted',status='old')

        !allocate the matrices with the fixed dimensions
        !there is probably a more elegant way to do this
        !however, it is unlikely for there to be more than 50 relativistic orbitals
        allocate(large_component(NUM_ORBITALS_FIXED,NUM_RAD_POINTS_FIXED))
        allocate(small_component(NUM_ORBITALS_FIXED,NUM_RAD_POINTS_FIXED))
        allocate(radial_grid(NUM_RAD_POINTS_FIXED))
    
        large_component = 0.0d0
        small_component = 0.0d0 

        READ (1) header
        write(*,1020)

        !I need to add an exception here for when the number of orbitals is larger
        !than NUM_ORBITALS_FIXED
        !either that or have a dynamic allocation should another be detected...
        

        !sketch of new algorithm:
        !do while iostat == True 
        !read NP(II),NAK(II),E(II),MFJ(ii)
        !extend allocation to allow for  another orbital
        !read orbital
        !repeat

        DO II = 1,NUM_ORBITALS_FIXED
            READ (1, iostat = j) NP(II),NAK(II),E(II),MFJ(ii)
            if (j .ne. -1) then 
                num_orbitals = num_orbitals + 1
                write(*,1010) NP(II),NAK(II),E(II),MFJ(ii)
                READ (1) PZ, (large_component(ii,I),I=1,MFJ(ii)), (small_component(ii,I),I=1,MFJ(ii))       !, (small_component(I),I=1,MFJ)
                READ (1) (radial_grid(I),I=1,MFJ(ii))
            else
                write(*,1030)num_orbitals
                exit
            end if
        end do 

        close(1)


        num_points = maxval(MFJ) 

        1010 format(I6,1X,I6,1X,F12.6,1X,I7)
        1020 format('PrincN',1X,' Kappa',1X,'  Eigenv(Ha)',1X,'NPoints')
        1030 format(I4,1X,'relativistic orbials found - exiting grasp2018 read')
    end subroutine

    subroutine read_occupation_numbers(num_orbitals)
        
        integer,intent(in) :: num_orbitals

        logical :: finished

        character*5 :: dummy_character
        real*8 :: dummy_real(5)
        !real*8, allocatable :: generalised_occupation(:)

        character*20 :: string_detector
        integer :: counter ,ii,jj
        integer :: subshell_mention_counter
        integer , parameter :: max_read = 100000

        logical :: sum_existence
            
        inquire(file="rmcdhf.sum",exist=sum_existence)

        if (.not. sum_existence) then 
            stop 'rmcdhf.sum not found - did you save it to something else? are you in the correct directory?'
        end if 
        open(3,file='rmcdhf.sum')


        finished=.false.
        counter = 0
        subshell_mention_counter = 0
        allocate(generalised_occupation(num_orbitals))

        do while (finished .eqv. .false.)
            counter = counter + 1

            read(3, *) string_detector
            !print*,string_detector


            !generalised occ is after the second mention of Subshell
            !in the 2018 version of grasp (Fischer version)
            if (string_detector .eq. 'Subshell') then 
                subshell_mention_counter = subshell_mention_counter + 1
            end if 

            if (subshell_mention_counter .gt. 1) then 
                exit
            end if 

            if (counter .gt. max_read) then 
                stop 'max read reached, check rmcdhf.sum file'
                exit
            end if 

        end do 

        !reading in the generalised occupation numbers.
        !dummy character, dummy reals to just read in the go
        do ii = 1,num_orbitals
            read(3,*) dummy_character ,  (dummy_real(jj),jj=1,5) ,generalised_occupation(ii)
        end do 
        
        write(*,1020) 
        do ii=1,num_orbitals
            write(*,1010) NP(II),NAK(II),generalised_occupation(ii)
        end do

        close(3)

        1010 format(I6,1X,I6,1X,F12.6)
        1020 format('PrincN',1X,' Kappa',1X,'    GenOCC')
    end subroutine

end program 