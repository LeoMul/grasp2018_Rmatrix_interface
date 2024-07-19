program orbsdump
    
    implicit none

    integer :: i,j,II,stat
    real*8 :: x,NUCCHARGE
    character*6 :: header

    real*8 :: RNT,H
    character*2 :: angular_symbols_neg(5)
    character*2 :: angular_symbols_pos(4)

    integer,parameter :: NUM_ORBITALS_FIXED   = 50
    integer,parameter :: NUM_RAD_POINTS_FIXED = 500

    integer :: NP(NUM_ORBITALS_FIXED),NAK(NUM_ORBITALS_FIXED),MFJ(NUM_ORBITALS_FIXED)
    real*8  :: E(NUM_ORBITALS_FIXED)
    real*8  ::PZ(NUM_ORBITALS_FIXED)
    real*8  ::QZ(NUM_ORBITALS_FIXED)
    character*2 :: NH_ARRAY(NUM_ORBITALS_FIXED)

    real*8,allocatable :: large_component(:,:),small_component(:,:),radial_grid(:)
    real*8,allocatable :: generalised_occupation(:)

    integer ::num_orbitals, num_points 

    logical :: end_of_file
    
    angular_symbols_neg = (/'S ','P ','D ','F ','G ' /)
    angular_symbols_pos = (/'P-','D-','F-','G-' /)
    
    write(*,100)
    print*,'Entering routine get_nuclear_charge'   
    call get_nuclear_charge(NUCCHARGE)
    write(*,100)
    print*,'Entering routine read_orbitals'
    call read_orbitals()    
    write(*,100)
    print*,'Entering routine write_out'
    call write_out()
    write(*,100)
    print*,'STOPPING NORMALLY'

    100 format('------------------------------------------------------------------------')

    contains

    subroutine write_out()

        real*8 :: hy,RNT,ey,zy,exph,ratio
        real*8 :: b,c
        character*2 :: NH

        open(2,file='ORBIN.DAT',form='UNFORMATTED')
        print*,'WRITING TO: ORBIN.DAT'
      
        zy = NUCCHARGE 
        !print*,zy
        ratio = radial_grid(3) / radial_grid(2)

        b = -ratio 
        c = ratio - 1.0 

        exph = 0.5d0*(-b + sqrt(b*b - 4.0d0 * c ))
        hy = log(exph)

        rnt = radial_grid(2) / (exph -1.0d0)

        do i = 1,num_orbitals


            WRITE(2)     NH_ARRAY(i),NP(i),NAK(i),MFJ(i),ZY,HY,rnt,E(I),PZ(i),QZ(i)
            WRITE (2) (large_component(I,j),j=1,MFJ(i)),(small_component(I,j),j=1,MFJ(i))

        end do 

        close(2)


    end subroutine

    subroutine read_orbitals()
        logical :: rwfn_existence
        REAL*8 :: QZII
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
                READ (1) PZ(II), (large_component(ii,I),I=1,MFJ(ii)), (small_component(ii,I),I=1,MFJ(ii))       !, (small_component(I),I=1,MFJ)
                READ (1) (radial_grid(I),I=1,MFJ(ii))
                !print*,NUCCHARGE
                QZ(II) = qzfunc(PZ(II),NAK(II),NUCCHARGE)
                
                if     (NAK(ii) .ge. 0) then 
                    NH_ARRAY(ii) = angular_symbols_pos(NAK(ii))
                elseif (NAK(ii) .le. 0) then
                    NH_ARRAY(ii) = angular_symbols_neg(abs(NAK(ii)))
                else   
                    print*,'invalid kappa in orbital',ii
                    stop 
                end if 

               write(*,1010) NP(II),NH_ARRAY(II),NAK(II),E(II),MFJ(ii),PZ(ii), QZ(II)

            else
                write(*,1030)num_orbitals
                exit
            end if
        end do 

        close(1)


        num_points = maxval(MFJ) 

        1010 format(I6,2X,A3,1X,I7,3X,es10.4,1X,I7,1X,F13.4,1X,F13.4)
        1020 format(' PrincN',' L(-)',1X,' Kappa',1X,'  Eigenv(Ha)',1X,'NPoints',12X,'PZ',12X,'QZ')
        1030 format(I4,1X,'relativistic orbials found - exiting grasp2018 read')
    end subroutine
    
    subroutine get_nuclear_charge(nuclear_charge)
        real*8,intent(inout) :: nuclear_charge
        logical :: isodata_existence
        inquire(file="isodata",exist=isodata_existence)
        
        if (.not. isodata_existence) then 
            stop 'isodata not found - did you save it to something else? are you in the correct directory?'
        end if

        open(10,file='isodata')
        read(10,*)
        read(10,*) nuclear_charge
        close(10)

        write(*,100) nuclear_charge

        100 format(' Found nuclear charge: ',F12.2)


    end subroutine get_nuclear_charge


    function qzfunc(pz,kappa,Z)
        real*8 :: pz ,qzfunc,factor
        integer :: kappa
        real*8 :: Z 
        real*8,parameter :: c = 137.03599d0
        real *8 :: gama 

        gama = sqrt(kappa*kappa - (Z/C)**2 )

        if (kappa .ge. 0.0d0) then 
            factor =  (C/Z) * 1.0d0 / (kappa + gama)
        else  
            factor =  (Z/C) * 1.0d0 / (kappa - gama)
        end if 
        !print*,kappa,gama,Z/C,factor
        qzfunc = pz * factor

    end function


end program orbsdump