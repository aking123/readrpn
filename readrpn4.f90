   module get_cfgs

    character(len=150)                :: outrep, outname
    character(len=4), dimension(100)  :: out_dyn_L, out_phy_L, out_chm_L
    integer                           :: initStep, finaStep, Step_rsti
    integer                           :: out_freq
    character(len=10)                 :: out_dateo
    logical                           :: do_zonal_ave

    logical                           :: do_2D, do_3D, not_done


    integer                           :: maxVars, ndyn, nphy, nchm

! Basic RPN file attributes 
    integer          :: ni, nj, nk, ip1, ip2, ip3
    integer          :: ig1, ig2, ig3, ig4
    integer          :: dateo, deet, npas, nbits, datyp, swa, lng, dltf,     &
                        ubc, extra1, extra2, extra3
    character        :: typvar, grtyp
    character(len=4) :: nomvar 
    character(len=8) :: etiket

!****** Global variables
    type :: myvars
      character(len=4)                :: varname
      integer                         :: nlevs
      integer                         :: unf
      real, dimension(:,:,:), pointer :: val
    end type myvars

    type(myvars), allocatable, dimension(:)  :: all_vars
    target all_vars

! Vertical levels ip1 lists and vertical hybrid coordinates
    integer, dimension(:), pointer    :: ip1t
    real, dimension(:), pointer       :: hybt

! Geographic latitude and longitude
    real, dimension(:), pointer       :: lat, lon

! Total number of distinct variables 
    integer nvars

   contains

!*************************************************************
      integer function read_nml()
       implicit none

       namelist /read_cfgs/ outrep    ! Directory containing GEM output folders
       namelist /read_cfgs/ outname   ! Filename (prefix) for retrieved data
       namelist /read_cfgs/ out_dyn_L ! List of Dynamics variables to read
       namelist /read_cfgs/ out_phy_L ! List of Physics variable to read
       namelist /read_cfgs/ out_chm_L ! List of Chemistry variable to read
       namelist /read_cfgs/ initStep  ! The model step number to start retrival 
       namelist /read_cfgs/ finaStep  ! The last model step number to retrieve
       namelist /read_cfgs/ Step_rsti ! Number of step before restart
       namelist /read_cfgs/ out_freq
       namelist /read_cfgs/ out_dateo
       namelist /read_cfgs/ do_zonal_ave

       integer nrec, iun, ier
       integer i

       integer  fnom
       external fnom

       read_nml = -1

! Initialize the namelist variables
       outrep       = ''
       outname      = 'outputname'
       initStep     = 0
       finaStep     = 0
       Step_rsti    = 10
       out_freq    = 1
       out_dateo    = '2009000000'       !date and time of origin
 
       out_dyn_L    = ''
       out_phy_L    = ''
       out_chm_L    = ''

       do_zonal_ave = .false.

! Also initialize the and do_ variables here
       do_2D = .false. ; do_3D = .false. ; not_done = .true.

       iun = 0
       ier = fnom (iun, 'read_settings.nml', 'SEQ+OLD' , nrec)
       if ( ier == 0 ) then
          rewind(iun)
          read (iun, nml=read_cfgs, end=9120, err=9120)
          call fclos (iun)
       else
          write(*,*) 'ALERT !! The settings file for controlling the output does not exist '
       endif

! Get the number of non-null values in the variable set
       ndyn = 0
       nphy = 0
       nchm = 0
       do i = 1, 100
         if (out_dyn_L(i) /= '') ndyn = ndyn + 1
         if (out_phy_L(i) /= '') nphy = nphy + 1
         if (out_chm_L(i) /= '') nchm = nchm + 1
       enddo
       maxVars = ndyn + nphy + nchm
       if (maxVars == 0) then
         write(*,*) 'No model output variable requested. '
         read_nml = -1
         return
       endif  

       read_nml = 1

       return
 9120  call fclos (iun) 

      end function read_nml

!*************************************************************
!** This routine is used to sort the variables obtained from the namelist,
!** and produce only distinct that are also really exists in the RPN filelist.
!** It also returns the latitude and longitudes values in model output
!*************************************************************
      integer function compact(flnom, iun, vars)
       implicit none

       integer iun
       character(len=200)                    :: flnom  ! RPN fullpath filename
       character(len=4), dimension(maxVars)  :: vars, res
       integer, dimension(maxVars)           :: ires
       integer i, j, k, nVar

       integer, parameter                    :: nkmax = 200
       integer, dimension(nkmax)             :: list
       integer  ier, nlist, nx, ny, nz
       logical  file_exist

       integer  fstinl, fstprm, fstlir, fnom, fstouv, fstfrm
       external fstinl, fstprm, fstlir

       compact = -1

       inquire(file=flnom, exist=file_exist)
       if ( .not. file_exist ) then
         write(*,*)'ALERT !! The model output file ', trim(flnom), ' does not exist'
         return
       endif
       ier = fnom(iun, flnom, 'STD+RND+R/O', 0)
       if (ier < 0) print *, 'error while opening'
       ier = fstouv(iun, 'STD+RND')

! Find the first non-null and valid variable in the set
       k = 1
       nVar = 0
       do i = 1, maxVars
         k = k + 1
         if (vars(i) /= '') then

           ier = fstinl(iun, nx, ny, nz, -1, ' ', -1, -1, -1, '',       &
                                vars(i), list, nlist, nkmax)
           if (nlist == 0 ) then
             write(*,*)'ALERT !! ',vars(i),' is not present in the model output'
             cycle
           endif

           nVar = 1
           res(1)  = vars(i)
           ires(1) = nlist
           ier     = fstprm(list(1), dateo, deet, npas, ni, nj, nk, nbits,    &
                            datyp, ip1, ip2, ip3, typvar, nomvar, etiket,     &
                            grtyp, ig1, ig2, ig3, ig4, swa, lng, dltf,        &
                            ubc, extra1, extra2, extra3)
           exit

         endif
       enddo

       if (nVar > 0) then   !  None of the requested variable is present

       outer: do i = k, maxVars
         do j = 1, nVar
           if (res(j) == vars(i)) then
           ! Found a match so start looking again
             cycle outer
           end if
         end do
       ! No match found, check if it is not null, then check if it actually
       ! exists in the file and if so add it to the output
         if ( vars(i) /= '' ) then

           ier = fstinl(iun, nx, ny, nz, -1, ' ', -1, -1, -1, '',       &
                                vars(i), list, nlist, nkmax)
           if (nlist == 0 ) then
             write(*,*)'ALERT !! ',vars(i),' is not present in the model output'
             cycle outer
           endif
           nVar       = nVar + 1
           res(nVar)  = vars(i)
           ires(nVar) = nlist
         endif
       end do outer

! Now copy the valid sets and their properties to the output variable 
       do i = 1, nVar
         k = nvars + i 
         all_vars(k)%varname = res(i)
         all_vars(k)%nlevs   = ires(i)
         all_vars(k)%unf     = iun

! Allocate space to store the retrieved variables
         allocate( all_vars(k)%val(ni,nj, all_vars(k)%nlevs) )
       enddo
 
! Retrieve geographic latitude and longitude values.
       if ( nvars == 0 ) then     ! nvar == 0 ensures this is done just once
         allocate( lat(nj), lon(ni) )
         ier = fstlir(lat, iun, nx, ny, nz, -1, '', -1, -1, -1, 'X', '^^')
         ier = fstlir(lon, iun, nx, ny, nz, -1, '', -1, -1, -1, 'X', '>>')
       endif

! Get the ip1 list and vertical coordinates (if necessary)
       ier  = get_ip1_list(iun, all_vars((nvars+1):(nvars+nVar)), nVar)
       if (ier < 0) then
        write(*,*) ' Retrieving IP1 lists unsuccessful. '
        return
       endif
       nvars = nvars + nVar

       compact = 1
       endif

! Close the openned RPN files
       ier = fstfrm(iun)
       call fclos(iun)

       return
      end function compact

!*************************************************************
!** Get the Level ip1 lists for and the thermo/ momentum staggered vertical
!** coordinates, as well as the ip1 list for the regolith levels
!*************************************************************
      integer function get_ip1_list(iun, vars, vsize)

      use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,VGD_OK
       implicit none 

       integer iun, vsize
       type(myvars), dimension(vsize)  :: vars

       integer i, k, ier
       integer, dimension(100)         :: list
       integer  fstinl, fstprm, nlist, nx, ny, nz
       external fstinl, fstprm

       type(vgrid_descriptor)          :: vgd

!*** The 'do_*' logical variables used here are global variables defined in the
!*** get_cfgs module. 
!*** The ip1 levels produced are also globally defined in the module variables
!*** declaration

       get_ip1_list = -1

       do i = 1, vsize
! The do_* are global variables whose values are 'saved' in get_cfgs module

! for 2D fields
         if ( vars(i)%nlevs == 1 ) then
            do_2D  = .true.
            vars(i)%stag = .false.
         endif

! for (Atmospheric) 3D fields
         if ( vars(i)%nlevs > 14 ) then
            if ( .not. do_3D ) then
!           ier = fstinl(iun, ni, nj, nk, -1, ' ', -1, -1, -1, '',       &
!                                '!!  ', list, nlist, nkmax)
!           write(*,*) ' found vcood with ',vars(i), nlist
            ier = vgd_new(vgd,unit=iun,format='fst',ip1=-1,ip2=-1)
            if (ier /= VGD_OK) then
               print*, 'compact error: cannot access the coordinate descriptor'
               return
            endif
            if (vgd_get(vgd,'VIPM - ip1 list (m)', ip1m) /= VGD_OK) ier = -1
            if (vgd_get(vgd,'VIPT - ip1 list (t)', ip1t) /= VGD_OK) ier = -1
            if (vgd_get(vgd,'VCDM - vert. coord. (m)',hybm) /= VGD_OK) ier = -1
            if (vgd_get(vgd,'VCDT - vert. coord. (t)',hybt) /= VGD_OK) ier = -1

            do_3D  = .true.
            endif

            if ( vars(i)%nlevs == size(ip1t) ) vars(i)%stag = .false.
            if ( vars(i)%nlevs == size(ip1m) ) vars(i)%stag = .true.
         endif
       enddo
       get_ip1_list = 1

      end function get_ip1_list

!*************************************************************
!** Create a grads program description file for all binary output of readrpn
!*************************************************************
      subroutine grads_description(vars, nvars, ct)
        integer nvars, ct
        type(myvars), dimension(nvars) :: vars

        integer  n, k
        integer  nu, nz, nv
        integer  n_2d, n_3d
        real, allocatable, dimension(:):: levels
        character(len=4)               :: nom
        character(len=150)             :: filename
        logical                        :: once_2d = .false., once_3d = .false.,&
                                          done = .false.  


        n_2d = 0 ; n_3d = 0
        do n = 1, nvars
! For 2-D fields
          if (vars(n)%nlevs == 1) then
             n_2d = n_2d + 1
             if ( .not. once_2d ) then
               nu = 23
               nz = vars(n)%nlevs
               allocate(levels(nz))
               levels = (/1000.0/)
               filename = trim(outname)//'_2d'
               open(unit=nu, file=trim(filename)//'.ctl', form='formatted')
               once_2d = .true.
             else
               cycle
             endif
! For 3-D fields
          elseif (vars(n)%nlevs >= size(ip1m)) then
             n_3d = n_3d + 1
             if ( .not. once_3d ) then
               nu = 24
               nz = size(ip1t)
               allocate(levels(nz))
               levels = hybt * 1000.0
               filename = trim(outname)//'_3d' 
               open(unit=nu, file=trim(filename)//'.ctl', form='formatted')
               once_3d = .true.
             else
               cycle
             endif
          endif
          write (nu,'(A8)')        'DSET %ch'
          write (nu,'(A8,I5,A20)') 'CHSUB 1 ',ct,trim(filename)//'.dat'
          write (nu,'(A12)')       'UNDEF -999.0'
          write (nu,'(A16)')       'OPTIONS template'
          write (nu,'(A18)')       'OPTIONS sequential'
          write (nu,'(A14)')       'TITLE variable'
          write (nu,105)           'XDEF', ni, ' LINEAR',lon(1), 360.0/(ni-1)
          write (nu,105)           'YDEF', nj, ' LINEAR',lat(1), 180.0/nj
          write (nu,'(A5,I3,A7,1X,F15.6)')  'ZDEF ',nz, ' LEVELS', levels(nz)
          do k = nz-1, 1, -1
            write (nu,'(5X,F15.6)') levels(k)
          enddo
          if (allocated(levels)) deallocate(levels)
          write (nu,'(A5,I5,A22)')  'TDEF ',ct,' LINEAR 0JAN0001 15MN'
        enddo
        if ( once_2d  )  write (23,'(A4,1X,I2)') 'VARS', n_2d
        if ( once_3d  )  write (24,'(A4,1X,I2)') 'VARS', n_3d

        do n = 1, nvars
           nom = vars(n)%varname
! For 2-D fields
          if (vars(n)%nlevs == 1) then
             nu = 23 ; nv = 0
! For 3-D fields
          elseif (vars(n)%nlevs >= size(ip1t)) then
             nu = 24 ; nv = size(ip1t)
          endif
          write (nu,'(A5,2X,I3,2X,I2,A15)') nom, nv, 99, ' VARIABLE_'//nom
        enddo
        if ( once_2d  )  write (23,'(A7)') 'ENDVARS'
        if ( once_3d  )  write (24,'(A7)') 'ENDVARS'
    

 105    format (A5,I3,A8,2(1X,F5.1))
      return

      end subroutine grads_description

   end module get_cfgs


!*************************************************************
!** Main program 
!*************************************************************
   program readrpn
    use get_cfgs
    implicit none

    character(len=200)                    :: flnom  ! RPN fullpath filename
    character(len=30)                     :: rpn_name    ! Filename suffix


! Unit numbers associated with the RPN file
    integer, parameter                    :: iun_D=14, iun_P=16, iun_K=18
! Arrays of variables to read from the rpn file(s)
    type(myvars), pointer, dimension(:)   :: fields

! Work fields
    real, allocatable, dimension(:,:)     :: ps  !Alternate surf. presurre
    real, allocatable, dimension(:,:)     :: fvalue !2-D field
    real, allocatable, dimension(:,:)     :: fld_av !Zonally average field
    integer i, j, k, n, l, m, ier, iun, nx, ny, nz, jk
    character(len=4)  :: nom
    logical file_exist, get_dyn_S, get_phy_S, get_chm_S

    integer ct, ct0, ct3
    character(len=12)  :: cnum

    integer  fnom, fstouv, fstfrm, fstinl, fstlir, fstinf
    external fnom, fstouv, fstfrm, fstinl, fstlir, fstinf, fclos

! Dcst_tcdk_8    | .27315e+3            ; conversion from celsius to kelvin |
! Dcst_knams_8   | .514791              ; conversion from knots to m/s |
    real, parameter :: tcdk = 0.27315e+3, knams = 0.514791

    ier = read_nml()
    if (ier < 0) write(*,*) ' Error in reading the retrieval configuration ' 

    get_dyn_S = .false. ; get_phy_S = .false. ; get_chm_S = .false.

    ct    = 0
    mainl: do ct0 = 0, finaStep, out_freq
        write(cnum, '(i3.3)') ct0
        rpn_name   = out_dateo//hr//'_'//trim(cnum)

!** Get the list of distinct variables, and the numbers of levels each contains.
!** This is done once (for the first file in the output folder), with the 
!** assumption that all the files in output folder ('outrep') have the same properties
        if (ct == 0) then
          nvars = 0     ! A global variable for total number of variables read 
          allocate(all_vars(maxVars))          

          flnom = trim(outrep)//'dm'//trim(rpn_name)
          if (ndyn > 0) then
            ier = compact(flnom, iun_D, out_dyn_L)
            if (ier > 0)  get_dyn_S = .true.
          endif

          flnom = trim(outrep)//'pm'//trim(rpn_name)
          if (nphy > 0) then
            ier = compact(flnom, iun_P, out_phy_L)
            if (ier > 0)  get_phy_S = .true.
          endif

          flnom = trim(outrep)//'km'//trim(rpn_name)
          if (nchm > 0) then
            ier = compact(flnom, iun_K, out_chm_L)
            if (ier > 0)  get_chm_S = .true.
          endif

          if (nvars > 0) then
            allocate(fields(nvars))
            fields => all_vars(1:nvars)

            allocate(fvalue(ni,nj))
            if ( do_2D )  &
               open(unit=31, file=trim(outname)//'_2d.dat', form='unformatted')
            if ( do_3D )  &
               open(unit=32, file=trim(outname)//'_3d.dat', form='unformatted')
          else
            stop    ! Model files or the requested variables cannot be found
          endif

        endif

!* Associate unit numbers with the files to be read
        if ( get_dyn_S ) then
          flnom = trim(outrep)//'dm'//trim(rpn_name)
          inquire(file=flnom, exist=file_exist)
          if ( .not. file_exist ) then
             write(*,*) 'The dynamics output file ',flnom , ' does not exist'
             stop
          endif
          ier = fnom(iun_D, flnom, 'STD+RND+R/O', 0)
          if (ier < 0) print *, 'error while opening'
          ier = fstouv(iun_D, 'STD+RND')

          
        endif

        if ( get_phy_S ) then
          flnom = trim(outrep)//'pm'//trim(rpn_name)
          inquire(file=flnom, exist=file_exist)
          if ( .not. file_exist ) then
             write(*,*) 'The physics output file ',flnom, ' does not exist'
             stop
          endif
          ier = fnom(iun_P, flnom, 'STD+RND+R/O', 0)
          if (ier < 0) print *, 'error while opening'
          ier = fstouv(iun_P, 'STD+RND')
        endif

        if ( get_chm_S ) then
          flnom = trim(outrep)//'km'//trim(rpn_name)
          inquire(file=flnom, exist=file_exist)
          if ( .not. file_exist ) then
             write(*,*) 'The chemistry output file ',flnom, ' does not exist'
             stop
          endif
          ier = fnom(iun_K, flnom, 'STD+RND+R/O', 0)
          if (ier < 0) print *, 'error while opening'
          ier = fstouv(iun_K, 'STD+RND')
        endif

        do n = 1, nvars
          iun = fields(n)%unf
          nom = fields(n)%varname

! For 2-D fields
          if ( do_2D .and. fields(n)%nlevs == 1) then
            ier = fstlir(fields(n)%val(:,:,1), iun, nx, ny, nz,           &
                          -1, etiket, -1, -1, -1, typvar, nom)
            if (trim(nom) == 'P0') then
               fields(n)%val(:,:,1) = fields(n)%val(:,:,1) * 100.00
            endif

            write(31) ((fields(n)%val(i,j,1),i=1,ni),j=1,nj )
            cycle
          endif

! For Atmospheric 3-D fields
          nk = size(ip1t)
          if ( do_3D ) then
            do k = nk, 1, -1
              ier = fstlir(fvalue, iun, nx, ny, nz,         &
                          -1, etiket, ip1t(k), -1, -1, typvar, nom)
              if (trim(nom) == 'TT') fvalue = fvalue + tcdk
              if (trim(nom) == 'UU') fvalue = fvalue * knams
              if (trim(nom) == 'VV') fvalue = fvalue * knams
              write(32) ((fvalue(i,j),i=1,ni),j=1,nj )
              fields(n)%val(:,:,k) = fvalue
            enddo

          if ( do_zonal_ave ) then
             if (ct == 0 .and. (.not. allocated(fld_av))) then
                allocate( fld_av(nj,nk) )
             endif

             do k = 1, nk
               do j = 1, nj
                  ct3 = 0
                  fld_av(j,k) = 0.0
                  do i = 1, ni
                     pres(:) = ps(i,j) * hybt(:)
                     if (pres(nk) > slev(k)) then
                        ct3 = ct3 + 1
                        fld_av(j,k) = fld_av(j,k) + fld_ez(i,j,k)
                     endif 
                  enddo
                  if (ct3 == 0) then
                    fld_av(j,k) = -999.0
                  else
                    fld_av(j,k) = fld_av(j,k) / ct3
                  endif
               enddo
             enddo
          endif    ! End do_zonal_ave

          cycle
          endif

        ct = ct + 1

! Close the openned RPN files
        if ( get_dyn_S ) then
          ier = fstfrm(iun_D)
          call fclos(iun_D)
        endif
        if ( get_phy_S ) then
          ier = fstfrm(iun_P)
          call fclos(iun_P)
        endif
        if ( get_chm_S ) then
          ier = fstfrm(iun_K)
          call fclos(iun_K)
        endif

    enddo mainl

    call grads_description(fields, nvars, ct)

   end program readrpn
