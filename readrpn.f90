   module get_cfgs

    character(len=150)                :: outrep, outname
    character(len=4), dimension(100)  :: out_dyn_L, out_phy_L, out_chm_L
    integer                           :: initStep, finaStep, Step_rsti
    integer                           :: out_freq
    character(len=10)                 :: out_dateo
    logical                           :: do_interp, do_zonal_ave

    logical                           :: do_2D, do_3D, do_reg, not_done


    integer                           :: maxVars, ndyn, nphy, nchm

! Basic RPN file attributes 
    integer          :: ni, nj, nk, ip1, ip2, ip3
    integer          :: ig1, ig2, ig3, ig4
    integer          :: dateo, deet, npas, nbits, datyp, swa, lng, dltf,     &
                        ubc, extra1, extra2, extra3
    character        :: typvar, grtyp
    character(len=2) :: nomvar 
    character(len=8) :: etiket

!****** Global variables
    type :: myvars
      character(len=4)                       :: varname
      integer                                :: nlevs
      integer                                :: unf
      integer, allocatable, dimension(:)     :: ip1s
      logical                                :: stag
    end type myvars

    type(myvars), allocatable, dimension(:)  :: all_vars
    target all_vars

! Staggerred levels ip1 lists and vertical coordinates
    integer, dimension(:), pointer    :: ip1m => null(), ip1t => null()
    real, dimension(:), pointer       :: hybm => null(), hybt => null()

! Geographic latitude and longitude
    real, dimension(:), allocatable   :: lat, lon

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
       namelist /read_cfgs/ do_interp
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

       do_interp    = .false.
       do_zonal_ave = .false.

! Also initialize the and do_ variables here
       do_2D = .false. ; do_3D = .false. ; do_reg = .false. ; not_done = .true.

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
!** and produce only distinct items. It also check that variable really exists
!** in the RPN filelist.
!** It also returns the latitude and longitudes values in model output
!*************************************************************
      integer function compact(flnom, iun, vars)
      use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,VGD_OK
       implicit none

       integer iun
       character(len=200)                    :: flnom  ! RPN fullpath filename
       character(len=4), dimension(maxVars)  :: vars, res
       integer, dimension(maxVars)           :: ires
       integer i, j, k, jk, nVar

       integer, parameter                    :: nkmax = 200
       integer, dimension(nkmax)             :: list
       integer, dimension(maxVars,nkmax)     :: handle
       integer  ier, nlist, nx, ny, nz
       logical  file_exist

       type(vgrid_descriptor)          :: vgd

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
           handle(1,:) = list
           
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
           handle(nVar,:) = list
         endif
       end do outer

! Now store the valid sets and their properties 
       do i = 1, nVar
         k = nvars + i

! Allocate space to store generic ip1 levels
         allocate(all_vars(k)%ip1s(ires(i)))

!** Get the Level ip1 lists for generic vertical levels
!** (This is needed for 3D variables that are not from dynamics output)
         do jk = 1, ires(i)
            ier = fstprm(handle(i,jk), dateo, deet, npas, ni, nj, nk, nbits,  &
                         datyp, ip1, ip2, ip3, typvar, nomvar, etiket,        &
                         grtyp, ig1, ig2, ig3, ig4, swa, lng, dltf,           &
                         ubc, extra1, extra2, extra3)
            all_vars(k)%ip1s(jk) = ip1
         enddo

         all_vars(k)%varname = res(i)
         all_vars(k)%nlevs   = ires(i)
         all_vars(k)%unf     = iun
       enddo
 
       if ( nvars == 0 ) then     ! nvars == 0 ensures this is done just once

!** Retrieve geographic latitude and longitude values.
         allocate( lat(nj), lon(ni) )
         ier = fstlir(lat, iun, nx, ny, nz, -1, '', -1, -1, -1, 'X', '^^')
         ier = fstlir(lon, iun, nx, ny, nz, -1, '', -1, -1, -1, 'X', '>>')

!** Get the coordinate descriptors for the thermo/ momentum staggered vertical levels
!        ier = fstinl(iun, ni, nj, nk, -1, ' ', -1, -1, -1, '',       &
!                             '!!  ', list, nlist, nkmax)
!        write(*,*) ' found vcood with ',vars(i), nlist
         ier = vgd_new(vgd,unit=iun,format='fst',ip1=-1,ip2=-1)
         if (ier /= VGD_OK) then
            print*, 'compact error: cannot access the coordinate descriptor'
            return
         endif
         if (vgd_get(vgd,'VIPM - ip1 list (m)', ip1m) /= VGD_OK) ier = -1
         if (vgd_get(vgd,'VIPT - ip1 list (t)', ip1t) /= VGD_OK) ier = -1
         if (vgd_get(vgd,'VCDM - vert. coord. (m)',hybm) /= VGD_OK) ier = -1
         if (vgd_get(vgd,'VCDT - vert. coord. (t)',hybt) /= VGD_OK) ier = -1

         if (ier < 0) then
          write(*,*) ' Retrieving IP1 lists unsuccessful. '
          return
         endif

       endif

       do i = 1, nVar
         k = nvars + i 
! for 2D fields
         if ( all_vars(k)%nlevs == 1 ) then
            do_2D  = .true.
            all_vars(k)%stag = .false.
         endif

! Regolith (hardcoded to 14 levels)
         if ( all_vars(k)%nlevs == 14 ) then
            do_reg = .true.
            all_vars(k)%stag = .false.
         endif

! for (Atmospheric) 3D fields
         if ( all_vars(k)%nlevs > 14 ) then
            do_3D  = .true.

            if ( all_vars(k)%nlevs == size(ip1t) ) all_vars(k)%stag = .false.
            if ( all_vars(k)%nlevs == size(ip1m) ) all_vars(k)%stag = .true.
         endif
       enddo

       nvars = nvars + nVar

       compact = 1
       endif

! Close the openned RPN files
       ier = fstfrm(iun)
       call fclos(iun)

       return
      end function compact

!*************************************************************
!** Create a grads program description file for all binary output of readrpn
!*************************************************************
      subroutine grads_description(vars, nvars, ct)
        integer nvars, ct
        type(myvars), dimension(nvars) :: vars

        integer  n, k
        integer  nu, nz, nv
        integer  n_2d, n_3d, n_reg
        real, allocatable, dimension(:):: levels
        character(len=4)               :: nom
        character(len=150)             :: filename
        logical                        :: once_2d = .false., once_3d = .false.,&
                                          once_reg = .false., done = .false.  


        n_2d = 0 ; n_3d = 0; n_reg = 0
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
! For Regolith fields
          elseif (vars(n)%nlevs == 14) then
             n_reg = n_reg + 1
             if ( .not. once_reg ) then
               nu = 25
               nz = size(vars(n)%ip1s)
               allocate(levels(nz))
               levels = vars(n)%ip1s * 1.0
               filename = trim(outname)//'_reg'
               open(unit=nu, file=trim(filename)//'.ctl', form='formatted')
               once_reg = .true.
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
        if ( once_reg )  write (25,'(A4,1X,I2)') 'VARS', n_reg

        do n = 1, nvars
           nom = vars(n)%varname
! For 2-D fields
          if (vars(n)%nlevs == 1) then
             nu = 23 ; nv = 0
! For 3-D fields
          elseif (vars(n)%nlevs >= size(ip1m)) then
             nu = 24 ; nv = size(ip1t)
! For Regolith fields
          elseif (vars(n)%nlevs == 14) then
             nu = 25 ; nv = size(vars(n)%ip1s)
          endif
          write (nu,'(A5,2X,I3,2X,I2,A15)') nom, nv, 99, ' VARIABLE_'//nom
        enddo
        if ( once_2d  )  write (23,'(A7)') 'ENDVARS'
        if ( once_3d  )  write (24,'(A7)') 'ENDVARS'
        if ( once_reg )  write (25,'(A7)') 'ENDVARS'
    

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
    character(len=180)                    :: dir_prefix  ! Path to flnom
    character(len=180)                    :: rpn_name    ! Filename suffix


! Unit numbers associated with the RPN file
    integer, parameter                    :: iun_D=14, iun_P=16, iun_K=18
! Arrays of variables to read from the rpn file(s)
    type(myvars), pointer, dimension(:)   :: fields

! Work fields
    real, allocatable, dimension(:,:)     :: ps  !Alternate surf. presurre
    real, allocatable, dimension(:,:)     :: fvalue !2-D field
    real, allocatable, dimension(:,:,:)   :: zfvalue !3D fvalue
    real, allocatable, dimension(:,:)     :: fld_av !Zonally average field
    real, allocatable, dimension(:,:,:)   :: fld_ez !Vertically interpolated fl
    real, allocatable, dimension(:)       :: pres, sg, slev
    real, allocatable, dimension(:)       :: wk1, wk2, wk3, wk4, wk5
    integer, allocatable, dimension(:)    :: khi, klo
    integer i, j, k, n, l, m, ier, iun, nx, ny, nz, jk
    character(len=4)  :: nom
    logical file_exist, get_dyn_S, get_phy_S, get_chm_S
    real, parameter                          :: bign = 1.0e31

    integer ct, ct0, ct1, ct3
    character(len=12)  :: cnum0, cnum

    integer  fnom, fstouv, fstfrm, fstinl, fstlir, fstinf
    external fnom, fstouv, fstfrm, fstinl, fstlir, fstinf, fclos

! Dcst_tcdk_8    | .27315e+3            ; conversion from celsius to kelvin |
! Dcst_knams_8   | .514791              ; conversion from knots to m/s |
    real, parameter :: tcdk = 0.27315e+3, knams = 0.514791

    ier = read_nml()
    if (ier < 0) write(*,*) ' Error in reading the retrieval configuration ' 

    get_dyn_S = .false. ; get_phy_S = .false. ; get_chm_S = .false.

    ct    = 0
    ct1   = initStep
! Start the main loop
    mainl: do ct0 = 0, finaStep, Step_rsti
!     The output sub-directory (GEMv4 style)
      write(cnum0, '(i10.10)') ct0 + Step_rsti
      dir_prefix = trim(outrep)//'OUT_laststep_'//trim(cnum0)//'/input/000-000/'

!*     Loop through the RPN files in the current output sub-folder
      do while (ct1 <= ct0+Step_rsti)
        if (ct1 > finaStep) exit mainl
        write(cnum, '(i6.6)') ct1
        rpn_name   = out_dateo//'-00-00_'//trim(cnum)//'p'

!** Get the list of distinct variables, and the numbers of levels each contains.
!** This is done once (for the first file in the output folder), with the 
!** assumption that all the files in output folder ('outrep') have the same properties
        if (ct == 0) then
          nvars = 0     ! A global variable for total number of variables read 
          allocate(all_vars(maxVars))          

          flnom = trim(dir_prefix)//'dm'//trim(rpn_name)
          if (ndyn > 0) then
            ier = compact(flnom, iun_D, out_dyn_L)
            if (ier > 0)  get_dyn_S = .true.
          endif

          flnom = trim(dir_prefix)//'pm'//trim(rpn_name)
          if (nphy > 0) then
            ier = compact(flnom, iun_P, out_phy_L)
            if (ier > 0)  get_phy_S = .true.
          endif

          flnom = trim(dir_prefix)//'km'//trim(rpn_name)
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
            if ( do_reg )  &
               open(unit=33, file=trim(outname)//'_reg.dat', form='unformatted')
          else
            stop    ! Model files or the requested variables cannot be found
          endif

          if ( .not. do_3D) do_interp = .false.
        endif

        if ( do_interp )  get_dyn_S = .true.  ! Because 'P0' is needed for
                                              ! vertical interpolation
!* Associate unit numbers with the files to be read
        if ( get_dyn_S ) then
          flnom = trim(dir_prefix)//'dm'//trim(rpn_name)
          inquire(file=flnom, exist=file_exist)
          if ( .not. file_exist ) then
             write(*,*) 'The dynamics output file ',flnom , ' does not exist'
             stop
          endif
          ier = fnom(iun_D, flnom, 'STD+RND+R/O', 0)
          if (ier < 0) print *, 'error while opening'
          ier = fstouv(iun_D, 'STD+RND')

!* If performing vertical interpolation, check for surface pressure in the
!* model output. If not abort; if yes allocate local variables needed for
!* performing the interpolation
          if ( do_interp ) then
            if (ct == 0) then
             ier = fstinf(iun_D, nx, ny, nz, -1, ' ', -1, -1, -1, '', 'P0')
             if (ier <= 0) then
                write(*,*) 'ALERT!! Surface pressure data are not present in '
                write(*,*) 'the model output. They are needed to perform the '
                write(*,*) 'requested vertical interpolation on pressure levels'
                write(*,*) 'Aborting!!!  '
                do_interp = .false.
                stop
             endif
             allocate (ps(ni,nj))

             open(unit=34, file=trim(outname)//'_it.dat', form='unformatted')

            endif
            ier = fstlir(ps, iun_D, nx, ny, nz, -1, '', -1, -1, -1, '', 'P0')
            ps = ps * 100.0

          endif
          
        endif

        if ( get_phy_S ) then
          flnom = trim(dir_prefix)//'pm'//trim(rpn_name)
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
          flnom = trim(dir_prefix)//'km'//trim(rpn_name)
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
            ier = fstlir(fvalue, iun, nx, ny, nz,           &
                          -1, etiket, -1, -1, -1, typvar, nom)
            if (trim(nom) == 'P0') then
               fvalue = fvalue * 100.00
            endif

            write(31) ((fvalue(i,j),i=1,ni),j=1,nj )
            cycle
          endif

! For Regolith fields
          if ( do_reg .and. fields(n)%nlevs == 14) then
            do k = 14, 1, -1
              ier = fstlir(fvalue, iun, nx, ny, nz,         &
                          -1, etiket, fields(n)%ip1s(k), -1, -1, typvar, nom)
               write(33) ((fvalue(i,j),i=1,ni),j=1,nj )
            enddo
            cycle
          endif

! For Atmospheric 3-D fields
          nk = size(ip1t)
          if ( .not. allocated(zfvalue) ) allocate (zfvalue(ni,nj,nk))
          if ( do_3D .and. fields(n)%nlevs >= (nk - 1)) then
          if ( fields(n)%stag ) then             ! Staggered levels first
            do jk = nk, 2, -1
              k   = jk - 1
              ier = fstlir(fvalue, iun, nx, ny, nz,         &
                          -1, etiket, fields(n)%ip1s(k), -1, -1, typvar, nom)
              if (trim(nom) == 'UU') fvalue = fvalue * knams
              if (trim(nom) == 'VV') fvalue = fvalue * knams
              write(32) ((fvalue(i,j),i=1,ni),j=1,nj )
              zfvalue(:,:,k) = fvalue
            enddo
            write(32) ((fvalue(i,j),i=1,ni),j=1,nj )  ! Dummy write
            zfvalue(:,:,k) = fvalue
!                                                  Thermo levels next
          else
            do k = nk, 1, -1
              ier = fstlir(fvalue, iun, nx, ny, nz,         &
                          -1, etiket, fields(n)%ip1s(k), -1, -1, typvar, nom)
              if (trim(nom) == 'TT') fvalue = fvalue + tcdk
              write(32) ((fvalue(i,j),i=1,ni),j=1,nj )
              zfvalue(:,:,k) = fvalue
            enddo
          endif

          if ( do_interp ) then
            if (ct == 0 .and. (.not. allocated(fld_ez))) then
             allocate( fld_ez(ni,nj,nk), fld_av(nj,nk), pres(nk), slev(nk),  &
                       sg(nk), klo(nk), khi(nk), wk1(nk), wk2(nk), wk3(nk),  &
                       wk4(nk), wk5(nk) )
 
             slev = hybt * 1000.0    ! Pre-defined interpolated pressure levels
            endif

            do j = 1, nj
              do i = 1, ni
                 pres(:) = ps(i,j) * hybt(:)
                 call spline1(pres, zfvalue(i,j,:), nk, bign, bign, sg)
                 call splint(pres, nk, slev, wk1, wk2, wk3, wk4, wk5, klo, khi)

                 do k = 1, nk
                    if (pres(nk) > slev(k)) then
                       l = klo(k)
                       m = khi(k)
                       fld_ez(i,j,k) = wk1(k) * zfvalue(i,j,l)       &
                                      +  wk3(k) * zfvalue(i,j,m)     &
                                   + (wk2(k) * sg(l) + wk4(k) * sg(m)) * wk5(k)
                    else
                       fld_ez(i,j,k) = -999.0
                    endif
                 enddo

              enddo
            enddo

            do k = nk, 1, -1
              write(34) ((fld_ez(i,j,k),i=1,ni),j=1,nj )
            enddo
 
          endif    ! End do_interp

          if ( do_zonal_ave ) then
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

        enddo
        ct = ct + 1
        ct1 = ct1 + out_freq

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

      enddo
    enddo mainl

    call grads_description(fields, nvars, ct)

   end program readrpn
