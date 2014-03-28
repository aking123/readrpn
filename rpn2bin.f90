! Transfer RPN data to binary or other data
! nlev,list,flnom, may be changed.

 program rpn2grads 

   implicit none
   integer i, j, k, ct, ct3
   integer :: ni, nj, nk, nk2, nk3
   integer, parameter :: nlev = 102, nlev2 = 14, nlev3 = 12
   real, parameter :: deg2rad = 57.295779515, grav = 3.69, rgasd = 188.95
   real, parameter :: big_n = 1.0e31, dt = 616.49484, bz = 1.3807E-23
! nlev is the number of vertical levels requested for output in outcfg.out 
! It may or may not be equal to numbers of vertical layers in the model.
! ntim, is the number of output files.
!_____________________________________________________

   character(10) :: cnum, cnum0
   character(20) :: rpn_prefix
   integer       :: ct0, ct1, tmstep

   real, dimension(nlev) :: sig, sig2, sg2, sg3, sg4, sg5, sg6, sg7

   integer :: l, m, nlev_out
   integer, dimension(nlev) :: khi, klo
   real, dimension(nlev)    :: a, a2, b, b2, h2
   real :: dz

   real, allocatable, dimension(:,:,:) :: te, aq, pr, uw, vw, ic, op, hgt, tnd
   real, allocatable, dimension(:,:,:) :: tr_col, o3v, o2v, o1d, o3p, cov, aqc
   real, allocatable, dimension(:,:,:) :: h1v, h2v, ohv, ho2v, h2o2v, o2dv
   real, allocatable, dimension(:,:,:) :: ite, iaq, iuw, ivw, iic, iop
   real, pointer,     dimension(:,:,:) :: pte, paq, puw, pvw, pic, pop
   target :: te, aq, uw, vw, ic, op
   target :: ite, iaq, iuw, ivw, iic, iop
   real, allocatable, dimension(:,:,:) :: tsl, rwv, rad, ric, rtw, tep
   real, allocatable, dimension(:,:)   :: phi, work, sn02
   real, allocatable, dimension(:,:)   :: ps, ts, cw, sn, wo, hb, tm, qs, qa
   real, allocatable, dimension(:,:)   :: itez, iuwz, ivwz,iaqz, iicz, iopz
   real, allocatable, dimension(:)     :: lat, lon

   logical :: vert_pro = .false., interp = .false., calc_hgt = .false.,      &
              soil_vars = .false., column = .false., Schm_chem = .false.,    &
              surf_params = .false., Schm_phys = .false., chem_col = .false.,&
              dyna_outp = .false.

! --- RPN file input/output variables

   character :: typvar, grtyp
   character(2) :: nomvar 
   character(8) :: etiket
   character(120) :: flnom, flnom2, flnom3, direc, direc2,              &
                     dir_prefix, dir_prefix2
   integer :: ier, iun = 16, iun2 = 18, iun3 = 19, key, ier2
   integer :: ip1, ip2, ip3, ip1b, ip1c
   integer :: ig1, ig2, ig3, ig4
   integer :: kin, nx, ny, nz
   integer :: dateo, deet, npas, nz2, nbits, datyp, swa, lng, dltf,     &
                     ubc, extra1, extra2, extra3, nliste, nliste2, nliste3
   integer, dimension(nlev)  :: liste, Level_ip1
   integer, dimension(nlev2) :: liste2, Level_ip1b
   integer, dimension(nlev3) :: liste3, Level_ip1c

   integer fnom, fstouv, fstfrm, fstlir, fclos, fstinl, fstinf, fstprm

   integer :: kount

!___________________________________________________________
     open(unit=13, file='testgem414_0.dat', form='unformatted')
!     open(unit=14, file='teuv90d.dat', form='unformatted')

     direc = '../exp0/RUNMOD/output/' 
     direc2 = '../exp0/RUNMOD/output/'
      
!     surf_params = .true.
     vert_pro = .true.
     if (vert_pro) then
!        interp = .true.
!        calc_hgt = .true.
!        column   = .true.
     endif
!     soil_vars = .true.
     Schm_phys = .true.
!     Schm_chem = .true.
     dyna_outp = .true.
     if (Schm_chem) then
        chem_col = .true.
     endif 
     chem_col = .false.

     nlev_out = 102     !# of levels to write out

     sg7 = 0.0
! Main loop
   rpn_prefix = '2009100609-00-00_'  !'2009100609-00-00_'
   ct = 0
   ct1 = 0 ; tmstep = 1

   do ct0 = 2000, 2200, 2000
      write(cnum0, '(i10.10)') 0+5 !ct0 + 2000
      dir_prefix = trim(direc)//'laststep_'//trim(cnum0)
      dir_prefix2 = trim(direc2)//'laststep_'//trim(cnum0)

      do while (ct1 <= ct0+2000)
!         if (vert_pro) then
!          if (ct1 < 33416) cycle  !ls=90, 2100 ;ls=180, 19500 ;ls=270, 33300 ; 
!          if (ct1 > 33512) exit   !ls=90, 4300 ;ls=180, 21300 ;ls=270, 34900 ;
!         endif
!         if (ct1 < 26316) then
!            ct1 = ct1 + 6
!            cycle
!         endif
         if (ct1 > 5) exit
         write(cnum, '(i6.6)') ct1
!___________________________________________________
! flnom, the file name from which data are read. 
!___________________________________________________
        if (Schm_phys) then
          flnom = trim(dir_prefix)//'/000-000/pm'//trim(rpn_prefix)//trim(cnum)//'p'
          ier = fnom(iun, flnom, 'STD+RND', 0)
          if (ier < 0) print *, 'error while opening'
          ier = fstouv(iun, 'STD+RND')
        endif

        if (Schm_chem) then
          flnom2 = trim(dir_prefix)//'/km'//trim(rpn_prefix)//trim(cnum)//'p'
          ier = fnom(iun2, flnom2, 'STD+RND', 0)
          if (ier < 0) print *, 'error while opening'
          ier = fstouv(iun2, 'STD+RND')
        endif

        if (dyna_outp) then
          flnom3 = trim(dir_prefix2)//'/000-000/dm'//trim(rpn_prefix)//trim(cnum)//'p'
          ier = fnom(iun3, flnom3, 'STD+RND', 0)
          if (ier < 0) print *, 'error while opening'
          ier = fstouv(iun3, 'STD+RND')
        endif

        if (ct == 0) then
           if (vert_pro) then
             ier = fstinl(iun3, ni, nj, nk, -1, ' ', -1, -1, -1, 'P',          &
                           'TT', liste, nliste, nlev)
             if (nliste < nlev) then
                write(*,*) 'nlev must be equal to (or less than) nliste for   &
     &                   this to work properly, check-up', nlev, nliste, liste
                call exit
             endif
             do k = nlev, 1, -1 
                ier = fstprm(liste(k), dateo, deet, npas, ni, nj, nk, nbits,  &
                            datyp, ip1, ip2, ip3, typvar, nomvar, etiket,     &
                            grtyp, ig1, ig2, ig3, ig4, swa, lng, dltf,        &
                            ubc, extra1, extra2, extra3)
                Level_ip1(k) = ip1
                call convip(ip1, sig(k), kin, -1, etiket, .false.)
                sig2(k) = 1000.0 * sig(k)
             enddo
           elseif (Schm_phys) then
             ier = fstinf(iun, ni, nj, nk, -1, ' ', -1, -1, -1, 'P', 'J8')
             ier2 = fstprm(ier, dateo, deet, npas, ni, nj, nk, nbits,         &
                            datyp, ip1, ip2, ip3, typvar, nomvar, etiket,     &
                            grtyp, ig1, ig2, ig3, ig4, swa, lng, dltf,        &
                            ubc, extra1, extra2, extra3)
 
           endif
           kount = npas

! For the subsoil
           if (soil_vars) then
           ier = fstinl(iun, ni, nj, nk2, -1, ' ', -1, -1, -1, 'P',          &
                         'TL', liste2, nliste2, nlev2)
           do k = nlev2, 1, -1 
              ier = fstprm(liste2(k), dateo, deet, npas, ni, nj, nk2, nbits,  &
                          datyp, ip1b, ip2, ip3, typvar, nomvar, etiket,      &
                          grtyp, ig1, ig2, ig3, ig4, swa, lng, dltf,          &
                          ubc, extra1, extra2, extra3)
              Level_ip1b(k) = ip1b
           enddo
           endif

! Chemical species column abundance
           if (chem_col) then
           ier = fstinl(iun2, ni, nj, nk3, -1, ' ', -1, -1, -1, 'P',          &
                         'TOTP', liste3, nliste3, nlev3)
           do k = 1, nlev3
              ier = fstprm(liste3(k), dateo, deet, npas, ni, nj, nk3, nbits,  &
                          datyp, ip1c, ip2, ip3, typvar, nomvar, etiket,      &
                          grtyp, ig1, ig2, ig3, ig4, swa, lng, dltf,          &
                          ubc, extra1, extra2, extra3)
              Level_ip1c(k) = ip1c
           enddo
           endif
!*******

           ALLOCATE(ps(ni,nj), ts(ni,nj), cw(ni,nj), sn(ni,nj), lat(nj),      &
                    te(ni,nj,nlev), aq(ni,nj,nlev), pr(ni,nj,nlev), lon(ni),  &
                    uw(ni,nj,nlev), vw(ni,nj,nlev), ic(ni,nj,nlev),           &
                    ite(ni,nj,nlev), iaq(ni,nj,nlev), iic(ni,nj,nlev),        &
                    iuw(ni,nj,nlev), ivw(ni,nj,nlev), iop(ni,nj,nlev),        &
                    hgt(ni,nj,nlev), tsl(ni,nj,nlev2), rwv(ni,nj,nlev2),      &
                    rad(ni,nj,nlev2), ric(ni,nj,nlev2), rtw(ni,nj,nlev2),     &
                    hb(ni,nj), wo(ni,nj), tm(ni,nj), op(ni,nj,nlev),          &
                    itez(nj,nlev), iuwz(nj,nlev), iaqz(nj,nlev),              &
                    iicz(nj,nlev), iopz(nj,nlev), ivwz(nj,nlev),              &
                    phi(nj,nlev),  tep(ni,nj,nlev), tnd(ni,nj,nlev),          &
                    work(ni,nj),                                &
                    sn02(ni,nj), qs(ni,nj), qa(ni,nj), stat = ier2)
            if (Schm_chem) then
              ALLOCATE(o3v(ni,nj,nlev), o2v(ni,nj,nlev), o1d(ni,nj,nlev),     &
                    cov(ni,nj,nlev), aqc(ni,nj,nlev), h1v(ni,nj,nlev),        &
                    ohv(ni,nj,nlev), ho2v(ni,nj,nlev), h2o2v(ni,nj,nlev),     &
                    o3p(ni,nj,nlev), h2v(ni,nj,nlev),  o2dv(ni,nj,nlev),      &
                    tr_col(ni,nj,12), stat = ier2)
            endif

           if (interp) then
              paq => iaq ; pte => ite ;  pic => iic ; puw => iuw ; pvw => ivw
              pop => iop
           else
              paq => aq ; pte => te ;  pic => ic ; puw => uw ; pvw => vw
              pop => op
           endif

!           key = fstlir(lat, iun, nx, ny, nz, -1, etiket, -1, -1, -1, 'X', '^^')
!           key = fstlir(lon, iun, nx, ny, nz, -1, etiket, -1, -1, -1, 'X', '>>')
        endif

!        ip1 = Level_ip1(nlev)
!        key = fstlir(ps, iun3, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'P0')
        key = fstlir(ps, iun3, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'P0')
        key = fstlir(ts, iun, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'J8')

       if (surf_params) then
        key = fstlir(ps, iun, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'P8')

        key = fstlir(sn, iun, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'I5')

        key = fstlir(ts, iun, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'TG') 

        key = fstlir(cw, iun, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'IH')

        key = fstlir(wo, iun, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'W4')

        key = fstlir(qs, iun, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'J9')

!        key = fstlir(hb, iun, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'SN')

!        key = fstlir(hb, iun, nx, ny, nz, -1, etiket, -1, -1, -1, typvar, 'H')

!       key = fstlir(tm, iun, nx, ny, nz, -1, etiket, ip1, -1, -1, typvar, 'T8')
       endif

        if (vert_pro) then
         do k = nlev,1,-1
            ip1 = Level_ip1(k)

           key = fstlir(te(:,:,k), iun3, nx, ny, nz,                     &
                                   -1, etiket, ip1, -1, -1, typvar, 'TT')

!           key = fstlir(uw(:,:,k), iun3, nx, ny, nz,                     &
!                                   -1, etiket, ip1, -1, -1, typvar, 'UU')

!           key = fstlir(vw(:,:,k), iun3, nx, ny, nz,                     &
!                                   -1, etiket, ip1, -1, -1, typvar, 'VV')

           key = fstlir(aq(:,:,k), iun3, nx, ny, nz,                     &
                                   -1, etiket, ip1, -1, -1, typvar, 'HU')

!           key = fstlir(ic(:,:,k), iun3, nx, ny, nz,                     &
!                                   -1, etiket, ip1, -1, -1, typvar, 'GZ')

!           key = fstlir(op(:,:,k), iun3, nx, ny, nz,                     &
!                                   -1, etiket, ip1, -1, -1, typvar, 'NIT0')

           if (schm_chem) then

           key = fstlir(o3v(:,:,k), iun2, nx, ny, nz,                       &
                                   -1, etiket, ip1, -1, -1, typvar, 'O3')

           key = fstlir(o2v(:,:,k), iun2, nx, ny, nz,                       &
                                   -1, etiket, ip1, -1, -1, typvar, 'O2')

           key = fstlir(o1d(:,:,k), iun2, nx, ny, nz,                       &
                                   -1, etiket, ip1, -1, -1, typvar, 'O1D')

           key = fstlir(o3p(:,:,k), iun2, nx, ny, nz,                       &
                                   -1, etiket, ip1, -1, -1, typvar, 'O3P')

           key = fstlir(cov(:,:,k), iun2, nx, ny, nz,                       &
                                   -1, etiket, ip1, -1, -1, typvar, 'CO')

           key = fstlir(h1v(:,:,k), iun2, nx, ny, nz,                        &
                                   -1, etiket, ip1, -1, -1, typvar, 'H1')

           key = fstlir(h2v(:,:,k), iun2, nx, ny, nz,                       &
                                   -1, etiket, ip1, -1, -1, typvar, 'H2')

           key = fstlir(ohv(:,:,k), iun2, nx, ny, nz,                       &
                                   -1, etiket, ip1, -1, -1, typvar, 'OH')

           key = fstlir(ho2v(:,:,k), iun2, nx, ny, nz,                      &
                                   -1, etiket, ip1, -1, -1, typvar, 'HO2')

           key = fstlir(aqc(:,:,k), iun2, nx, ny, nz,                       &
                                   -1, etiket, ip1, -1, -1, typvar, 'H2O')

           key = fstlir(h2o2v(:,:,k), iun2, nx, ny, nz,                     &
                                   -1, etiket, ip1, -1, -1, typvar, 'H2O2')

           key = fstlir(o2dv(:,:,k), iun2, nx, ny, nz,                      &
                                   -1, etiket, ip1, -1, -1, typvar, 'O21D')
           endif
         enddo
        endif

        if (soil_vars) then
          do k = nlev2,1,-1
             ip1b = Level_ip1b(k)
             key = fstlir(tsl(:,:,k), iun, nx, ny, nz2,                      &
                                     -1, etiket, ip1b, -1, -1, typvar, 'TL')

             key = fstlir(rtw(:,:,k), iun, nx, ny, nz2,                      &
                                     -1, etiket, ip1b, -1, -1, typvar, 'J6')

             key = fstlir(rwv(:,:,k), iun, nx, ny, nz2,                      &
                                     -1, etiket, ip1b, -1, -1, typvar, 'I1')

             key = fstlir(rad(:,:,k), iun, nx, ny, nz2,                      &
                                     -1, etiket, ip1b, -1, -1, typvar, 'I2')

             key = fstlir(ric(:,:,k), iun, nx, ny, nz2,                      &
                                     -1, etiket, ip1b, -1, -1, typvar, 'I4')
          enddo
        endif

        if (chem_col) then
          do k = 1, nlev3
             ip1c = Level_ip1c(k)
             key = fstlir(tr_col(:,:,k), iun2, nx, ny, nz2,                 &
                                     -1, etiket, ip1c, -1, -1, typvar, 'TOTP')
          enddo
        endif 

!**********Interpolation*********************
        if (vert_pro) then

! calculate height
        if (calc_hgt) then
           do j = 1,nj
             do i = 1,ni
                hgt(i,j,nlev)   = 5.0
                hgt(i,j,nlev-1) = rgasd * (te(i,j,nlev-1))           &
                                   * alog(1.0/sig(nlev-1)) / grav
                do k = nlev-2,1,-1
                   hgt(i,j,k)   = hgt(i,j,k+1) + rgasd*0.5*(te(i,j,k)   &
                                 + te(i,j,k+1)) * alog(sig(k+1)/sig(k)) / grav
                enddo
             enddo
           enddo
        endif

        if (column) then
           do j = 1,nj
              do i = 1,ni
                 tm(i,j) = 0.0 ; sn(i,j) = 0.0 ; cw(i,j) = 0.0
                 hb(i,j) = 0.0 ; wo(i,j) = 0.0 ; ts(i,j) = 0.0
                 do k = 2,nlev
                   pr(i,j,k)  = ps(i,j)*sig(k)
                   tnd(i,j,k) = pr(i,j,k) / (bz * te(i,j,k))
                   dz         = (hgt(i,j,k-1) - hgt(i,j,k)) 

                   tm(i,j) = tm(i,j) + o3v(i,j,k) * tnd(i,j,k) * dz
                   sn(i,j) = sn(i,j) + o2v(i,j,k) * tnd(i,j,k) * dz
                   cw(i,j) = cw(i,j) + o1d(i,j,k) * tnd(i,j,k) * dz
                   hb(i,j) = hb(i,j) + aqc(i,j,k) * tnd(i,j,k) * dz
                   wo(i,j) = wo(i,j) + cov(i,j,k) * tnd(i,j,k) * dz
                   ts(i,j) = ts(i,j) + o2dv(i,j,k) * tnd(i,j,k) * dz
                 enddo
              enddo
           enddo
        endif

        if (interp) then
           do j = 1,nj
              do i = 1,ni
                 do k = 1,nlev
                   pr(i,j,k)  = ps(i,j)*sig(k)*100.0

!                   o3v(i,j,k) = o3v(i,j,k) * tnd(i,j,k)
!                   o2v(i,j,k) = o2v(i,j,k) * tnd(i,j,k)
!                   o1d(i,j,k) = o1d(i,j,k) * tnd(i,j,k)
!                   aqc(i,j,k) = aqc(i,j,k) * tnd(i,j,k)
!                   cov(i,j,k) = cov(i,j,k) * tnd(i,j,k)
!                   o2dv(i,j,k) = o2dv(i,j,k) * tnd(i,j,k)
                 enddo

!        if (interp) then
                 call spline1(pr(i,j,:), te(i,j,:), nlev, big_n, big_n, sg2)
!                 call spline1(pr(i,j,:), uw(i,j,:), nlev, big_n, big_n, sg3)
!                 call spline1(pr(i,j,:), vw(i,j,:), nlev, big_n, big_n, sg4)
!                 call spline1(pr(i,j,:), aq(i,j,:), nlev, big_n, big_n, sg5)
!                 call spline1(pr(i,j,:), ic(i,j,:), nlev, big_n, big_n, sg6)
!                 call spline1(pr(i,j,:), op(i,j,:), nlev, big_n, big_n, sg7)

                 call splint(pr(i,j,:), nlev, sig2, a, a2, b, b2, h2, klo, khi)

                 do k = 1,nlev
                    if (pr(i,j,nlev) > sig2(k)) then
                       l = klo(k)
                       m = khi(k)
                       ite(i,j,k) = a(k) * te(i,j,l) + b(k) * te(i,j,m)   &
                                    + (a2(k) * sg2(l) + b2(k) * sg2(m)) * h2(k)
!                       iuw(i,j,k) = a(k) * uw(i,j,l) + b(k) * uw(i,j,m)   &
!                                    + (a2(k) * sg3(l) + b2(k) * sg3(m)) * h2(k)
!                       ivw(i,j,k) = a(k) * vw(i,j,l) + b(k) * vw(i,j,m)   &
!                                    + (a2(k) * sg4(l) + b2(k) * sg4(m)) * h2(k)
!                       iaq(i,j,k) = a(k) * aq(i,j,l) + b(k) * aq(i,j,m)   &
!                                    + (a2(k) * sg5(l) + b2(k) * sg5(m)) * h2(k)
!                       iic(i,j,k) = a(k) * ic(i,j,l) + b(k) * ic(i,j,m)   &
!                                    + (a2(k) * sg6(l) + b2(k) * sg6(m)) * h2(k)
!                       iop(i,j,k) = a(k) * op(i,j,l) + b(k) * op(i,j,m)   &
!                                    + (a2(k) * sg7(l) + b2(k) * sg7(m)) * h2(k)
                    else
                       ite(i,j,k) = -999.0
!                       iuw(i,j,k) = -999.0
!                       ivw(i,j,k) = -999.0
!                       iaq(i,j,k) = -999.0
!                       iic(i,j,k) = -999.0
!                       iop(i,j,k) = -999.0
                    endif
                 enddo
!        endif

              enddo
           enddo

           do k = 1, nlev
              do j = 1, nj
                 ct3 = 0
                 itez(j,k) = 0.0
!                 iuwz(j,k) = 0.0
!                 ivwz(j,k) = 0.0
!                 iaqz(j,k) = 0.0
!                 iicz(j,k) = 0.0
!                 iopz(j,k) = 0.0
                 do i = 1, ni
                    if (pr(i,j,nlev) > sig2(k)) then
                       ct3 = ct3 + 1
                       itez(j,k) = itez(j,k) + te(i,j,k)
!                       iuwz(j,k) = iuwz(j,k) + uw(i,j,k)
!                       ivwz(j,k) = ivwz(j,k) + vw(i,j,k)
!                       iaqz(j,k) = iaqz(j,k) + aq(i,j,k)
!                       iicz(j,k) = iicz(j,k) + ic(i,j,k)
!                       iopz(j,k) = iopz(j,k) + op(i,j,k)
                    endif 
                 enddo
                 if (ct3 == 0) then
                   itez(j,k) = -999.0
!                   iuwz(j,k) = -999.0
!                   ivwz(j,k) = -999.0
!                   iaqz(j,k) = -999.0
!                   iicz(j,k) = -999.0
!                   iopz(j,k) = -999.0
                 else
                   itez(j,k) = itez(j,k) / ct3
!                   iuwz(j,k) = iuwz(j,k) / ct3
!                   ivwz(j,k) = ivwz(j,k) / ct3
!                   iaqz(j,k) = iaqz(j,k) / ct3
!                   iicz(j,k) = iicz(j,k) / ct3
!                   iopz(j,k) = iopz(j,k) / ct3
                 endif
              enddo
           enddo

! Calculate the streamfunction
!        call ccmp_zm_mspi (ni, nj, nlev, ivw, lat, sig2, ps, -999.0, phi)

        endif

        endif    !vertical profile

!********************
!        read(22)  ((work(i,j),i=1,ni),j=1,nj )
!        sn02 = sn / work

!* write output
        if (surf_params) then
           write(13) ((ts(i,j),i=1,ni),j=1,nj )
           write(13) ((ps(i,j),i=1,ni),j=1,nj )
!           write(13) ((sn02(i,j),i=1,ni),j=1,nj )
           write(13) ((sn(i,j),i=1,ni),j=1,nj )
           write(13) ((cw(i,j),i=1,ni),j=1,nj )
           write(13) ((wo(i,j),i=1,ni),j=1,nj )
           write(13) ((qs(i,j),i=1,ni),j=1,nj )
!           write(13) ((hb(i,j),i=1,ni),j=1,nj )
!           write(13) ((tm(i,j),i=1,ni),j=1,nj )
!            write(15) ((ts(i,j), i=1,ni),j=1,nj)
        endif

        ps = ps * 100.0
        te = te + 273.16
           write(13) ((ts(i,j),i=1,ni),j=1,nj )
           write(13) ((ps(i,j),i=1,ni),j=1,nj )
           do k = 1, nlev_out
              write(13) ((te(i,j,k),i=1,ni),j=1,nj )
           enddo
        if (vert_pro .and. schm_chem) then
           do k = 1, nlev_out
              write(13) ((te(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(13) ((uw(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(13) ((vw(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(13) ((aq(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(13) ((ic(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(13) ((op(i,j,k),i=1,ni),j=1,nj )
           enddo
        endif

        if (schm_chem) then
           do k = 1, nlev_out
             write(14) ((o3v(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((o2v(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((o1d(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((o3p(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((cov(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((h1v(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((h2v(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((ohv(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((ho2v(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((aqc(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((h2o2v(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1, nlev_out
              write(14) ((o2dv(i,j,k),i=1,ni),j=1,nj )
           enddo
        endif

        if (interp) then
           do k = nlev,1,-1
              write(14) (itez(j,k),j=1,nj )
           enddo
!           do k = nlev,1,-1
!              write(14) (iuwz(j,k),j=1,nj )
!           enddo
!           do k = nlev,1,-1
!              write(14) (ivwz(j,k),j=1,nj )
!           enddo
!           do k = nlev,1,-1
!              write(14) (iaqz(j,k),j=1,nj )
!           enddo
!           do k = nlev,1,-1
!              write(14) (iicz(j,k),j=1,nj )
!           enddo
!           do k = nlev,1,-1
!              write(14) (iopz(j,k),j=1,nj )
!           enddo
        endif    !vertical profile

        if (soil_vars) then
           do k = 1,nlev2
              write(13) ((tsl(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1,nlev2
              write(13) ((rtw(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1,nlev2
              write(13) ((rwv(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1,nlev2
              write(13) ((rad(i,j,k),i=1,ni),j=1,nj )
           enddo
           do k = 1,nlev2
              write(13) ((ric(i,j,k),i=1,ni),j=1,nj )
           enddo
        endif

        if (chem_col) then
           write(13) ((tr_col(i,j,1),i=1,ni),j=1,nj )
           write(13) ((tr_col(i,j,2),i=1,ni),j=1,nj )
           write(13) ((tr_col(i,j,5),i=1,ni),j=1,nj )
           write(13) ((tr_col(i,j,10),i=1,ni),j=1,nj )
           write(13) ((tr_col(i,j,12),i=1,ni),j=1,nj )
           write(13) ((tr_col(i,j,4),i=1,ni),j=1,nj )
!           do k = 1,nlev3
!              write(13) ((tr_col(i,j,k),i=1,ni),j=1,nj )
!           enddo
        endif


        if (Schm_phys) then
          ier = fstfrm(iun)
          ier = fclos(iun)
        endif
        if (Schm_chem) then
          ier = fstfrm(iun2)
          ier = fclos(iun2)
        endif
        if (dyna_outp) then
          ier = fstfrm(iun3)
          ier = fclos(iun3)
        endif

        ct = ct + 1
        ct1 = ct1 + tmstep
        enddo  
     enddo

!     close(13)
!     close(14)

     write(*,*) 'All done at ', ct

     DEALLOCATE(ps, ts, cw, sn, lat, lon, te, aq, pr, uw, vw, ic, op,     &
                    ite, iaq, iic, iuw, ivw, iop, itez, iuwz, iaqz,       &
                    iicz, iopz, ivwz,                                     &
                    phi, tnd, tep, hgt, tsl, rwv, rad, ric, rtw,          & 
                    hb, tm, work, sn02, wo, qs, qa, stat = ier2)
     if (Schm_chem) then
       DEALLOCATE(o3v, o2v, o1d, cov, aqc, h1v, ohv, ho2v, h2o2v,         &
                    o3p, h2v,  o2dv, tr_col, stat = ier2)
     endif

  stop
  end
