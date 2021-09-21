program main
   
   ! Read original high-res data
   ! Read coarse-grained high-res data
   ! Calculate air-sea fluxes for both using Aerobulk
   !
   ! Author: Joakim Kjellsson, Sept 2021
   ! Notes: Larry...Larry...Larry...you dummy
   !
   
   use netcdf
   
   ! Load Aerobulk from Laurent Brodeau
   ! https://github.com/brodeau/aerobulk
   use mod_aerobulk
   
   implicit none
   
   integer, parameter :: sp = selected_real_kind(6, 37)
   integer, parameter :: dp = selected_real_kind(15, 307)
   
   integer, parameter :: nx_hires = 1024
   integer, parameter :: ny_hires = 768
   integer, parameter :: nt_hires = 360
   integer, parameter :: nl_hires = 1
   
   integer, parameter :: nx_lores = 192
   integer, parameter :: ny_lores = 96
   integer, parameter :: nt_lores = 360
   integer, parameter :: nl_lores = 1
   
   integer :: ncid, x_dimid, y_dimid, t_dimid, &
               & lon_varid, lat_varid, time_varid, &
               & qe_varid, qh_varid, taux_varid,  &
               & tauy_varid, evap_varid, dimids(3), &
               & start(3), count(3)
   
   integer(kind=sp) :: irec
   
   real(kind=sp), dimension(nx_hires) :: lon_hires
   real(kind=sp), dimension(ny_hires) :: lat_hires
   real(kind=sp), dimension(nl_hires) :: time_hires
   
   real(kind=sp), dimension(nx_lores) :: lon_lores
   real(kind=sp), dimension(ny_lores) :: lat_lores
   real(kind=sp), dimension(nl_lores) :: time_lores
   
   real(kind=sp), dimension(nx_hires, ny_hires, nl_hires) :: uas_hires
   real(kind=sp), dimension(nx_hires, ny_hires, nl_hires) :: vas_hires
   real(kind=sp), dimension(nx_hires, ny_hires, nl_hires) :: tas_hires
   real(kind=sp), dimension(nx_hires, ny_hires, nl_hires) :: huss_hires
   real(kind=sp), dimension(nx_hires, ny_hires, nl_hires) :: psl_hires
   real(kind=sp), dimension(nx_hires, ny_hires, nl_hires) :: sst_hires
   
   real(kind=sp), dimension(nx_lores, ny_lores, nl_lores) :: uas_lores
   real(kind=sp), dimension(nx_lores, ny_lores, nl_lores) :: vas_lores
   real(kind=sp), dimension(nx_lores, ny_lores, nl_lores) :: tas_lores
   real(kind=sp), dimension(nx_lores, ny_lores, nl_lores) :: huss_lores
   real(kind=sp), dimension(nx_lores, ny_lores, nl_lores) :: psl_lores
   real(kind=sp), dimension(nx_lores, ny_lores, nl_lores) :: sst_lores
   
   integer(kind=sp) :: year, month, day
   integer(kind=sp) :: jt, jl
   
   real(kind=sp), allocatable, dimension(:,:,:) :: Qe_hires
   real(kind=sp), allocatable, dimension(:,:,:) :: Qh_hires
   real(kind=sp), allocatable, dimension(:,:,:) :: Tau_x_hires
   real(kind=sp), allocatable, dimension(:,:,:) :: Tau_y_hires
   real(kind=sp), allocatable, dimension(:,:,:) :: Evap_hires
   
   real(kind=sp), allocatable, dimension(:,:,:) :: Qe_lores
   real(kind=sp), allocatable, dimension(:,:,:) :: Qh_lores
   real(kind=sp), allocatable, dimension(:,:,:) :: Tau_x_lores
   real(kind=sp), allocatable, dimension(:,:,:) :: Tau_y_lores
   real(kind=sp), allocatable, dimension(:,:,:) :: Evap_lores
   
   real(kind=dp), allocatable, dimension(:,:) :: sst
   real(kind=dp), allocatable, dimension(:,:) :: t_zt
   real(kind=dp), allocatable, dimension(:,:) :: hum_zt
   real(kind=dp), allocatable, dimension(:,:) :: U_zu
   real(kind=dp), allocatable, dimension(:,:) :: V_zu
   real(kind=dp), allocatable, dimension(:,:) :: SLP
   real(kind=dp), allocatable, dimension(:,:) :: Qe
   real(kind=dp), allocatable, dimension(:,:) :: Qh
   real(kind=dp), allocatable, dimension(:,:) :: Tau_x
   real(kind=dp), allocatable, dimension(:,:) :: Tau_y
   real(kind=dp), allocatable, dimension(:,:) :: Evap
   real(kind=dp), allocatable, dimension(:,:) :: mask
   
   character(len=200) :: indir, outdir, indir_lores
   
   character(len=10), parameter :: x_name = "lon"
   character(len=10), parameter :: y_name = "lat"
   character(len=10), parameter :: t_name = "time"
   character(len=200) :: outfile
   character(len=50)  :: start_date, time_units
   
   logical :: lhires = .true.
   logical :: llores = .true.
   
   indir = "/Users/joakimkjellsson/Downloads/highresmip/hadgem-hh/1day/"
   indir_lores = "/Users/joakimkjellsson/Downloads/highresmip/hadgem-hh/1day/crs_N48/"
   
   if (lhires) then
      
      year = 1950
      month = 1
      day = 1
      
      print*," --- STARTING high-res calculations --- "
      
      allocate ( Qe_hires     (nx_hires, ny_hires, nl_hires),  &
                 Qh_hires     (nx_hires, ny_hires, nl_hires),  & 
                 Tau_x_hires  (nx_hires, ny_hires, nl_hires),  &
                 Tau_y_hires  (nx_hires, ny_hires, nl_hires),  &
                 Evap_hires   (nx_hires, ny_hires, nl_hires) )
         
      allocate ( sst    (nx_hires, ny_hires), &
                 t_zt   (nx_hires, ny_hires), & 
                 hum_zt (nx_hires, ny_hires), &
                 U_zu   (nx_hires, ny_hires), &
                 V_zu   (nx_hires, ny_hires), &
                 SLP    (nx_hires, ny_hires), &
                 Qe     (nx_hires, ny_hires), &
                 Qh     (nx_hires, ny_hires), &
                 Tau_x  (nx_hires, ny_hires), &
                 Tau_y  (nx_hires, ny_hires), &
                 Evap   (nx_hires, ny_hires), &
                 mask   (nx_hires, ny_hires) )
      
      
      print*," --- OPENING output file --- "
      outfile = "fluxes_original.nc"
      call check( nf90_create( trim(outfile), nf90_clobber, ncid ) )
      
      call check( nf90_def_dim( ncid, x_name, nx_hires, x_dimid ) )
      call check( nf90_def_dim( ncid, y_name, ny_hires, y_dimid ) )
      call check( nf90_def_dim( ncid, t_name, nf90_unlimited, t_dimid ) )
      
      call check( nf90_def_var( ncid, x_name, nf90_real, x_dimid, lon_varid ) )
      call check( nf90_def_var( ncid, y_name, nf90_real, y_dimid, lat_varid ) )
      call check( nf90_def_var( ncid, t_name, nf90_real, t_dimid, time_varid ) )
      
      write(start_date,"(I0.4,A1,I0.2,A1,I0.2)") year,"-",month,"-",day
      time_units = "days since "//trim(start_date)
      print*,time_units
      call check( nf90_put_att(ncid, lon_varid, "units", "degrees_east" ) )
      call check( nf90_put_att(ncid, lat_varid, "units", "degrees_north" ) )
      call check( nf90_put_att(ncid, time_varid, "units", trim(time_units) ) )
      call check( nf90_put_att(ncid, time_varid, "calendar", "360_day" ) )
      call check( nf90_put_att(ncid, lon_varid, "standard_name", "longitude" ) )
      call check( nf90_put_att(ncid, lat_varid, "standard_name", "latitude" ) )
      call check( nf90_put_att(ncid, time_varid, "standard_name", "time" ) )
      
      dimids = (/ x_dimid, y_dimid, t_dimid /)
      
      call check( nf90_def_var( ncid, "Qe",    nf90_real, dimids, qe_varid ) )
      call check( nf90_def_var( ncid, "Qh",    nf90_real, dimids, qh_varid ) )
      call check( nf90_def_var( ncid, "Tau_x", nf90_real, dimids, taux_varid ) )
      call check( nf90_def_var( ncid, "Tau_y", nf90_real, dimids, tauy_varid ) )
      call check( nf90_def_var( ncid, "Evap",  nf90_real, dimids, evap_varid ) )
      
      call check( nf90_put_att( ncid, qe_varid,   "units", "W/m2" ) )
      call check( nf90_put_att( ncid, qh_varid,   "units", "W/m2" ) )
      call check( nf90_put_att( ncid, taux_varid, "units", "N/m2" ) )
      call check( nf90_put_att( ncid, tauy_varid, "units", "N/m2" ) )
      call check( nf90_put_att( ncid, evap_varid, "units", "mm/s" ) )
      
      call check( nf90_put_att( ncid, qe_varid,   "standard_name", "surface_upward_latent_heat_flux" ) )
      call check( nf90_put_att( ncid, qh_varid,   "standard_name", "surface_upward_sensible_heat_flux" ) )
      call check( nf90_put_att( ncid, taux_varid, "standard_name", "surface_downward_eastward_stress" ) )
      call check( nf90_put_att( ncid, tauy_varid, "standard_name", "surface_downward_northward_stress" ) )
      call check( nf90_put_att( ncid, evap_varid, "standard_name", "evaporation" ) )
      
      call check( nf90_enddef(ncid) )
      
      do jt = 1, nt_hires
         
         jl = 1
         irec = jt
         
         print*," --- READING high-res data --- "
         call read_hires(indir, year, month, day, jt)
         
         if (jt == 1) then
            ! write lon, lat on first time step
            call check( nf90_put_var( ncid, lon_varid, lon_hires ) )
            call check( nf90_put_var( ncid, lat_varid, lat_hires ) )
            
         end if
         
         print*," --- CALL AEROBULK for high-res data, step = ",jt," -- "
         ! jt = current time step
         ! nt = number of time steps
         ! algo = algorithm to use
         ! zt = height for temp, hum
         ! zu = height for wind
         ! sst = SST
         ! t_zt = air temp at zt
         ! hum_zt = air temp at zt
         ! U_zu = zonal wind at zu
         ! V_zu = meridional wind at zu
         ! slp = sea-lev pressure
         
         sst(:,:)    = sst_hires(:,:,jl) + 273.15
         t_zt(:,:)   = tas_hires(:,:,jl)
         hum_zt(:,:) = huss_hires(:,:,jl)
         U_zu(:,:)   = uas_hires(:,:,jl)
         V_zu(:,:)   = vas_hires(:,:,jl)
         SLP(:,:)    = psl_hires(:,:,jl)
         
         call aerobulk_model( 1, 1, 'coare3p0', 2.d0, 10.d0, sst, t_zt, hum_zt,   &
           &                  U_zu, V_zu, SLP, Qe, Qh,          & 
           &                  Tau_x, Tau_y, Evap, Niter=10)
         
         ! Qe = lat heat flux
         ! Qh = sens heat flux
         ! tau_x = zonal wind str
         ! tau_y = meridional wind str
         ! evap = evaporation 
         
         ! find land-sea mask using sst
         where (sst > 333.0 .OR. sst < 263.0)  
            mask = 0.0
         elsewhere
            mask = 1.0
         end where
         
         Qe(:,:) = Qe(:,:) * mask(:,:)
         Qh(:,:) = Qh(:,:) * mask(:,:)
         Tau_x(:,:) = Tau_x(:,:) * mask(:,:)
         Tau_y(:,:) = Tau_y(:,:) * mask(:,:)
         Evap(:,:) = Evap(:,:) * mask(:,:)
         
         count = (/ nx_hires, ny_hires, 1 /)
         start = (/        1,        1, 1 /)
         start(3) = irec
         
         print*," --- WRITE high-res data --- "
         call check( nf90_put_var( ncid, time_varid, time_hires, start=(/irec/), count=(/1/) ) )
         call check( nf90_put_var( ncid, qe_varid,   Qe,    start = start, count = count) )
         call check( nf90_put_var( ncid, qh_varid,   Qh,    start = start, count = count) )
         call check( nf90_put_var( ncid, taux_varid, Tau_x, start = start, count = count) )
         call check( nf90_put_var( ncid, tauy_varid, Tau_y, start = start, count = count) )
         call check( nf90_put_var( ncid, evap_varid, Evap,  start = start, count = count) )
         
      
      end do
      
      call check( nf90_close(ncid) )
      
      deallocate ( Qe_hires ,Qh_hires, Tau_x_hires, Tau_y_hires, Evap_hires )
     
      deallocate ( sst, t_zt, hum_zt, U_zu, V_zu, SLP, Qe, Qh, Tau_x, Tau_y, Evap, mask )
      
      print*," --- SUCCESS for high-res data --- "
      
   end if 
   
   
   if (llores) then
   
      print*," --- STARTING low-res calculations --- "
      
      allocate ( Qe_lores     (nx_lores, ny_lores, nl_lores),  &
                 Qh_lores     (nx_lores, ny_lores, nl_lores),  & 
                 Tau_x_lores  (nx_lores, ny_lores, nl_lores),  &
                 Tau_y_lores  (nx_lores, ny_lores, nl_lores),  &
                 Evap_lores   (nx_lores, ny_lores, nl_lores) )
     
      allocate ( sst    (nx_lores, ny_lores), &
                 t_zt   (nx_lores, ny_lores), & 
                 hum_zt (nx_lores, ny_lores), &
                 U_zu   (nx_lores, ny_lores), &
                 V_zu   (nx_lores, ny_lores), &
                 SLP    (nx_lores, ny_lores), &
                 Qe     (nx_lores, ny_lores), &
                 Qh     (nx_lores, ny_lores), &
                 Tau_x  (nx_lores, ny_lores), &
                 Tau_y  (nx_lores, ny_lores), &
                 Evap   (nx_lores, ny_lores), &
                 mask   (nx_lores, ny_lores) )
      
      
      print*," --- OPENING output file --- "
      outfile = "fluxes_lowres.nc"
      call check( nf90_create( trim(outfile), nf90_clobber, ncid ) )
      
      call check( nf90_def_dim( ncid, x_name, nx_lores, x_dimid ) )
      call check( nf90_def_dim( ncid, y_name, ny_lores, y_dimid ) )
      call check( nf90_def_dim( ncid, t_name, nf90_unlimited, t_dimid ) )
      
      call check( nf90_def_var( ncid, x_name, nf90_real, x_dimid, lon_varid ) )
      call check( nf90_def_var( ncid, y_name, nf90_real, y_dimid, lat_varid ) )
      call check( nf90_def_var( ncid, t_name, nf90_real, t_dimid, time_varid ) )
      
      write(start_date,"(I0.4,A1,I0.2,A1,I0.2)") year,"-",month,"-",day
      time_units = "days since "//trim(start_date)
      print*,time_units
      call check( nf90_put_att(ncid, lon_varid, "units", "degrees_east" ) )
      call check( nf90_put_att(ncid, lat_varid, "units", "degrees_north" ) )
      call check( nf90_put_att(ncid, time_varid, "units", time_units ) )
      call check( nf90_put_att(ncid, time_varid, "calendar", "360_day" ) )
      call check( nf90_put_att(ncid, lon_varid, "standard_name", "longitude" ) )
      call check( nf90_put_att(ncid, lat_varid, "standard_name", "latitude" ) )
      call check( nf90_put_att(ncid, time_varid, "standard_name", "time" ) )
      
      dimids = (/ x_dimid, y_dimid, t_dimid /)
      
      call check( nf90_def_var( ncid, "Qe",    nf90_real, dimids, qe_varid ) )
      call check( nf90_def_var( ncid, "Qh",    nf90_real, dimids, qh_varid ) )
      call check( nf90_def_var( ncid, "Tau_x", nf90_real, dimids, taux_varid ) )
      call check( nf90_def_var( ncid, "Tau_y", nf90_real, dimids, tauy_varid ) )
      call check( nf90_def_var( ncid, "Evap",  nf90_real, dimids, evap_varid ) )
      
      call check( nf90_put_att( ncid, qe_varid,   "units", "W/m2" ) )
      call check( nf90_put_att( ncid, qh_varid,   "units", "W/m2" ) )
      call check( nf90_put_att( ncid, taux_varid, "units", "N/m2" ) )
      call check( nf90_put_att( ncid, tauy_varid, "units", "N/m2" ) )
      call check( nf90_put_att( ncid, evap_varid, "units", "mm/s" ) )
      
      call check( nf90_put_att( ncid, qe_varid,   "standard_name", "surface_upward_latent_heat_flux" ) )
      call check( nf90_put_att( ncid, qh_varid,   "standard_name", "surface_upward_sensible_heat_flux" ) )
      call check( nf90_put_att( ncid, taux_varid, "standard_name", "surface_downward_eastward_stress" ) )
      call check( nf90_put_att( ncid, tauy_varid, "standard_name", "surface_downward_northward_stress" ) )
      call check( nf90_put_att( ncid, evap_varid, "standard_name", "evaporation" ) )
      
      call check( nf90_enddef(ncid) )
      
      jl = 1
      
      year = 1950
      month = 1
      day = 1
      
      do jt = 1, nt_lores
   
         print*," --- READING low-res data --- "
         call read_lores(year, month, day, jt)
         
         if (jt == 1) then
         
            call check( nf90_put_var( ncid, lon_varid, lon_lores ) )
            call check( nf90_put_var( ncid, lat_varid, lat_lores ) )
            
         end if
            
         
         print*," --- CALL AEROBULK for low-res data, step = ",jt," -- "
         ! jt = current time step
         ! nt = number of time steps
         ! algo = algorithm to use
         ! zt = height for temp, hum
         ! zu = height for wind
         ! sst = SST
         ! t_zt = air temp at zt
         ! hum_zt = air temp at zt
         ! U_zu = zonal wind at zu
         ! V_zu = meridional wind at zu
         ! slp = sea-lev pressure
         
         sst(:,:)    = sst_lores(:,:,jl) + 273.15
         t_zt(:,:)   = tas_lores(:,:,jl)
         hum_zt(:,:) = huss_lores(:,:,jl)
         U_zu(:,:)   = uas_lores(:,:,jl)
         V_zu(:,:)   = vas_lores(:,:,jl)
         SLP(:,:)    = psl_lores(:,:,jl)
         
         call aerobulk_model( 1, 1, 'coare3p0', 2.d0, 10.d0, sst, t_zt, hum_zt,   &
           &                  U_zu, V_zu, SLP, Qe, Qh,          & 
           &                  Tau_x, Tau_y, Evap, Niter=10)
         
         ! Qe = lat heat flux
         ! Qh = sens heat flux
         ! tau_x = zonal wind str
         ! tau_y = meridional wind str
         ! evap = evaporation 
         
         ! find land-sea mask using sst
         where (sst > 333.0 .OR. sst < 263.0)  
            mask = 0.0
         elsewhere
            mask = 1.0
         end where
         
         Qe(:,:) = Qe(:,:) * mask(:,:)
         Qh(:,:) = Qh(:,:) * mask(:,:)
         Tau_x(:,:) = Tau_x(:,:) * mask(:,:)
         Tau_y(:,:) = Tau_y(:,:) * mask(:,:)
         Evap(:,:) = Evap(:,:) * mask(:,:)
         
         print*," --- WRITE low-res data --- "
         irec = jt
         count = (/ nx_lores, ny_lores, 1 /)
         start = (/        1,        1, 1 /)
         start(3) = irec
         call check( nf90_put_var( ncid, time_varid, time_lores, start=(/irec/), count=(/1/) ) )
         call check( nf90_put_var( ncid, qe_varid,   Qe,    start = start, count = count) )
         call check( nf90_put_var( ncid, qh_varid,   Qh,    start = start, count = count) )
         call check( nf90_put_var( ncid, taux_varid, Tau_x, start = start, count = count) )
         call check( nf90_put_var( ncid, tauy_varid, Tau_y, start = start, count = count) )
         call check( nf90_put_var( ncid, evap_varid, Evap,  start = start, count = count) )
         
      end do
      
      call check( nf90_close(ncid) )
      
      deallocate ( Qe_lores ,Qh_lores, Tau_x_lores, Tau_y_lores, Evap_lores )
     
      deallocate ( sst, t_zt, hum_zt, U_zu, V_zu, SLP, Qe, Qh, Tau_x, Tau_y, Evap, mask )
      
      print*," --- SUCCESS for low-res data --- "
   
   end if 
   
   

   contains
   
   
   subroutine read_field(filename, field_name, nx, ny, nt, lon, lat, time, irec, field)
      
      ! take filename as input
      character (len=*), intent(in) :: filename
      
      character (len=*), intent(in) :: field_name
      
      integer(kind=sp), intent(in) :: nx, ny, nt, irec
      
      ! return grid lon, lat, time
      real(kind=sp), dimension(nx), intent(out) :: lon
      real(kind=sp), dimension(ny), intent(out) :: lat
      real(kind=sp), dimension(nt), intent(out) :: time
      ! return u10, v10, t2m, q2m as output
      ! sizes are set in the main routine
      real(kind=sp), dimension(nx,ny,nt), intent(out) :: field
      
      character (len=*), parameter :: lon_name = "lon"
      character (len=*), parameter :: lat_name = "lat"
      character (len=*), parameter :: time_name = "time"
      
      integer(kind=sp) :: ncid
      integer(kind=sp) :: lon_varid, lat_varid, time_varid, field_varid
      integer(kind=sp) :: start(3), count(3), jrec
      
      ! open the file in non-writable mode
      call check( nf90_open(filename, nf90_nowrite, ncid) )
      
      ! get ids for lon/lat coordinates and time
      call check( nf90_inq_varid(ncid, lat_name, lat_varid) )
      call check( nf90_inq_varid(ncid, lon_name, lon_varid) )
      call check( nf90_inq_varid(ncid, time_name, time_varid) )
      
      ! get id for field
      call check( nf90_inq_varid(ncid, field_name, field_varid) )
      
      ! read coordinates
      call check( nf90_get_var(ncid, lon_varid, lon, start=(/1/), count=(/nx/)) )
      call check( nf90_get_var(ncid, lat_varid, lat, start=(/1/), count=(/ny/)) )
      call check( nf90_get_var(ncid, time_varid, time, start=(/irec/), count=(/1/)) )
      
      ! read field
      start = (/ 1, 1, 1 /)
      count = (/ nx, ny, 1/)
      
      !do jrec = 1, nt
         
      start(3) = irec
      call check( nf90_get_var(ncid, field_varid, field(:,:,1), start=start, count=count) )
      print*,'done reading'
      
      !enddo 
      
      ! close file
      call check( nf90_close(ncid) )
      
      ! reassume the user 
      print*, " --- SUCCESS in reading field ",field_name," at step ",irec," from ",filename
      
      
   end subroutine read_field
   
   subroutine read_hires(indir, year, month, day, step)
      
      integer, intent(in) :: year, month, day, step
      character(len=200) :: indir
      
      character(len=200) :: ufile
      character(len=200) :: vfile
      character(len=200) :: pfile
      character(len=200) :: tfile
      character(len=200) :: qfile
      character(len=200) :: sfile
      
      ufile = trim(indir)//"uas_day_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230.nc"
      vfile = trim(indir)//"vas_day_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230.nc"
      tfile = trim(indir)//"tas_day_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230.nc"
      qfile = trim(indir)//"huss_day_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230.nc"
      pfile = trim(indir)//"psl_day_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230.nc"
      sfile = trim(indir)//"tos_Oday_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230_tasgrid.nc"
      
      call read_field(ufile, "uas", nx_hires, ny_hires, nt_hires, lon_hires, lat_hires, time_hires, step, uas_hires)
      call read_field(vfile, "vas", nx_hires, ny_hires, nt_hires, lon_hires, lat_hires, time_hires, step, vas_hires)
      call read_field(tfile, "tas", nx_hires, ny_hires, nt_hires, lon_hires, lat_hires, time_hires, step, tas_hires)
      call read_field(qfile, "huss", nx_hires, ny_hires, nt_hires, lon_hires, lat_hires, time_hires, step, huss_hires)
      call read_field(pfile, "psl", nx_hires, ny_hires, nt_hires, lon_hires, lat_hires, time_hires, step, psl_hires)
      call read_field(sfile, "tos", nx_hires, ny_hires, nt_hires, lon_hires, lat_hires, time_hires, step, sst_hires)
      
   end subroutine read_hires
   
   subroutine read_lores(year, month, day, step)
      
      integer, intent(in) :: year, month, day, step
      
      character(len=200) :: ufile
      character(len=200) :: vfile
      character(len=200) :: pfile
      character(len=200) :: tfile
      character(len=200) :: qfile
      character(len=200) :: sfile
      
      ufile = trim(indir_lores)//"uas_day_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230.nc"
      vfile = trim(indir_lores)//"vas_day_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230.nc"
      tfile = trim(indir_lores)//"tas_day_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230.nc"
      qfile = trim(indir_lores)//"huss_day_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230.nc"
      pfile = trim(indir_lores)//"psl_day_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230.nc"
      sfile = trim(indir_lores)//"tos_Oday_HadGEM3-GC31-HH_control-1950_r1i1p1f1_gn_19500101-19501230_tasgrid.nc"
      
      call read_field(ufile, "uas", nx_lores, ny_lores, nt_lores, lon_lores, lat_lores, time_lores, step, uas_lores)
      call read_field(vfile, "vas", nx_lores, ny_lores, nt_lores, lon_lores, lat_lores, time_lores, step, vas_lores)
      call read_field(tfile, "tas", nx_lores, ny_lores, nt_lores, lon_lores, lat_lores, time_lores, step, tas_lores)
      call read_field(qfile, "huss", nx_lores, ny_lores, nt_lores, lon_lores, lat_lores, time_lores, step, huss_lores)
      call read_field(pfile, "psl", nx_lores, ny_lores, nt_lores, lon_lores, lat_lores, time_lores, step, psl_lores)
      call read_field(sfile, "tos", nx_lores, ny_lores, nt_lores, lon_lores, lat_lores, time_lores, step, sst_lores)
      
   end subroutine read_lores
   
   
   subroutine check(status)
      integer, intent ( in) :: status
      
      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if
   end subroutine check  


end program main
