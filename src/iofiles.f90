module iofiles

implicit none
public

! define output file name, i/o unit number and output directory

! work character array for all output purposes
character (len=160) :: fname
!***tak 10-25-2012: add input file
! input file (in)
character(len=32), parameter :: input_file='input.dat'
integer, parameter :: iunit_input=7

! run debug file (out)
character(len=32), parameter :: res_file='output.dat'
integer, parameter :: iunit_res=4

! restart file (out)
character(len=32), parameter :: restart_file='start.dat'
integer, parameter :: iunit_restart=10
character(len=32), parameter :: restart_log_file='start.log'
integer, parameter :: iunit_restart_log=11

! max/min values of u,v,w,rho (out)
character(len=32), parameter :: minmax_file='minmaxtest.dat'
integer, parameter :: iunit_minmax=50

! cfl # data (out)
character(len=32), parameter ::cfl_file='cfl.dat'
integer, parameter :: iunit_cfl=51

!***tak 9-13-2012: uncomment - not used
! energy data (out)
character(len=32), parameter :: energy_file='energy.dat'
integer, parameter :: iunit_energy=60

! time step record (out)
character(len=32), parameter :: dtchange_file='dtchange.dat'
integer, parameter :: iunit_dtchange=65

! field solution output file list (out)
character(len=32), parameter ::sol_file='uvw.dat'
integer, parameter :: iunit_sol=68

! reference solution output file list (out)
character(len=32), parameter ::refsol_file='refsol.dat'
integer, parameter :: iunit_refsol=686

! vorticity field output file list (out)
character(len=32), parameter ::vort_file='vort.dat'
integer, parameter :: iunit_vort=619

! parallel data (in/out)
character(len=32), parameter :: paral_file='paral.dat'
integer, parameter :: iunit_paral=679

! grid data (in/out)
character(len=32), parameter :: grid_file='grid.dat'
integer, parameter :: iunit_grid=680

! x-grid data (out)
character(len=32), parameter :: xgrid_file='xgrid.dat'
integer, parameter :: iunit_xgrid=681

! y-grid data (out)
character(len=32), parameter :: ygrid_file='ygrid.dat'
integer, parameter :: iunit_ygrid=682

! z-grid data (out)
character(len=32), parameter :: zgrid_file='zgrid.dat'
integer, parameter :: iunit_zgrid=683

! probe data file name (out)
character(len=30), parameter :: prb_out_1='probe01.dat'
character(len=30), parameter :: prb_out_2='probe02.0dat'
character(len=30), parameter :: prb_out_3='probe03.dat'
character(len=30), parameter :: prb_out_4='probe04.dat'

end module iofiles