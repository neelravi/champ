module contrl
  !> Arguments: idump, irstar, isite, nconf, nblk, nblkeq, nconf_new, nstep

   integer  :: idump
   integer  :: irstar
   integer  :: isite
   integer  :: nconf
   integer  :: nblk
   integer  :: nblkeq
   integer  :: nconf_new
   integer  :: nstep

   private 
   public :: idump, irstar, isite, nconf, nblk, nblkeq, nconf_new, nstep 
   save
 end module contrl

 module contr2
   !> Arguments: i3body, ianalyt_lap, iaver, icusp, icusp2, ifock, ijas, irewgt, isc, istrch

    integer  :: i3body
    integer  :: ianalyt_lap
    integer  :: iaver
    integer  :: icusp
    integer  :: icusp2
    integer  :: ifock
    integer  :: ijas
    integer  :: irewgt
    integer  :: isc
    integer  :: istrch

    private
    public :: i3body, ianalyt_lap, iaver, icusp, icusp2, ifock, ijas, irewgt, isc, istrch
    save
 end module contr2

 module contr3
  !> Arguments: mode

  character*12 :: mode

  private
  public :: mode
  save
end module contr3

 module contrl_per
   !> Arguments: iperiodic, ibasis 

   integer  :: iperiodic, ibasis

   private
   public   :: iperiodic, ibasis
   save
 end module contrl_per
 
 module contrldmc
   !> Arguments: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq, itau_eff, nfprod, rttau, tau, taueff, tautot
   use force_mod, only: MFORCE
   use precision_kinds, only: dp

   integer  :: iacc_rej
   integer  :: icross
   integer  :: icuspg
   integer  :: icut_br
   integer  :: icut_e
   integer  :: idiv_v
   integer  :: idmc
   integer  :: ipq
   integer  :: itau_eff
   integer  :: nfprod
   real(dp) :: rttau
   real(dp) :: tau
   real(dp) :: taueff(MFORCE)
   real(dp) :: tautot

   private
   public :: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq, itau_eff, nfprod, rttau, tau, taueff, tautot
   save
 end module contrldmc