module move_walker
      use age,     only: iage
      use system,    only: ncent
      use branch,  only: eold,nwalk,pwt,wt,wthist
      use config,  only: d2o,peo_dmc,psido_dmc,psijo_dmc,vold_dmc
      use config,  only: xold_dmc
      use force_pth, only: PTH
      use m_force_analytic, only: iforce_analy
      use jacobsave, only: ajacold
      use mmpol_reduce_mod, only: mmpol_recv,mmpol_send
      use mpi
      use multiple_geo, only: nforce,nwprod
      use pcm_reduce_mod, only: pcm_recv,pcm_send
      use prop_reduce_mod, only: prop_recv,prop_send
      use system,  only: nelec
      use vd_mod, only: ehist, esnake, deriv_eold, dmc_ivd
      use velratio, only: fratio
      use pathak_mod, only: pold
      contains
      subroutine send_walker(irecv)
! Written by Claudia Filippi

      implicit none

      integer :: ierr, ifr, ip, irecv, irequest, iph
      integer :: isend, itag
      integer, dimension(MPI_STATUS_SIZE) :: istatus




      call mpi_isend(wt(nwalk),1,mpi_double_precision,irecv,1 &
      ,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(iage(nwalk),1,mpi_integer,irecv,2 &
      ,MPI_COMM_WORLD,irequest,ierr)

      itag=2
      do ifr=1,nforce
        call mpi_isend(ajacold(nwalk,ifr),1,mpi_double_precision,irecv &
        ,itag+1,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(eold(nwalk,ifr),1,mpi_double_precision,irecv &
        ,itag+2,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(psido_dmc(nwalk,ifr),1,mpi_double_precision,irecv &
        ,itag+3,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(psijo_dmc(nwalk,ifr),1,mpi_double_precision,irecv &
        ,itag+4,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(peo_dmc(nwalk,ifr),1,mpi_double_precision,irecv &
        ,itag+5,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(d2o(nwalk,ifr),1,mpi_double_precision,irecv &
        ,itag+6,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(pwt(nwalk,ifr),1,mpi_double_precision,irecv &
        ,itag+7,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(fratio(nwalk,ifr),1,mpi_double_precision,irecv &
        ,itag+8,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(vold_dmc(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
        ,irecv,itag+9,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(xold_dmc(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
        ,irecv,itag+10,MPI_COMM_WORLD,irequest,ierr)
        itag=itag+10
        do ip=0,nwprod-1
        itag=itag+1
        call mpi_isend(wthist(nwalk,ip,ifr),1,mpi_double_precision,irecv &
        ,itag,MPI_COMM_WORLD,irequest,ierr)
        enddo
      enddo

!     call send_det(itag,irecv)
!     call send_jas(itag,irecv)

!     nwalk=nwalk-1

      if(iforce_analy.eq.1) then
        if(dmc_ivd.gt.0) then
          itag=itag+1
          call mpi_isend(pold(nwalk,1),PTH,mpi_double_precision,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
          itag=itag+1
          call mpi_isend(deriv_eold(1,1,nwalk),3*ncent,mpi_double_precision,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
          itag=itag+1
          call mpi_isend(esnake(1,1,nwalk,1),3*ncent*PTH,mpi_double_precision,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
          do  ip=0,nwprod-1
            itag=itag+1
            call mpi_isend(ehist(1,1,nwalk,ip,1),3*ncent*PTH,mpi_double_precision,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
          enddo
        else
          itag=itag+1
          call mpi_isend(pold(nwalk,1),PTH,mpi_double_precision,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
        endif
      endif

      call prop_send(irecv,itag)
      call pcm_send(irecv,itag)
      call mmpol_send(irecv,itag)

      end subroutine

      subroutine recv_walker(isend)
      implicit none

      integer :: ierr, ifr, ip, irecv, irequest
      integer :: isend, itag
      integer, dimension(MPI_STATUS_SIZE) :: istatus

!     nwalk=nwalk+1

      call mpi_recv(wt(nwalk),1,mpi_double_precision,isend,1 &
      ,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(iage(nwalk),1,mpi_integer,isend,2 &
      ,MPI_COMM_WORLD,istatus,ierr)

      itag=2
      do ifr=1,nforce
        call mpi_recv(ajacold(nwalk,ifr),1,mpi_double_precision,isend &
        ,itag+1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(eold(nwalk,ifr),1,mpi_double_precision,isend &
        ,itag+2,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(psido_dmc(nwalk,ifr),1,mpi_double_precision,isend &
        ,itag+3,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(psijo_dmc(nwalk,ifr),1,mpi_double_precision,isend &
        ,itag+4,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(peo_dmc(nwalk,ifr),1,mpi_double_precision,isend &
        ,itag+5,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(d2o(nwalk,ifr),1,mpi_double_precision,isend &
        ,itag+6,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(pwt(nwalk,ifr),1,mpi_double_precision,isend &
        ,itag+7,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(fratio(nwalk,ifr),1,mpi_double_precision,isend &
        ,itag+8,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(vold_dmc(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
        ,isend,itag+9,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(xold_dmc(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
        ,isend,itag+10,MPI_COMM_WORLD,istatus,ierr)
        itag=itag+10
        do ip=0,nwprod-1
        itag=itag+1
        call mpi_recv(wthist(nwalk,ip,ifr),1,mpi_double_precision,isend &
        ,itag,MPI_COMM_WORLD,istatus,ierr)
        enddo
      enddo

!     call recv_det(itag,isend)
!     call recv_jas(itag,isend)

      if(iforce_analy.eq.1) then
        if(dmc_ivd.gt.0) then
          itag=itag+1
          call mpi_recv(pold(nwalk,1),PTH,mpi_double_precision,isend,itag,MPI_COMM_WORLD,istatus,ierr)
          itag=itag+1
          call mpi_recv(deriv_eold(1,1,nwalk),3*ncent,mpi_double_precision,isend,itag,MPI_COMM_WORLD,istatus,ierr)
          itag=itag+1
          call mpi_recv(esnake(1,1,nwalk,1),3*ncent*PTH,mpi_double_precision,isend,itag,MPI_COMM_WORLD,istatus,ierr)
          do ip=0,nwprod-1
            itag=itag+1
            call mpi_recv(ehist(1,1,nwalk,ip,1),3*ncent*PTH,mpi_double_precision,isend,itag,MPI_COMM_WORLD,istatus,ierr)
          enddo
        else
          itag=itag+1
          call mpi_recv(pold(nwalk,1),PTH,mpi_double_precision,isend,itag,MPI_COMM_WORLD,istatus,ierr)
        endif
      endif

      call prop_recv(isend,itag)
      call pcm_recv(isend,itag)
      call mmpol_recv(isend,itag)

      return
      end
end module
