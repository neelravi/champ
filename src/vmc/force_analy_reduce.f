      module force_analy_reduce_mod
      contains
      subroutine force_analy_reduce
      use atom, only: ncent, ncent_tot
      use vd_mod, only: dmc_ivd, da_branch_cum
      use mpiconf, only: wid
      use da_energy_sumcum, only: da_energy_cm2, da_energy_cum, da_psi_cum
      use force_analy, only: iforce_analy
      use force_pth, only: PTH
      use mpi

      use precision_kinds, only: dp
      implicit none

      integer :: ic, ierr, ii, k, iph

      real(dp), dimension(3*ncent_tot*PTH) :: collect



      if(iforce_analy.eq.0) return

      call mpi_reduce(da_energy_cum,collect,3*ncent*PTH
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            ii=ii+1
            da_energy_cum(k,ic,iph)=collect(ii)
          enddo
        enddo
        enddo
       else
        do iph=1,PTH               
        do ic=1,ncent
          do k=1,3
            da_energy_cum(k,ic,iph)=0.d0
          enddo
        enddo
        enddo
      endif

      call mpi_reduce(da_psi_cum,collect,3*ncent*PTH
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            ii=ii+1
            da_psi_cum(k,ic,iph)=collect(ii)
          enddo
        enddo
        enddo
       else
        do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            da_psi_cum(k,ic,iph)=0.d0
          enddo
        enddo
        enddo
      endif

      if (dmc_ivd.gt.0) then
      call mpi_reduce(da_branch_cum,collect,3*ncent*PTH
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            ii=ii+1
            da_branch_cum(k,ic,iph)=collect(ii)
          enddo
        enddo
        enddo
       else
        do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            da_branch_cum(k,ic,iph)=0.d0
          enddo
        enddo
        enddo
      endif
      endif

      call mpi_reduce(da_energy_cm2,collect,3*ncent*PTH
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            ii=ii+1
            da_energy_cm2(k,ic,iph)=collect(ii)
          enddo
        enddo
        enddo
       else
        do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            da_energy_cm2(k,ic,iph)=0.d0
          enddo
        enddo
        enddo
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
      end module
