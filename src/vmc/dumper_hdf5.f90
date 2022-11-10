module dumper_hdf5_mod
        contains
        subroutine dumper_hdf5(restart_filename)
        !> @brief Dumps the data to a HDF5 file
        !> @details This subroutine dumps the data to a HDF5 file for restarting purposes
        !> @author Ravindra Shinde
        !> @date 2022-11-03
        !> @email r.l.shinde@utwente.nl

        use age,     only: iage,ioldest,ioldestmx
        use basis,   only: ns, np, nd, nf, ng, zex
        use branch,  only: eest,eigv,ff,fprod,nwalk,wdsumo,wgdsumo,wt
        use branch,  only: wtgen
        use coefs,   only: nbasis
        use csfs,    only: ncsf, nstates, ccsf
        use config,  only: xold_dmc
        use constants, only: hb
        use contrl_file, only: ounit
        use contrldmc, only: idmc,nfprod,rttau,tau
        use control, only: mode
        use control_dmc, only: dmc_nconf
        use custom_broadcast, only: bcast
        use denupdn, only: rprobdn,rprobup
        use derivest, only: derivcm2,derivcum,derivtotave_num_old
        use dmc_mod, only: MWALK
        ! use dumper_gpop_mod, only: dumper_gpop
        use est2cm,  only: ecm21_dmc,ecm2_dmc,efcm2,efcm21,egcm2,egcm21
        use est2cm,  only: ei1cm2,ei2cm2,ei3cm2,pecm2_dmc,r2cm2_dmc,ricm2
        use est2cm,  only: tjfcm_dmc,tpbcm2_dmc,wcm2,wcm21,wdcm2,wdcm21
        use est2cm,  only: wfcm2,wfcm21,wgcm2,wgcm21,wgdcm2
        use estcum,  only: ecum1_dmc,ecum_dmc,efcum,efcum1,egcum,egcum1
        use estcum,  only: ei1cum,ei2cum,ei3cum,iblk,ipass,pecum_dmc
        use estcum,  only: r2cum_dmc,ricum,taucum,tjfcum_dmc,tpbcum_dmc
        use estcum,  only: wcum1,wcum_dmc,wdcum,wdcum1,wfcum,wfcum1,wgcum
        use estcum,  only: wgcum1,wgdcum
        use jacobsave, only: ajacob
        use mmpol,   only: mmpol_dump
        use mpi
        use mpiblk,  only: iblk_proc
        use mpiconf, only: idtask,nproc,wid
        use multiple_geo, only: fgcm2,fgcum,nforce,pecent
        use pcm_mod, only: pcm_dump
        use precision_kinds, only: dp
        use properties_mod, only: prop_dump
        use pseudo,  only: nloc
        use qua,     only: nquad,wq,xq,yq,zq
        use rannyu_mod, only: savern
        use slater,  only: cdet,coef,ndet,norb
        use stats,   only: acc,dfus2ac,dfus2un,dr2ac,dr2un,nacc,nbrnch
        use stats,   only: nodecr,trymove
        use step,    only: rprob
        use strech_mod, only: strech
        use system,  only: cent,iwctype,ncent,nctype,ndn,nelec,newghostype
        use system,  only: nghostcent,nup,znuc
        use velratio, only: fratio
        use vmc_mod, only: nrad
        use hdf5
        use hdf5_utils, only: hdf5_file_create, hdf5_file_close, hdf5_file_open
        use hdf5_utils, only: hdf5_group_create, hdf5_group_close, hdf5_group_open
        use hdf5_utils, only: hdf5_write!, hdf5_read


        implicit none

        ! HDF5 related variables
        character(len=*), intent(in)   ::  restart_filename
        integer(hid_t)                 ::  file_id
        integer(hid_t)                 ::  group_id
        integer(hid_t)                 ::  plist_id
        integer(hid_t)                 ::  dataset_id
        integer(hid_t)                 ::  dataspace_id
        integer(hid_t)                 ::  memspace_id
        integer(hid_t)                 ::  filespace
        integer                        ::  error
        character(len=20)              ::  author = "Ravindra Shinde", read_author

        integer :: i, ib, ic, id, ierr
        integer :: ifr, irequest, iw, j
        integer :: k, nscounts
        integer, dimension(4, 0:nproc) :: irn
        integer, dimension(MPI_STATUS_SIZE) :: istatus
        integer, dimension(4, 0:nproc) :: irn_tmp

        real(dp), parameter             :: zero = 0.0d0
        real(dp), parameter             :: one  = 1.0d0
        real(dp), parameter             :: small = 1.0d-6

        if(nforce.gt.1) call strech(xold_dmc,xold_dmc,ajacob,1,0)

        call savern(irn(1,idtask))

        nscounts=4
        call mpi_gather(irn(1,idtask),nscounts,mpi_integer,irn_tmp,nscounts,mpi_integer,0,MPI_COMM_WORLD,ierr)

        ! Bradcast the data which is spread across the processors
        call bcast(nwalk)
        call bcast(xold_dmc)
        call bcast(wt)
        call bcast(ff(0))
        call bcast(fprod)
        call bcast(fratio)
        call bcast(eigv)
        call bcast(eest)
        call bcast(wdsumo)
        call bcast(iage)
        call bcast(ioldest)
        call bcast(ioldestmx)
        call bcast(xq)
        call bcast(yq)
        call bcast(zq)

        ! Only the master process will write the data to the HDF5 file
        if (wid) then
        ! Open the HDF5 file
        call hdf5_file_create(restart_filename, file_id)

        call hdf5_group_create(file_id, "Metadata", group_id)
        call hdf5_group_open(file_id, "Metadata", group_id)

        call get_environment_variable ("USER", author)
        call hdf5_write(file_id, group_id, "Author", trim(author))
        call hdf5_write(file_id, group_id, " Code Compilation Date ", __DATE__)
        call hdf5_write(file_id, group_id, " Code Compilation Time ", __TIME__)

#if defined(GIT_HEAD_BRANCH)
        call hdf5_write(file_id, group_id, " Git Branch ", GIT_HEAD_BRANCH)
#endif

#if defined(GIT_REVISION_HASH)
        call hdf5_write(file_id, group_id, " Git Commit Hash ", GIT_REVISION_HASH)
#endif

#if defined(CMAKE_Fortran_COMPILER)
        call hdf5_write(file_id, group_id, " Compiler ", CMAKE_Fortran_COMPILER)
#endif

#if defined(CMAKE_Fortran_COMPILER_VERSION)
        call hdf5_write(file_id, group_id, " Compiler Version ", CMAKE_Fortran_COMPILER_VERSION)
#endif

#if defined(TARGET_ARCHITECTURE)
        call hdf5_write(file_id, group_id, " Vectorization Instructions ", TARGET_ARCHITECTURE)
#endif

#if defined(HDF5_VERSION)
        call hdf5_write(file_id, group_id, " HDF5 Version ", HDF5_VERSION)
#endif
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "Electrons", group_id)
        call hdf5_group_open(file_id, "Electrons", group_id)
        call hdf5_write(file_id, group_id, "Number of Up-Spin Electrons", nup)
        call hdf5_write(file_id, group_id, "Number of Down-Spin Electrons", ndn)
        call hdf5_write(file_id, group_id, "Total Number of Electrons", nelec)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "System", group_id)
        call hdf5_group_open(file_id, "System", group_id)
        call hdf5_write(file_id, group_id, "Number of Center Types", nctype)
        call hdf5_write(file_id, group_id, "Number of Centers", ncent)
        call hdf5_write(file_id, group_id, "Center Coordinates", cent)
        call hdf5_write(file_id, group_id, "Number of Ghost Center Types", newghostype)
        call hdf5_write(file_id, group_id, "Number of Ghost Centers", nghostcent)
        call hdf5_write(file_id, group_id, "Index of Which Center Type", iwctype)
        call hdf5_write(file_id, group_id, "PE Centers", pecent)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "ECP", group_id)
        call hdf5_group_open(file_id, "ECP", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "Basis", group_id)
        call hdf5_group_open(file_id, "Basis", group_id)
        call hdf5_write(file_id, group_id, "Zex", zex)
        call hdf5_write(file_id, group_id, "nquad", nquad)
        call hdf5_write(file_id, group_id, "xq", xq)
        call hdf5_write(file_id, group_id, "yq", yq)
        call hdf5_write(file_id, group_id, "zq", zq)
        call hdf5_write(file_id, group_id, "wq", wq)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "AO", group_id)
        call hdf5_group_open(file_id, "AO", group_id)
        call hdf5_write(file_id, group_id, "Number of Basis", nbasis)
        call hdf5_write(file_id, group_id, "Number of S Type AOs", ns)
        call hdf5_write(file_id, group_id, "Number of P Type AOs", np)
        call hdf5_write(file_id, group_id, "Number of D Type AOs", nd)
        call hdf5_write(file_id, group_id, "Number of F Type AOs", nf)
        call hdf5_write(file_id, group_id, "Number of G Type AOs", ng)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "MO", group_id)
        call hdf5_group_open(file_id, "MO", group_id)
        call hdf5_write(file_id, group_id, "MO Coefficients", coef)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "Determinants", group_id)
        call hdf5_group_open(file_id, "Determinants", group_id)
        call hdf5_write(file_id, group_id, "Number of Determinants", ndet)
        call hdf5_write(file_id, group_id, "Determinant Coefficients", cdet)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "CSFs", group_id)
        call hdf5_group_open(file_id, "CSFs", group_id)
        call hdf5_write(file_id, group_id, "Number of CSFs", ncsf)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "States", group_id)
        call hdf5_group_open(file_id, "States", group_id)
        call hdf5_write(file_id, group_id, "Number of States", nstates)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "UnitCell", group_id)
        call hdf5_group_open(file_id, "UnitCell", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "Periodic", group_id)
        call hdf5_group_open(file_id, "Periodic", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "QMC", group_id)
        call hdf5_group_open(file_id, "QMC", group_id)
        call hdf5_write(file_id, group_id, "Number of Processors", nproc)
        call hdf5_write(file_id, group_id, "Number of Walkers", nwalk)
        call hdf5_write(file_id, group_id, "xold_dmc", xold_dmc)
        call hdf5_write(file_id, group_id, "nfprod", nfprod)
        call hdf5_write(file_id, group_id, "ff", ff)
        call hdf5_write(file_id, group_id, "wt", wt)
        call hdf5_write(file_id, group_id, "fprod", fprod)
        call hdf5_write(file_id, group_id, "eigv", eigv)
        call hdf5_write(file_id, group_id, "eest", eest)
        call hdf5_write(file_id, group_id, "wdsumo", wdsumo)
        call hdf5_write(file_id, group_id, "iage", iage)
        call hdf5_write(file_id, group_id, "ioldest", ioldest)
        call hdf5_write(file_id, group_id, "ioldestmx", ioldestmx)
        call hdf5_write(file_id, group_id, "nforce", nforce)
        call hdf5_write(file_id, group_id, "fratio", fratio)
        call hdf5_group_close(group_id)


        call hdf5_file_close(file_id)
        ! Close the HDF5 file
        endif

        end subroutine dumper_hdf5
end module dumper_hdf5_mod
