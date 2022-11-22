module dumper_hdf5_mod
        contains
        subroutine vmc_dumper_hdf5(restart_filename)
        !> @brief Dumps the data to a HDF5 file
        !> @details This subroutine dumps the data to a HDF5 file for restarting purposes
        !> @author Ravindra Shinde
        !> @date 2022-11-03
        !> @email r.l.shinde@utwente.nl

        use mpi
        use hdf5
        use custom_broadcast, only: bcast
        use hdf5_utils, only: hdf5_file_create, hdf5_file_close, hdf5_file_open
        use hdf5_utils, only: hdf5_group_create, hdf5_group_close, hdf5_group_open
        use hdf5_utils, only: hdf5_write, hdf5_read

        ! init routines
        use mmpol,   only: mmpol_init
        use optci_mod, only: optci_init
        use optjas_mod, only: optjas_init
        use optorb_f_mod, only: optorb_init
        use pcm_mod, only: pcm_init
        use properties_mod, only: prop_init

        ! union of all the required arrays
        use basis,   only: ns, np, nd, nf, ng, zex
        use coefs,   only: nbasis
        use config,  only: eold,nearesto,psi2o,psido,psijo,rmino,rvmino
        use config,  only: tjfo,vold,xnew,xold
        use constants, only: hb
        use contrl_file, only: errunit,ounit
        use control, only: mode
        use control_vmc, only: vmc_idump, vmc_irstar, vmc_isite, vmc_nconf, vmc_nblk, vmc_nblk_max
        use control_vmc, only: vmc_nblkeq, vmc_nconf_new, vmc_nstep, vmc_icharged_atom, vmc_nblk_ci
        use csfs,    only: ncsf, nstates, ccsf
        use denupdn, only: rprobdn,rprobup
        ! use determinant_psig_mod, only: determinant_psig
        use determinante_mod, only: compute_determinante_grad
        use error,   only: fatal_error
        use est2cm,  only: ecm2,ecm21,pecm2,r2cm2,tjfcm2,tpbcm2
        use estcum,  only: ecum,ecum1,iblk,pecum,r2cum,tjfcum,tpbcum
        use estsig,  only: ecm21s,ecum1s
        use estsum,  only: acc,esum,pesum,r2sum,tjfsum,tpbsum
        use force_analytic, only: force_analy_dump,force_analy_rstrt
        use forcewt, only: wcum,wsum
        use hpsi_mod, only: hpsi
        use inputflags, only: eps_node_cutoff,node_cutoff
        use metropolis, only: delta,deltar,deltat
        use mstates_ctrl, only: iguiding
        use mstates_mod, only: MSTATES
        use multiple_geo, only: fcm2,fcum,fgcm2,fgcum,iwftype,nforce,nwftype,pecent
        ! use multiple_states, only: efficiency_dump,efficiency_rstrt
        use mpiconf, only: idtask,nproc,wid
        use nodes_distance_mod, only: nodes_distance,rnorm_nodes_num
        use optci_mod, only: optci_dump,optci_rstrt,optci_save
        use optjas_mod, only: optjas_dump,optjas_rstrt,optjas_save
        use optorb_cblock, only: ns_current
        use optorb_f_mod, only: optorb_dump,optorb_rstrt,optorb_save
        use optwf_control, only: ioptorb
        use optx_jas_ci, only: optx_jas_ci_dump,optx_jas_ci_rstrt
        use optx_jas_orb, only: optx_jas_orb_dump,optx_jas_orb_rstrt
        use optx_orb_ci, only: optx_orb_ci_dump,optx_orb_ci_rstrt
        use pcm_mod, only: pcm_dump,pcm_rstrt
        use precision_kinds, only: dp
        ! use prop_vmc, only: prop_save
        use properties_mod, only: prop_dump,prop_rstrt
        use pseudo,  only: nloc
        use qua,     only: nquad,wq,xq,yq,zq
        use rannyu_mod, only: rannyu,savern,setrn
        use slater,  only: cdet,coef,ndet,norb
        use stats,   only: rejmax
        use step,    only: ekin,ekin2,rprob,suc,trunfb,try
        use system,  only: cent,iwctype,ncent,ncent_tot,nctype,nctype_tot,ndn,nelec,newghostype
        use system,  only: nghostcent,nup,znuc
        use strech_mod, only: setup_force,strech
        use vmc_mod, only: norb_tot,nrad


        implicit none

        ! HDF5 related variables
        character(len=*), intent(in)   ::  restart_filename
        integer(hid_t)                 ::  file_id
        integer(hid_t)                 ::  group_id, group_id1, group_id2
        integer(hid_t)                 ::  plist_id
        integer(hid_t)                 ::  dataset_id
        integer(hid_t)                 ::  dataspace_id
        integer(hid_t)                 ::  memspace_id
        integer(hid_t)                 ::  filespace
        character(len=20)              ::  author = "Ravindra Shinde", read_author

        integer :: i, ib, ic, id, idfrom, idget, ierr
        integer :: ifr, istate, j, k
        integer :: nelecx, nforcex, nlocx, nproco
        integer :: nq_id, nqd_id, nqx, nscounts
        integer, dimension(4,0:nproc) :: irn
        integer, dimension(MPI_STATUS_SIZE) :: istatus
        integer, dimension(4,0:nproc) :: irn_tmp
        integer, dimension(0:nproc) :: ircounts
        integer, dimension(0:nproc) :: idispls
        integer :: irequest, iw
        real(dp) :: rnd, wq_id, x_id, xq_id, yq_id, zq_id


        do i=0,nproc-1
            ircounts(i)=4
            idispls(i)=i*4
        enddo

        idispls(nproc)=4*nproc
        nscounts=ircounts(idtask)

        call savern(irn(1,idtask))

        call mpi_gatherv(irn(1,idtask),nscounts,mpi_integer,irn_tmp,ircounts,idispls,mpi_integer,0,MPI_COMM_WORLD,ierr)

        ! Bradcast the data which is spread across the processors
        call bcast(xold)
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
        call hdf5_write(file_id, group_id, "Nforce", nforce)
        call hdf5_write(file_id, group_id, "Nloc", nloc)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "ECP", group_id)
        call hdf5_group_open(file_id, "ECP", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "Basis", group_id)
        call hdf5_group_open(file_id, "Basis", group_id)
        call hdf5_write(file_id, group_id, "Zex", zex)
        if (nloc .gt. 0) then
            call hdf5_write(file_id, group_id, "nquad", nquad)
            call hdf5_write(file_id, group_id, "xq", xq(1:nquad))
            call hdf5_write(file_id, group_id, "yq", yq(1:nquad))
            call hdf5_write(file_id, group_id, "zq", zq(1:nquad))
            call hdf5_write(file_id, group_id, "wq", wq(1:nquad))
        endif
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
        call hdf5_write(file_id, group_id, "Number of Orbitals", norb)
        call hdf5_write(file_id, group_id, "Number of Orbitals Total", norb_tot)
        call hdf5_write(file_id, group_id, "MO Coefficients", coef(1:nbasis,1:norb,1))
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "Determinants", group_id)
        call hdf5_group_open(file_id, "Determinants", group_id)
        call hdf5_write(file_id, group_id, "Number of Determinants", ndet)
        call hdf5_write(file_id, group_id, "Determinant Coefficients", cdet(1:ndet,1,1))
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
        call hdf5_write(file_id, group_id, "Mode", mode)
        call hdf5_write(file_id, group_id, "Number of Processors", nproc)
        call hdf5_write(file_id, group_id, "Random Numbers Each Processor", irn_tmp(1:4,0:nproc-1))
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "VMC", group_id)
        call hdf5_group_open(file_id, "VMC", group_id)
        call hdf5_write(file_id, group_id, "Number of VMC Blocks", vmc_nblk)
        call hdf5_write(file_id, group_id, "Number of Equilibration VMC Blocks", vmc_nblkeq)
        call hdf5_write(file_id, group_id, "Maximum Number of VMC Blocks", vmc_nblk_max)
        call hdf5_write(file_id, group_id, "Number of VMC Steps per Block", vmc_nstep)
        call hdf5_write(file_id, group_id, "Number of VMC Configurations", vmc_nconf)
        call hdf5_write(file_id, group_id, "Number of New VMC Configurations", vmc_nconf_new)
        call hdf5_write(file_id, group_id, "xold", xold)
        call hdf5_write(file_id, group_id, "Delta", delta)
        call hdf5_write(file_id, group_id, "Deltar", deltar)
        call hdf5_write(file_id, group_id, "Deltat", deltat)


        call hdf5_write(file_id, group_id, "ecum1", ecum1(1:nstates))
        call hdf5_write(file_id, group_id, "ecum", ecum(1:nstates,1:nforce))
        call hdf5_write(file_id, group_id, "pecum", pecum(1:nstates))
        call hdf5_write(file_id, group_id, "tpbcum", tpbcum(1:nstates))
        call hdf5_write(file_id, group_id, "tjfcum", tjfcum(1:nstates))
        call hdf5_write(file_id, group_id, "r2cum", r2cum)
        call hdf5_write(file_id, group_id, "acc", acc)
        call hdf5_write(file_id, group_id, "ecm21", ecm21(1:nstates))
        call hdf5_write(file_id, group_id, "ecm2", ecm2(1:nstates,1:nforce))
        call hdf5_write(file_id, group_id, "pecm2", pecm2(1:nstates))
        call hdf5_write(file_id, group_id, "tpbcm2", tpbcm2(1:nstates))
        call hdf5_write(file_id, group_id, "tjfcm2", tjfcm2(1:nstates))
        call hdf5_write(file_id, group_id, "r2cm2", r2cm2)

        if (nforce .gt. 1) then
            call hdf5_write(file_id, group_id, "wcum", wcum(1:nstates,1:nforce))
            call hdf5_write(file_id, group_id, "fcum", fcum(1:nstates,1:nforce))
            call hdf5_write(file_id, group_id, "fcm2", fcm2(1:nstates,1:nforce))
        else
            call hdf5_write(file_id, group_id, "wcum", wcum(1:nstates,1))
        end if

        call hdf5_write(file_id, group_id, "ecum1s", ecum1s(1:nstates))
        call hdf5_write(file_id, group_id, "ecm21s", ecm21s(1:nstates))

        call hdf5_write(file_id, group_id, "try", try(1:nrad))
        call hdf5_write(file_id, group_id, "suc", suc(1:nrad))
        call hdf5_write(file_id, group_id, "trunfb", trunfb(1:nrad))
        call hdf5_write(file_id, group_id, "rprob", rprob(1:nrad))
        call hdf5_write(file_id, group_id, "rprobup", rprobup(1:nrad))
        call hdf5_write(file_id, group_id, "rprobdn", rprobdn(1:nrad))
        call hdf5_write(file_id, group_id, "ekin", ekin(1:nrad))
        call hdf5_write(file_id, group_id, "ekin2", ekin2(1:nrad))
        call hdf5_write(file_id, group_id, "rejmax", rejmax)


        ! call optorb_dump(10)
        ! call optci_dump(10)
        ! call prop_dump(10)
        ! call efficiency_dump(10)
        ! call pcm_dump(10)
        ! call force_analy_dump(10)
        ! call optjas_dump(10)
        ! call optx_jas_orb_dump(10)
        ! call optx_jas_ci_dump(10)
        ! call optx_orb_ci_dump(10)

        call hdf5_group_close(group_id)

        call hdf5_file_close(file_id)
        ! Close the HDF5 file

        ! Reading section trial

        ! Open hdf5 file for reading data
        call hdf5_file_open(restart_filename, file_id)
        ! Open the group
        call hdf5_group_open(file_id, "System", group_id)
        ! Read the data
        call hdf5_read(file_id, group_id, "Number of Centers", temporary_integer)
        print*, " Number of Centers read from hdf5 ", temporary_integer
        call hdf5_read(file_id, group_id, "Index of Which Center Type", temporary_integer_array)
        print*, " Number of Centers types read from hdf5 ", temporary_integer_array
        ! Close the group
        call hdf5_group_close(group_id)
        ! Close the file
        call hdf5_file_close(file_id)



        endif   ! master thread

        end subroutine vmc_dumper_hdf5
end module dumper_hdf5_mod
