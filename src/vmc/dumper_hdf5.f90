module dumper_hdf5_mod
        contains
        subroutine dumper_hdf5
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
        character(len=20), parameter   ::  restart_file_dmc = "restart_dmc.hdf5"
        character(len=20), parameter   ::  restart_file_vmc = "restart_dmc.hdf5"
        integer(hid_t)                 ::  file_id
        integer(hid_t)                 ::  group_id
        integer(hid_t)                 ::  plist_id
        integer(hid_t)                 ::  dataset_id
        integer(hid_t)                 ::  dataspace_id
        integer(hid_t)                 ::  memspace_id
        integer(hid_t)                 ::  filespace
        integer                        ::  error
        ! integer                        ::  nup, ndn, nelec
        character(len=20)              ::  author = "Ravindra Shinde", read_author

        ! integer :: i, ib, ic, id, ierr
        ! integer :: ifr, irequest, iw, j
        ! integer :: k, nscounts
        ! integer, dimension(4, 0:nproc) :: irn
        ! integer, dimension(MPI_STATUS_SIZE) :: istatus
        ! integer, dimension(4, 0:nproc) :: irn_tmp

        real, parameter :: zero = 0.0
        double precision, parameter :: one = 1.0d0
        real, dimension(2)  :: zero_array = (/0.0, 0.0/)
        integer, dimension(4)  :: int_array = (/-5, 4, 3, 10000/)
        double precision, dimension(2) :: one_array = (/1.0d0, 1.0d0/)
        real, dimension(3,3) :: real_2d_array = reshape((/ 1., 2., 3., 4., 5., 6., 7., 8., 9. /), shape(real_2d_array))
        logical, parameter :: debug = 0
        ! real(dp), parameter :: small = 1.e-6
        nup = 4
        ndn = 4
        nelec = 8

        call hdf5_file_create(restart_file_dmc, file_id)

        call hdf5_group_create(file_id, "Metadata", group_id)
        call hdf5_group_open(file_id, "Metadata", group_id)
        call hdf5_write(file_id, group_id, "Author", trim(author))
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "Electrons", group_id)
        call hdf5_group_open(file_id, "Electrons", group_id)
        call hdf5_write(file_id, group_id, "nup", nup)
        call hdf5_write(file_id, group_id, "ndn", ndn)
        call hdf5_write(file_id, group_id, "nelec", nelec)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "System", group_id)
        call hdf5_group_open(file_id, "System", group_id)
        call hdf5_write(file_id, group_id, "zero sp", zero)
        call hdf5_write(file_id, group_id, "one dp", one)
        call hdf5_write(file_id, group_id, "zero array", zero_array)
        call hdf5_write(file_id, group_id, "one array", one_array)
        call hdf5_write(file_id, group_id, "integer 1d array", int_array)
        call hdf5_write(file_id, group_id, "real 2d 1d array", real_2d_array)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "ECP", group_id)
        call hdf5_group_open(file_id, "ECP", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "Basis", group_id)
        call hdf5_group_open(file_id, "Basis", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "AO", group_id)
        call hdf5_group_open(file_id, "AO", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "MO", group_id)
        call hdf5_group_open(file_id, "MO", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "Determinants", group_id)
        call hdf5_group_open(file_id, "Determinants", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "CSFs", group_id)
        call hdf5_group_open(file_id, "CSFs", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "CSFmap", group_id)
        call hdf5_group_open(file_id, "CSFmap", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "States", group_id)
        call hdf5_group_open(file_id, "States", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "UnitCell", group_id)
        call hdf5_group_open(file_id, "UnitCell", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "Periodic", group_id)
        call hdf5_group_open(file_id, "Periodic", group_id)
        call hdf5_group_close(group_id)

        call hdf5_group_create(file_id, "QMC", group_id)
        call hdf5_group_open(file_id, "QMC", group_id)
        call hdf5_group_close(group_id)


        call hdf5_file_close(file_id)


        ! call hdf5_file_open(restart_file_dmc, file_id)
        ! call hdf5_group_open(file_id, "Metadata", group_id)
        ! call hdf5_read(file_id, group_id, "Author", read_author)
        ! call hdf5_group_close(group_id)


        ! call hdf5_file_close(file_id)

        end subroutine dumper_hdf5
end module dumper_hdf5_mod
