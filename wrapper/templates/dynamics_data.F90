module dynamics_data_mod

use atmosphere_mod, only: Atm, mygrid
use fv_grid_utils_mod, only: cubed_a2d, cubed_to_latlon
use mpp_domains_mod, only: DGRID_NE, mpp_update_domains
use fv_mp_mod,       only: start_group_halo_update, complete_group_halo_update, group_halo_update_type
use tracer_manager_mod, only: get_tracer_names, get_number_tracers, get_tracer_index
use field_manager_mod,  only: MODEL_ATMOS
use iso_c_binding

implicit none

contains

    pure function i_start() result(i)
        integer :: i
        i = Atm(mygrid)%bd%is
    end function i_start

    pure function i_end() result(i)
        integer :: i
        i = Atm(mygrid)%bd%ie
    end function i_end

    pure function j_start() result(j)
        integer :: j
        j = Atm(mygrid)%bd%js
    end function j_start

    pure function j_end() result(j)
        integer :: j
        j = Atm(mygrid)%bd%je
    end function j_end

    pure function nz() result(n)
        integer :: n
        n = Atm(mygrid)%flagstruct%npz
    end function nz

    subroutine get_centered_grid_dimensions(nx_out, ny_out, nz_out) bind(c)
        integer(c_int), intent(out) :: nx_out, ny_out, nz_out
        nx_out = i_end() - i_start() + 1 ! Fortran end index is inclusive
        ny_out = j_end() - j_start() + 1
        nz_out = nz()
    end subroutine get_centered_grid_dimensions

    subroutine get_n_ghost_cells_subroutine(n_ghost_out) bind(c)
        integer(c_int), intent(out) :: n_ghost_out
        n_ghost_out = Atm(mygrid)%bd%ng
    end subroutine get_n_ghost_cells_subroutine

{% for item in dynamics_properties %}
    subroutine set_{{ item.fortran_name }}({{ item.fortran_name }}_in) bind(c)
        real(c_double), intent(in), dimension({{ item.dim_ranges }}) :: {{ item.fortran_name }}_in
        Atm(mygrid)%{{ item.fortran_name }}({{ item.dim_ranges }}) = {{ item.fortran_name }}_in({{ item.dim_ranges }})
    end subroutine set_{{ item.fortran_name }}

    subroutine get_{{ item.fortran_name }}({{ item.fortran_name }}_out) bind(c)
        real(c_double), intent(out), dimension({{ item.dim_ranges }}) :: {{ item.fortran_name }}_out
        {{ item.fortran_name }}_out({{ item.dim_ranges }}) = Atm(mygrid)%{{ item.fortran_name }}({{ item.dim_ranges }})
    end subroutine get_{{ item.fortran_name }}
{% endfor %}

    subroutine get_tracer_count(n_prognostic_tracers, n_total_tracers) bind(c)
        integer(c_int), intent(out) :: n_prognostic_tracers, n_total_tracers
        call get_number_tracers(MODEL_ATMOS, num_tracers=n_total_tracers, num_prog=n_prognostic_tracers)
    end subroutine get_tracer_count

    subroutine get_tracer_breakdown(nwat, dnats, pnats) bind(c)
        integer(c_int), intent(out) :: nwat, dnats, pnats
        nwat = Atm(mygrid)%flagstruct%nwat
        dnats = Atm(mygrid)%flagstruct%dnats
        pnats = Atm(mygrid)%flagstruct%pnats
    end subroutine get_tracer_breakdown

    subroutine get_tracer(tracer_index, array_out) bind(c)
        ! get tracer at the given one-based index
        real(c_double), intent(out) :: array_out(i_start():i_end(), j_start():j_end(), nz())
        integer(c_int), intent(in) :: tracer_index
        integer(c_int) :: n_prognostic_tracers, n_total_tracers
        call get_tracer_count(n_prognostic_tracers, n_total_tracers)
        if (tracer_index <= n_prognostic_tracers) then
            array_out(:, :, :) = Atm(mygrid)%q(i_start():i_end(), j_start():j_end(), 1:nz(), tracer_index)
        else
            array_out(:, :, :) = Atm(mygrid)%qdiag(i_start():i_end(), j_start():j_end(), 1:nz(), tracer_index)
        end if
    end subroutine get_tracer

    subroutine set_tracer(tracer_index, array_in) bind(c)
        ! set tracer at the given one-based index
        real(c_double), intent(in) :: array_in(i_start():i_end(), j_start():j_end(), nz())
        integer(c_int), intent(in) :: tracer_index
        integer(c_int) :: n_prognostic_tracers, n_total_tracers
        call get_tracer_count(n_prognostic_tracers, n_total_tracers)
        if (tracer_index <= n_prognostic_tracers) then
            Atm(mygrid)%q(i_start():i_end(), j_start():j_end(), 1:nz(), tracer_index) = array_in
        else
            Atm(mygrid)%qdiag(i_start():i_end(), j_start():j_end(), 1:nz(), tracer_index) = array_in
        end if
    end subroutine set_tracer

    subroutine get_tracer_name(tracer_index, tracer_name_out, tracer_long_name_out, tracer_units_out) bind(c)
        integer(c_int), intent(in) :: tracer_index
        character(kind=c_char), dimension(64) :: tracer_name_out, tracer_long_name_out, tracer_units_out
        character(kind=c_char, len=64) :: tracer_name, tracer_long_name, tracer_units
        integer i
        call get_tracer_names(model=MODEL_ATMOS, n=tracer_index, name=tracer_name, longname=tracer_long_name, units=tracer_units)
        tracer_name = trim(tracer_name) // c_null_char
        tracer_long_name = trim(tracer_long_name) // c_null_char
        tracer_units = trim(tracer_units) // c_null_char
        do i= 1, 64
            tracer_name_out(i) = tracer_name(i:i)
            tracer_long_name_out(i) = tracer_long_name(i:i)
            tracer_units_out(i) = tracer_units(i:i)
        enddo
    end subroutine get_tracer_name

    subroutine transform_agrid_winds_to_dgrid_winds_subroutine(ua, va, u, v) bind(c)
        ! Wraps the internal cubed_a2d subroutine in a more convenient way, handling
        ! halo updates and the additional grid arguments within fortran for convenience.
        real, intent(in), dimension(i_start():i_end(),j_start():j_end(),1:nz()) :: ua, va
        real, intent(out), dimension(i_start():i_end(),j_start():j_end()+1,1:nz()) :: u
        real, intent(out), dimension(i_start():i_end()+1,j_start():j_end(),1:nz()) :: v

        real, allocatable, dimension(:,:,:) :: ua_halo, va_halo, u_halo, v_halo

        integer :: is, ie, js, je, isd, ied, jsd, jed, npx, npy, npz

        is = i_start()
        ie = i_end()
        js = j_start()
        je = j_end()
        isd = Atm(mygrid)%bd%isd
        ied = Atm(mygrid)%bd%ied
        jsd = Atm(mygrid)%bd%jsd
        jed = Atm(mygrid)%bd%jed
        npx = Atm(mygrid)%npx
        npy = Atm(mygrid)%npy
        npz = nz()

        allocate(ua_halo(isd:ied,jsd:jed,1:npz))
        allocate(va_halo(isd:ied,jsd:jed,1:npz))
        allocate(u_halo(isd:ied,jsd:jed+1,1:npz))
        allocate(v_halo(isd:ied+1,jsd:jed,1:npz))

        ! Note we do not need to do a halo update here, since cubed_a2d takes
        ! of this internally.
        ua_halo(is:ie,js:je,1:npz) = ua
        va_halo(is:ie,js:je,1:npz) = va

        call cubed_a2d(&
            npx, &
            npy, &
            npz, &
            ua_halo, &
            va_halo, &
            u_halo, &
            v_halo, &
            Atm(mygrid)%gridstruct, &
            Atm(mygrid)%domain, &
            Atm(mygrid)%bd &
        )

        u = u_halo(is:ie,js:je+1,1:npz)
        v = v_halo(is:ie+1,js:je,1:npz)
    end subroutine transform_agrid_winds_to_dgrid_winds_subroutine

    subroutine transform_dgrid_winds_to_agrid_winds_subroutine(u, v, ua, va) bind(c)
        ! Wraps the internal cubed_to_latlon subroutine in a more convenient way, 
        ! handling halo updates and the additional grid arguments within fortran
        ! for convenience.
        real, intent(in), dimension(i_start():i_end(),j_start():j_end()+1,1:nz()) :: u
        real, intent(in), dimension(i_start():i_end()+1,j_start():j_end(),1:nz()) :: v
        real, intent(out), dimension(i_start():i_end(),j_start():j_end(),1:nz()) :: ua, va

        real, allocatable, dimension(:,:,:) :: ua_halo, va_halo, u_halo, v_halo

        type(group_halo_update_type), save :: i_pack
        integer :: is, ie, js, je, isd, ied, jsd, jed, npx, npy, npz, mode, c2l_ord

        is = i_start()
        ie = i_end()
        js = j_start()
        je = j_end()
        isd = Atm(mygrid)%bd%isd
        ied = Atm(mygrid)%bd%ied
        jsd = Atm(mygrid)%bd%jsd
        jed = Atm(mygrid)%bd%jed
        npx = Atm(mygrid)%npx
        npy = Atm(mygrid)%npy
        npz = nz()
        mode = 1
        c2l_ord = 4  ! Recommended default

        allocate(ua_halo(isd:ied,jsd:jed,1:npz))
        allocate(va_halo(isd:ied,jsd:jed,1:npz))
        allocate(u_halo(isd:ied,jsd:jed+1,1:npz))
        allocate(v_halo(isd:ied+1,jsd:jed,1:npz))

        u_halo(is:ie,js:je+1,1:npz) = u
        v_halo(is:ie+1,js:je,1:npz) = v
        call start_group_halo_update(i_pack, u_halo, v_halo, Atm(mygrid)%domain, gridtype=DGRID_NE)
        call complete_group_halo_update(i_pack, Atm(mygrid)%domain)

        call cubed_to_latlon(&
            u_halo, &
            v_halo, &
            ua_halo, &
            va_halo, &
            Atm(mygrid)%gridstruct, &
            npx, &
            npy, &
            npz, &
            mode, &
            Atm(mygrid)%gridstruct%grid_type, &
            Atm(mygrid)%domain, &
            Atm(mygrid)%gridstruct%nested, &
            c2l_ord, &
            Atm(mygrid)%bd &
        )

        ua = ua_halo(is:ie,js:je,1:npz)
        va = va_halo(is:ie,js:je,1:npz)
    end subroutine transform_dgrid_winds_to_agrid_winds_subroutine

end module dynamics_data_mod
