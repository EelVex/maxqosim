! Main template
! for solving master equations
!
! Part of maxqosim Maxima package for open quantum systems
!
! Copyright 2020 by Kostas Blekos, eelvex@gmail.com
!
! Released under the terms of the GNU General Public License
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!/
program lindblad_template

    use integrators
    implicit none

    ! Defined in Maxima:
#include "parameters.F90"
    integer,parameter::N = mx_Nout
    real(dp),parameter::t_init = mx_t0
    real(dp),parameter::t_fin = mx_tfin
    real(dp),parameter::dt = mx_dt
    complex(dp),dimension(:,:),allocatable::rho0

    complex(dp),dimension(:,:,:),allocatable::output
    real(dp),dimension(:),allocatable::times

    integer::i,fu=2

    open(unit=fu, file=mx_file)

    rho0 = mx_rho0

    ! Warning: Older fortran standards require allocation of 'output' before the assignment
    allocate(output(size(rho0,1),size(rho0,2), N))
    output = RKH4_integrate(t_init, t_fin, rho0, rhodot, dt, N)
    times = space(t_init, t_fin, N)

    do i=1,N
        write(fu,*)times(i),abs(output(:,:,i))
    enddo

contains

    function rhodot(t, rho)
        real(dp),intent(in)::t
        complex(dp),dimension(:,:),intent(in)::rho
        complex(dp),dimension(size(rho,1),size(rho,2))::rhodot

#include "equations.F90"

    end function

    pure function space(t_ini, t_fin, N)
        ! Split the interval [t_ini, t_fin] to N linear-spaced values
        ! t_ini and t_fin are included
        real(dp),intent(in)::t_ini, t_fin
        integer,intent(in)::N
        real(dp),dimension(N)::space
        integer::i

        associate(dd => (t_fin - t_ini) / (N-1))
            space(:) = [ ((i-1)*dd + t_ini, i=1,N) ]
        end associate

    end function


end program lindblad_template
! EOF
