! IMPLEMENTS Runge Kutta 4-th order method of integration
! for solving master equations, i.e. systems NxN
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
module integrators

    implicit none

    integer,parameter::dp=kind(0.0d0)
    complex(dp),parameter:: I_c=(0._dp,1._dp) ! i

    private

    public RKH4_step, RKH4_integrate
    public dp, I_c

    interface
        function rhodot_if(t, rho)
            implicit none
            integer,parameter::dp=kind(0.0d0)
            real(dp),intent(in)::t
            complex(dp),dimension(:,:),intent(in)::rho
            complex(dp),dimension(size(rho,1),size(rho,2))::rhodot_if
        end function
    end interface

contains

    function RKH4_step(t,rho,rhodot,dt) result(rho_new)

        complex(dp),dimension(:,:),intent(in)::rho
        complex(dp),dimension(size(rho,1),size(rho,2))::rho_new
        complex(dp),dimension(size(rho,1),size(rho,2),4)::k
        complex(dp),dimension(size(rho,1),size(rho,2))::Delta
        real(dp),intent(in)::t, dt
        procedure(rhodot_if)::rhodot
        real(dp),dimension(4),parameter:: rk4_t = (/ 0._dp, 0.5_dp, 0.5_dp, 1._dp /)
        ! real(dp),dimension(4),parameter:: rk4_d = (/ 1._dp, 2.0_dp, 2.0_dp, 1._dp /)
        integer::i

        do i=1,4
            k(:,:,i) = rhodot(t + rk4_t(i)*dt, rho + rk4_t(i)*k(:,:,i-1)*dt)
        enddo

        Delta = (1._dp/6._dp)*dt*(k(:,:,1) + 2*k(:,:,2) + 2*k(:,:,3) + k(:,:,4))

        rho_new = rho + Delta

    end function

    function RKH4_integrate(t_ini, t_fin, rho0, rhodot, dt, N) result(output)
        real(dp),intent(in)::t_ini, t_fin, dt
        integer,intent(in)::N
        complex(dp),dimension(:,:),intent(in)::rho0
        procedure(rhodot_if)::rhodot
        complex(dp),dimension(size(rho0,1),size(rho0,2))::rho
        complex(dp),dimension(size(rho0,1),size(rho0,2),4)::k
        complex(dp),dimension(size(rho0,1),size(rho0,2))::Delta
        real(dp),dimension(4),parameter:: rk4_t = (/ 0._dp, 0.5_dp, 0.5_dp, 1._dp /)
        ! real(dp),dimension(4),parameter:: rk4_d = (/ 1._dp, 2.0_dp, 2.0_dp, 1._dp /)
        complex(dp),dimension(size(rho0,1),size(rho0,2),N)::output
        ! complex(dp),dimension(:,:,:),allocatable::output
        real(dp)::t, dd
        integer::i, iter

        dd = (t_fin - t_ini) / dt / (N-1)
        ! allocate(output(size(rho0,1),size(rho0,2),min(N, int(dd))))
        output(:,:,1) = rho0

        iter = 0
        t = t_ini; rho = rho0
        do while (t < t_fin)

            do i=1,4
                k(:,:,i) = rhodot(t + rk4_t(i)*dt, rho + rk4_t(i)*k(:,:,i-1)*dt)
            enddo

            Delta = (1._dp/6._dp)*dt*(k(:,:,1) + 2*k(:,:,2) + 2*k(:,:,3) + k(:,:,4))

            rho = rho + Delta
            output(:,:,int(2+floor(iter/dd))) = rho
            iter = iter + 1
            t = t + dt
        end do

        ! deallocate(output)

    end function


end module integrators
! EOF
