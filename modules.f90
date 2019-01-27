module system
real*8 delta   ! delta is the discretization lenght in z direction
integer  dimz  ! number of lattice sites in z direction
real*8 sigma
real*8 csalt
real*8 pKa
real*8 pHbulk
real*8 st
endmodule

module results
real*8, allocatable :: xflag(:)
real*8, allocatable :: xh(:)
real*8, allocatable :: psi(:)
real*8, allocatable :: xpos(:)
real*8, allocatable :: xneg(:)
real*8, allocatable :: xHplus(:)
real*8, allocatable :: xOHmin(:)
real*8, allocatable :: qtot(:)
real*8, allocatable :: pro(:)
real*8 q
end subroutine

module kai
integer Xulimit ! cutoff for poor sv interaction in lattice sites
real*8, allocatable :: Xu(:) ! poor-sv coefficients
endmodule

module chainsdat
integer cuantas  ! number of polyelectrolyte conformations
integer long     ! lenght of polyelectrolyte chain
real*8 lseg
integer*1, allocatable :: pz(:,:)
endmodule

module const
real*8, parameter :: pi = 3.14159 ! pi 
real*8, parameter :: Na = 6.02d23 ! Avogadro's number
real*8, parameter :: lb = 0.714   ! bjerrum lenght in water in nm
real*8 constq
real*8 pKw
endmodule


module solver
integer infile
real*8, parameter :: error = 1.0d-6 ! maximum kinsol norm
endmodule


module molecules
real*8 zpos, zneg, zpol ! charges of cation, anions and polyelectrolyte segment
real*8 vsalt, vpol      ! volume of salt and polyelectrolyte segments in units of vsol
real*8 vsal             ! solvent volume 
real*8 K0
endmodule

module bulk
real*8 xHplusbulk, xOHminbulk ! bulk volume fraction of H+ and OH-
real*8 xposbulk, xnegbulk     ! bulk volume fraction of cation and anion
real*8 xsolbulk               ! bulk volume fraction of solvent
real*8 expmupos, expmuneg, expmuHplus, expmuOHmin  ! exp(-beta*mu)*(bulk volume fraction), where mu is the chemical potential
endmodule

module rans
integer seed
endmodule


