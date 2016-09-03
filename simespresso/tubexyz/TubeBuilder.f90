!
!  GetAtomXYZByTube.f90
!
!  Free-Format Fortran Source File 
!  Generated by PGI Visual Fortran(R)
!  10/1/2011 9:11:25 PM
!
module TubeBuilder
use class_NanoTube
use HexagonalStructure
use PrismicStructure
use HexagonalHolesStructure
contains

subroutine GetAtomsXYZByTube(tube,xyz)
	type (NanoTube) :: tube
	real , allocatable :: xyz(:,:,:,:)

	select case (tube%TubeType)
		case (1)
			call GetHexAtomsXYZ(tube,xyz)
		case (2)
			call GetPrismicAtomsXYZ(tube,xyz)
		case (3)
			call GetHexHolesAtomsXYZ(tube,xyz)
	end select
	
end subroutine GetAtomsXYZByTube

end module TubeBuilder