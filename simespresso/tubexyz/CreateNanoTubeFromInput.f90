!
!  constructNanoTubeFromInput.f90
!
!  Free-Format Fortran Source File 
!  Generated by PGI Visual Fortran(R)
!  9/27/2011 9:37:15 PM
!
module CreateNanoTubeFromInput
use class_NanoTube
implicit none

contains

	subroutine Welcome()
		
			print "(A)","==== Welcome to tubeXYZ ===="
			print "(A)","==== Nanotube generator for solid state research ===="
			print "(A)","====================================================="
		
	end subroutine Welcome

	function HowManyNanoTubes() result(tubeTotal)
		integer :: tubeTotal

		integer :: choice

		do
			print "(A)","In this software you can construct multiple nanotubes in the same XYZ file"
			print "(A)","The tubes will be constructed one by one from the input you provide"
			print "(A/)","Please choose the mode of operation"
			print "(A)","1. construct one nanotube (default) "
			print "(A)","2. construct multiple nanotubes in the same file"


			read (*,"(I1)"), choice
			
			if (choice .eq. 1 .or. choice .eq. 0 ) then
					tubeTotal=1
					exit
			end if

			if (choice .eq. 2 ) then
				print "(A/)","Set the number of tubes to construct"
				read *, tubeTotal

				if (tubeTotal .gt. 0 ) exit
				
			end if

		end do




	end function HowManyNanoTubes

	function EnterTube() result(tube)
		type (NanoTube)::tube

		integer :: choiceForType=0,choice=0
		integer :: choiceForN = 0,ChoiceForM=0
		real :: bondLength = 0

		integer :: tubeLength = 0

	

		do
			print "(A/)","Choose the nanotube type"
			print "(A)","1. Carbon nanotube (default) "
			print "(A)","2. Silicon nanotube"
			print "(A)","3. Boron nanotube"

			read (*,"(I1)"),choiceForType
			
			if (choiceForType .eq. 0 ) choiceForType = 1

			if (choiceForType <=4) then
				exit
			end if

			print "(A/)","Please select 1-4"
			
		end do

		do
			print "(A)","Choose the nanotube chiral vector n,m default (5,5) "

			read (*,"(I3,I3)"),choiceForN,ChoiceForM

			if (choiceForN .gt. 0 .and. choiceForM .ge. 0) then
				exit
			end if

			if (choiceForN .eq. 0 .and. choiceForM .eq. 0) then
				choiceForN = 5
				choiceForM = 5
				exit
			end if 

			print "(A)", "Please enter the vector correctly (integers only)"
		end do

		tube = make_NanoTube(choiceForType,choiceForN,ChoiceForM)

		do
			print "(A/)","Choose the bond length (angstrom):"
			print "(A,F4.2,A)","1. Use default bond length : ",tube%Sigma, " (default)"
			print "(A)","2. Set the bond length"

			read (*,"(I1)"),choice

			if (choice .eq. 0 ) choice = 1

			if (choice.eq.2) then
				print "(A/)","Enter the bond length for the tube (angstrom)"
				read *,bondLength
					if (bondLength > 0 ) then
						exit
					end if
			end if
			
			if (choice .eq. 1) exit
			
			print "(A/)", "Please select 1 or 2"
		end do

		do
			print "(A/)","Choose the nanotube length in atom layers units"
			print "(A/)","this will determine how many helix turns the program will generate"
			print "(A)","1. Use default length of 30 (default)"
			print "(A)","2. Set the length"

			read (*,"(I1)"),choice

			if (choice .eq. 0 ) choice = 1

			if (choice.eq.2) then
				print "(A)","Enter the tube length"
				read *,tubeLength
					if (tubeLength > 0 ) then
						exit
					end if
			end if

			if (choice .eq. 1) exit
			
		end do

		if (bondLength > 0 ) then
			tube = SetBondLength(tube,bondLength)
		end if

		if (tubeLength > 0) then
			tube = SetLength(tube,tubeLength)
		end if

		end function

		function ChooseWhatKindOfFileToOutput() result(fileName)

			character (len=200) :: fileName

			print "(A)","We are now ready to write the tubes to a data file"
			print "(A)","An XYZ file will be constructed to be compatible with Aviz"
			print "(A)","An additional file with the same name but with an extension of xyz.info"
			print "(A)/","Will be constructed and will contain the tubes data for future reference"
			print "(A)","Please set the filename for the output data (default = tubes.xyz)"

			read "(A)",fileName

			if (fileName .eq. "" ) fileName = "tubes.xyz"

		end function ChooseWhatKindOfFileToOutput


		


end module