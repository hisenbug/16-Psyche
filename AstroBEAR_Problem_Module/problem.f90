!#########################################################################
!
!    Copyright (C) 2003-2012 Department of Physics and Astronomy,
!                            University of Rochester,
!                            Rochester, NY
!
!    problem.f90 of module Template is part of AstroBEAR.
!
!    AstroBEAR is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AstroBEAR is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AstroBEAR.  If not, see <http://www.gnu.org/licenses/>.
!
!#########################################################################
!> @dir Template
!! @brief Contains files necessary for the Template Calculation

!> @file problem.f90
!! @brief Main file for module Problem

!> @defgroup Template Template Module
!! @brief Module for calculating collapse of a uniform cloud
!! @ingroup Modules

!> Template Module
!! @ingroup Template
MODULE Problem
  USE GlobalDeclarations
  USE DataDeclarations
  USE Clumps
  USE Ambients
  USE Fields
  USE Refinements
  USE Shapes
  USE Clumps
  USE ProcessingDeclarations
  USE Winds
  USE Tests
  USE TemplateControl
  USE Outflows
  USE InternalBoundary
  USE PhysicsDeclarations
  USE SweepDeclarations
  IMPLICIT NONE
  SAVE
  PUBLIC ProblemModuleInit, ProblemGridInit, ProblemBeforeStep, &
       ProblemAfterStep, ProblemSetErrFlag, ProblemBeforeGlobalStep

  PRIVATE

  REAL(KIND=qPREC) :: pos(3), rPlanet, rhoBoundary, tBoundary, rhoPlanet, tPlanet, rhoAmb, tAmb, rhoWind, tWind, betaWind, machWind, Bdir(3), bwind, vwind
  REAL(KIND=qPREC) :: bodyRes=1d6, ambRes=1d-3, shellFrac=0.1, shellResFrac=0.1, coreFrac=0d0, coreRes=1d-5, criticalRes=-1.0, rExcl=0d0
  INTEGER :: refFlag = 1, wAxis = 1, rwFlag = 0, RESRAMP = 3, excludeBField = 0
  REAL(KIND=qPREC) :: wDir(3) = (/1d0,0d0,0d0/)
  REAL(KIND=qPREC) :: tRamp = 1.224d0, vAmb, delta_vel=1d10, delta_T, delta_rho, bAmb(3), shellRes, OKres, diff_alpha2_0, rampstart, ramptime
  REAL(KIND=qPREC) :: resSItoEtaComp = 9.74028491071365d0, softZone=0d0, rAtm=0d0, rP_b, psyB=0d0
  TYPE(InternalBoundaryDef), POINTER :: lInternalBoundary
  TYPE(WindDef), POINTER :: Wind
  INTEGER :: expPoints = 1	!! TURN OFF BEFORE RUNNING CODE NORMALLY

CONTAINS

  SUBROUTINE ProblemModuleInit() BIND(C)
    TYPE(AmbientDef), POINTER :: Ambient
    TYPE(ClumpDef), POINTER :: Planet
    TYPE(RefinementDef), POINTER :: Refinement			! For testing AMR (1 Aug'19)
    INTEGER :: itrace
    REAL(KIND=qPREC) ::cs, vxnorm=-1, vxnormR=1, bigR
    real, pointer :: a
    NAMELIST/ProblemData/ pos, rPlanet, rhoBoundary, tBoundary, rhoAmb, tAmb, rhoWind, tWind, bWind, Bdir, vwind, vxnorm, vxnormR, tRamp, &
	  & bodyRes, ambRes, shellFrac, shellResFrac, coreFrac, criticalRes, refFlag, wDir, wAxis, bAmb, rwFlag, rhoPlanet, tPlanet, vAmb, rExcl, &
	  & diff_alpha2_0, rampstart, ramptime, softZone, rAtm, rP_b

    OPEN(UNIT=PROBLEM_DATA_HANDLE, FILE='problem.data', STATUS="OLD")
    READ(PROBLEM_DATA_HANDLE,NML=ProblemData)
    CLOSE(PROBLEM_DATA_HANDLE)


	psyB = acos(bDir(2))
    rhoWind=rhoWind/nScale
    tWind=tWind/TempScale
    bWind=bWind*1d-9*tesla/BScale
    vWind=vWind*1d5/velScale
    vAmb=vAmb*1d5/velScale
    bAmb=bAmb*1d-9*tesla/BScale

    rhoAmb=rhoAmb/nScale
    tAmb=tAmb/TempScale
	rhoBoundary=rhoBoundary/nScale
	tBoundary=tBoundary/TempScale

    rhoPlanet=rhoPlanet/rScale
    tPlanet=tPlanet/TempScale
    rPlanet=rPlanet*1d5/lScale

    if(criticalRes <= 0) then
    	criticalRes = vWind*2d0*rPlanet	!! Crude Calculation
    endif

    bigR = rPlanet
    if(rPlanet<rP_b) bigR=rP_b

    bodyRes = bodyRes*resSItoEtaComp			!! Convert
    shellRes = shellResFrac*bodyRes
    OKres = criticalRes*0.2d0							!! 2% of critRes is sufficient
    delta_vel = vWind/tRamp ! gamma*tWind/(2d0*rPlanet) !! Want to ramp up to full in tRamp computational time
	delta_T = (tWind - tBoundary)/tRamp
	delta_rho = (rhoWind - rhoBoundary)/tRamp

    resistivity_function => ProblemResistivity

    CALL CreateAmbient(Ambient, rhoAmb, rhoAmb*tAmb)
    Ambient%velocity=vAmb*wDir
    Ambient%B = bAmb
    ! if (all(bdir == [0,1,0]) .or. all(bdir == [1,0,0]) .or. all(bdir == [0,0,1])) excludeBField = 1
    if((rExcl > 0 ) .or. rwFlag == 1) Ambient%B = [0d0,0d0,0d0]		! .and. excludeBField == 1  !! Always exclude field
    ! Make 0 Ambient for case with velocity Ramping up or Excluded Clump
    CALL UpdateAmbient(Ambient)

    CALL CreateWind(Wind, rhoWind, tWind, vWind)
    Wind%dir=wAxis			! Change to 1, 2, 3 (x, y, z)
    Wind%edge=1				! Left or Right edge
    Wind%B=bWind*Bdir
    CALL UpdateWind(Wind)

    if(rAtm > 1d0 .and. rAtm < levels(MaxLevel)%dx) rAtm = levels(MaxLevel)%dx

!    CALL CreateClump(Planet, rhoPlanet, rhoPlanet*tPlanet, pos)
!    Planet%radius=rPlanet
!    CALL UpdateClump(Planet)

	if(refFlag /= 0) then
	!! Right Hemisphere
		CALL CreateInternalBoundary(lInternalBoundary)
		lInternalBoundary%Shape%type=ELLIPSOID	! WAS SPHERE
		lInternalBoundary%Shape%position=pos
		lInternalBoundary%Shape%size_param=[rPlanet,rPlanet,rP_b]

		!lInternalBoundary%Shape%isinshape_userdef => ProblemShape

		CALL SetShapeOrientation(lInternalBoundary%Shape, 0d0, 0.9553166181245092d0, 3.14159265358979323846264d0/4d0) !creates rotation matrix
		lInternalBoundary%density=rhoBoundary
		lInternalBoundary%pressure=rhoBoundary*tBoundary
		lInternalBoundary%type=REFLECTING
		lInternalBoundary%vnormfactor=vxnormR
        lInternalBoundary%MaxLevel=MaxLevel
		CALL UpdateInternalBoundary(lInternalBoundary)

		!! Refine region behind body for reconnection
		CALL CreateRefinement(Refinement)
		CALL CreateShape(Refinement%Shape)
		Refinement%Shape%type=RECTANGULAR_PRISM
		Refinement%Shape%position=pos + bigR*(/3d0,0d0,0d0/)
		Refinement%Shape%size_param=bigR*(/1.9d0,1d0,1d0/)
		Refinement%MaxLevel = MaxLevel - 1
		Refinement%type=REFINE_INSIDE
		CALL SetShapeBounds(Refinement%Shape)

		!! Refine region in sphere around body
		if (.false.) then
			CALL CreateRefinement(Refinement)
			CALL CreateShape(Refinement%Shape)
			Refinement%Shape%type=RECTANGULAR_PRISM
			Refinement%Shape%position=pos
			Refinement%Shape%size_param=bigR*(/1.1d0,1.1d0,1.1d0/)
			Refinement%MaxLevel = MaxLevel
			Refinement%type=REFINE_INSIDE
			CALL SetShapeBounds(Refinement%Shape)
		endif

	!! MAKE LV. Max at interface (Failed completely)
		if (.false.) then
			CALL CreateRefinement(Refinement)
			CALL CreateShape(Refinement%Shape)
			Refinement%Shape%type=SPHERE
			Refinement%Shape%position=pos
			Refinement%Shape%size_param=(8d0*levels(MaxLevel)%dx+rPlanet)*(/1d0,1d0,1d0/)
			Refinement%MaxLevel = MaxLevel
			Refinement%type=GRADIENT
			Refinement%field=UserResistivity_Field
			Refinement%scale=LogScale
			CALL SetShapeBounds(Refinement%Shape)
		endif

    else
	!! For testing AMR Stability on 1st Aug 2019
		CALL CreateRefinement(Refinement)
		CALL CreateShape(Refinement%Shape)
		Refinement%Shape%type=SPHERE
		Refinement%Shape%position=pos ! + rPlanet*(/0d0,2.5d0,0d0/)
		Refinement%Shape%size_param=rPlanet*(/1d0,1d0,1d0/)
		! Refinement%Shape%vel=(/1d0,0d0,0d0/) ! Doesn't work
		CALL SetShapeBounds(Refinement%Shape)
		Refinement%type=REFINE_INSIDE
		 !Refinement%BufferCells=4
	endif
	!! For AMR

    IF (MPI_ID == 0) then
    	write(*,*) 'Beta = ', 2*rhoWind*tWind/bwind**2
    	write(*,*) 'Mach = ', vwind/sqrt(gamma*tWind)
    	write(*,*) 'Critical Resistivty (approx in computational scale) = ', criticalRes
    	Write(*,*) 'Body Resistivty (approx in computational scale) = ', bodyRes
    	write(*,*) 'Divergence Warning is turned off.'
	ENDIF

    CALL AddTracer(itrace,'tracer')

    CALL AddDiagnosticVar(UserResistivity_Field)
    CALL AddDiagnosticVar(ErrFlag_Field)

  END SUBROUTINE ProblemModuleInit

  SUBROUTINE ProblemGridInit(Info) BIND(C)
    TYPE(InfoDef) :: Info

    DIFF_ALPHA2 = diff_alpha2_0

    if (refflag /= 0 .and. rwFlag /= 1 .and. rExcl > 0) then	!! For Excluding Clump Field and Velocity
       call applyfunction(ambientVelocity, Info)  !! Apply velocity
       IF (lMHD) then	! .and. excludeBField == 1 !! Always want to exclude field
          call ApplyEMF(uniformExcludedPlanet, Info, lPeriodic_opt=(/.false.,.false.,.false./),subsample_opt=2**(MaxLevel-Info%level))
       end IF
    end if

    write(*,*) "Grid Init Was Called."	! For testing 3D AMR Bug, Aug'19
!!	call OutputDoubleArray(Info%q(0:Info%mx(1)+1,0:Info%mx(2)+1,1,iBx))
!!	call OutputDoubleArray(Info%q(0:Info%mx(1)+1,0:Info%mx(2)+1,1,iBy))
!!	call OutputDoubleArray(Info%q(0:Info%mx(1)+1,0:Info%mx(2)+1,1,iBz))
  END SUBROUTINE ProblemGridInit


  SUBROUTINE AmbientMagneticField(cellpos, t, dt, q, Ghost)	!! Better to introduce magnetic field here?
    REAL(KIND=qPREC), DIMENSION(:), INTENT(IN) :: cellpos
    REAL(KIND=qPREC), DIMENSION(:), INTENT(INOUT) :: q
    REAL(KIND=qPREC), INTENT(IN) :: t, dt
    LOGICAL, INTENT(IN) :: Ghost   !! Use?
    if (cellpos(1) < pos(1)-rplanet) then
       q(1)=rhoWind
       q(2)=get_velocity(t)
       if (iE /= 0) then
          q(iE)=rhoWind*twind
       end if
       if (lMHD) then
          !! q((/iBz,iBy,iBx/))=f(cellpos((/3,2,1/)))
          q(iBx:iBz)=bwind*bdir
       end if
    elseif (cellpos(1) > pos(1)+rplanet) then
       if (nDim == 1) then
          q(ivx)=0
       end if
    end if
  end SUBROUTINE AmbientMagneticField

  SUBROUTINE AmbientVelocity(cellpos, t, dt, q, Ghost)
    REAL(KIND=qPREC), DIMENSION(:), INTENT(IN) :: cellpos
    REAL(KIND=qPREC), DIMENSION(:), INTENT(INOUT) :: q
    REAL(KIND=qPREC), INTENT(IN) :: t, dt
    REAL(KIND=qPREC) :: r, tmpV(3), rvec(3)
    LOGICAL, INTENT(IN) :: Ghost   !! Use?

    if (rwFlag == 1) then ! Ramping up solar wind
    	q(2:4) = [0d0,0d0,0d0]
    	return
	endif

    if (excludeBField /= 1) then
    	if(r < rExcl) then
    		q(2:4) = [0d0,0d0,0d0]
    	else
    		q(2:4) = vAmb*wDir
    	endif
    	return
	endif

    r = sqrt(sum((cellpos-pos)**2))
    rvec = cellpos-pos

    if (vAmb /= 0) then		!! Why do I need this if?
    	if(r < rExcl) then
    		q(2:4) = [0d0,0d0,0d0]
		else
			!if(wAxis == 1) then
			!	tmpV = [cellpos(1)-pos(1),0d0,0d0]
			!elseif(wAxis == 2) then
			!	tmpV = [0d0,cellpos(2)-pos(2),0d0]
			!elseif(wAxis == 3) then
			!	tmpV = [0d0,0d0,cellpos(3)-pos(3)]
			!endif
       		!q(2:4) = vAmb*( wDir - 0.5d0*rPlanet**3*( DOT_PRODUCT(cellpos-pos,tmpV) *(cellpos-pos)/r**3 &
       		!	& - wDir ) /r**3)
       		!! Optimized for when wDir == (1,0,0)
   			q(2) = vAmb*(1d0-(rExcl/r)**2+(rvec(2)**2+rvec(3)**2)*rExcl**2/r**4)
   			q(3) = -vAmb*rvec(2)*rvec(1)*rExcl**2/r**4
   			q(4) = -vAmb*rvec(3)*rvec(1)*rExcl**2/r**4
       	endif
    end if
  end SUBROUTINE AmbientVelocity

  ! Initial function
  ! q is vector potential in computational units
  SUBROUTINE AmbientPotential(cellpos, t, dt, q)
    REAL(KIND=qPREC), DIMENSION(:), INTENT(IN) :: cellpos
    REAL(KIND=qPREC), DIMENSION(:), INTENT(INOUT) :: q
    REAL(KIND=qPREC), INTENT(IN) :: t, dt
    REAL(KIND=qPREC) :: Az
    Az=0
    if (cellpos(1) < pos(1)-rplanet) then
       Az=-(bWind-bamb(2))*(cellpos(1)-(pos(1)-rplanet))
    end if
    if (all(bdir == (/0,1,0/))) then
       q=(/0d0,0d0,Az/)
    elseif (all(bdir == (/0,0,1/))) then
       q=(/0d0,-Az,0d0/)
    end if
  END SUBROUTINE AmbientPotential


  SUBROUTINE ProblemBeforeStep(Info) BIND(C)
    TYPE(InfoDef) :: Info
  END SUBROUTINE ProblemBeforeStep

  SUBROUTINE ProblemAfterStep(Info) BIND(C)
    TYPE(InfoDef) :: Info
  END SUBROUTINE ProblemAfterStep

  SUBROUTINE ProblemSetErrFlag(Info) BIND(C)
    TYPE(InfoDef) :: Info
  END SUBROUTINE ProblemSetErrFlag

  ! Function to exclude field from the planet interior
  SUBROUTINE uniformExcludedPlanet(cellpos, t, dt, q)  !!pos - PlanetPos, PlanetRadius, StellarMagRadius, StellarB
    REAL(KIND=qPREC), DIMENSION(:), INTENT(IN) :: cellpos
    REAL(KIND=qPREC), DIMENSION(:), INTENT(INOUT) :: q
    REAL(KIND=qPREC), INTENT(IN) :: t, dt
    ! LOGICAL, INTENT(IN) :: Ghost   !! Use?

    REAL(KIND=qPREC), DIMENSION(3) :: rvec, thetavec, phivec
    REAL(KIND=qPREC) :: r, theta, phi, c1, s1
    REAL(KIND=qPREC) :: rotMat(3,3)=RESHAPE((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),(/3,3/))    !Rotates lab frame into object frame

	c1=cos(psyB)
    !c2=cos(theta_)
    !c3=cos(phi_)
    s1=sin(psyB)
    !s2=sin(theta_)
    !s3=sin(phi_)

    !Shape%RotationMatrix=RESHAPE((/ &
    !     c1*c2*c3-s1*s3, -c3*c2*s1-c1*s3, +c3*s2,&
    !     c1*c2*s3+c3*s1, c1*c3-c2*s1*s3, s2*s3, &
    !     -c1*s2, s1*s2, c2/),(/3,3/),(/0d0/),(/2,1/))
    rotMat=RESHAPE((/ &
         1d0, 0d0, 0d0,&
         0d0, c1, -s1, &
         0d0, s1, c1/),(/3,3/),(/0d0/),(/2,1/))

	rvec = MATMUL(rotMat,(cellpos-pos))
    r=sqrt(sum((cellpos-pos)**2))
    ! B = B0*((stellarRadius + sqrt(sum(stellarPos**2)))/sqrt(sum(planetPos**2)))**2

    if (r > rExcl) then
      ! Modified so field is along x axis

      if (all(bdir == [0,1,0])) then
      	q = bWind*(rExcl**3-r**3)*[-rvec(3),0d0,rvec(1)]/(2*r**3)
      elseif (all(bdir == [1,0,0])) then
      	q = bWind*(rExcl**3-r**3)*[0d0,rvec(3),-rvec(2)]/(2*r**3)
      elseif (all(bdir == [0,0,1])) then
      	q = bWind*(rExcl**3-r**3)*[rvec(2),-rvec(1),0d0]/(2*r**3)
      else
		  theta=acos(rvec(2)/r)
		  phi=atan2(rvec(1),rvec(3))
		  !thetavec = MATMUL(rotMat,[-sin(theta), cos(theta)*sin(phi), cos(theta)*cos(phi)])
		  phivec = MATMUL(TRANSPOSE(rotMat),[cos(phi), 0d0, -sin(phi)])	! rotMat needs to be transposed.

		  ! A (Vector Potential)
		  q = bWind*(rExcl**3-r**3)*sin(theta)*phivec/(2*r**2)
      endif
    else
      q = [0d0,0d0,0d0]
    end if
  END SUBROUTINE uniformExcludedPlanet

  !! Do not call this during normal run!
  subroutine ExportPoints(Info)
  	integer :: i,j,k
  	real(kind=qPrec) :: psn(3),r,rPos(3)
  	type(infodef), pointer :: info

  	do i=1,Info%mX(1)
  		do j=1,Info%mX(2)
  			do k=1,Info%mX(3)
  				psn=cellpos(Info,i,j,k)

  				rPos=TransformCoordsToShape(lInternalBoundary%Shape, psn)
       			r=sqrt(sum(rPos(:)**2/lInternalBoundary%Shape%size_param(:)**2))
       			if (r>=1d0-1d0*levels(MaxLevel)%dx/MINVAL(lInternalBoundary%Shape%size_param(:)) .and. &
       				r<=1d0+1d0*levels(MaxLevel)%dx/MINVAL(lInternalBoundary%Shape%size_param(:))) &
  					write(MPI_ID+3000,'(6E25.15)') psn, info%q(i,j,k,iBx:iBz)
  			enddo
  		enddo
  	enddo
  end subroutine ExportPoints

  SUBROUTINE ProblemBeforeGlobalStep(n)		 BIND(C)! Modify time variable parameter here
    INTEGER :: n, ierr
    REAL(KIND=qPREC) :: temp

	TYPE(NodeDef), POINTER :: node
	TYPE(NodeDefList), POINTER :: nodelist
	type(infodef), pointer :: info

	!! Needs to be disabled to make code work
	!! Run on a single CPU
    if (expPoints==1 .and. n==MaxLevel) then
		nodelist=>Nodes(n)%p
		DO WHILE (associated(nodelist))
		 node=>nodelist%self
		 info=>node%info
		 call ExportPoints(info)
		 nodelist=>nodelist%next
		END DO
		close(unit=MPI_ID+1000)
		call MPI_Barrier(MPI_COMM_WORLD,ierr)
		stop
	endif

    DIFF_ALPHA2 = diff_alpha2_0

    if (.false. .and. n == 0 .and. ramptime > 100d0) then
       DIFF_ALPHA2 = diff_alpha2_0 * min(1d0, max(0d0,(ramptime - (levels(n)%tnow - rampstart))/ramptime))
    endif

	!! Ramping up solar wind
    if (rwFlag == 1) then
    	Wind%velocity=get_velocity(levels(n)%tnow)

	!! Ramping down planet rho and P. Check!! See delta definitions above!
	elseif (rwFlag == 2 .and. levels(n)%tnow <= tRamp) then
    	lInternalBoundary%density = rhoWind - levels(n)%tnow*delta_rho
		lInternalBoundary%pressure = lInternalBoundary%density*(tWind - levels(n)%tnow*delta_T)
		CALL UpdateInternalBoundary(lInternalBoundary)
	! elseif (rwFlag == RESRAMP .and. levels(n)%tnow <= tRamp) then ! ProblemResistivity takes care of time

	end if
  END SUBROUTINE ProblemBeforeGlobalStep

  FUNCTION get_velocity(t)  !! For ramping up wind velocity
  	real(kind=qPREC) :: get_velocity
  	real(kind=qPREC), INTENT(IN) :: t
  	get_velocity = min(vWind, delta_vel*t)
  	!! defined delta_vel = gamma*tWind/(2d0*rPlanet) in comp units for 1 mach/per crossing time (TOO LONG!)
  end function

  function ProblemShape(cellpos)
	logical :: ProblemShape
	real(kind=qprec), intent(IN), dimension(:) :: cellpos
	ProblemShape = .false.
	! Define the geometry of the shape.
  end function ProblemShape

  function ProblemResistivity(cellpos, q, t)
    real(kind=qPREC), dimension(3) :: ProblemResistivity
    real(kind=qPREC), dimension(3), INTENT(IN) :: cellpos
    real(kind=qPREC), dimension(:), INTENT(IN) :: q
    real(kind=qPREC), INTENT(IN) :: t
    real(kind=qPREC) :: res(3), rPos(3)
    real(kind=qPREC) :: r, tmpR, shRes

    if(lResistive  == .false.) then
    	ProblemResistivity=ambres
    	return
	endif

	if (shellFrac < 1d-4 .and. refFlag /= 11) shellRes = bodyRes !! For fixed Res case
    shRes = shellRes
    if (rwFlag == RESRAMP .and. t < tRamp) shRes = OKres-(OKres-shellRes)*t/tRamp

    r_soft = softZone*levels(MaxLevel)%dx + rPlanet + rAtm
	!r_soft = 2d0*(GxBounds(1,2)-GxBounds(1,1)) / (GlInternalBoundarymX(1)) + rPlanet
    !r_soft = 2d0*(GxBounds(1,2)-GxBounds(1,1)) / ((2d0**3)*GmX(1)) + rPlanet  ! Assuming 3 levels. For softening resistivity

    if (refFlag == 0) then
       res=ambres
	else
       rPos=TransformCoordsToShape(lInternalBoundary%Shape, cellpos, t)
       r=sqrt(sum((rPos)**2))

       if (r < coreFrac*rP_B) then
          res = coreRes
       elseif (SUM((rpos(:))**2/lInternalBoundary%Shape%size_param(:)**2) <= 1d0) then
          res = bodyRes
       ! Can add a 3rd layer which is ellipse
       !elseif (r < r_soft) then
       !   res = shRes * ( (ambRes/shRes)**( (r-rPlanet)/(r_soft-rPlanet) ) )	!! Doesn't implement properly for no shell (fixed).
       ! elseif (nDim == 1 .and. cellpos(1) > pos(1)+rPlanet) then
       !   res=criticalRes   ! Should be criticalRes ?
       else
          res=ambRes
       end if
    end if
    ProblemResistivity=res
  end function ProblemResistivity

END MODULE Problem

