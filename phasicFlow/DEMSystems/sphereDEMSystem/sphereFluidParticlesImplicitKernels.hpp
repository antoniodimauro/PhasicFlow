/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): Antonio Di Mauro
  email: antoniodimauro03@gmail.com
------------------------------------------------------------------------------
Licence:
  This file is part of PhasicFlow, a CFD-DEM stack built on phasicFlow and
  PhasicFlowPlus (www.cemf.ir). It is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License v3 or
  any later version.

  It is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.
-----------------------------------------------------------------------------*/

#ifndef __sphereFluidParticlesImplicitKernels_hpp__
#define __sphereFluidParticlesImplicitKernels_hpp__

namespace pFlow::sphereFluidParticlesImplicitKernels
{

using rpImplicitDrag = Kokkos::RangePolicy<
			DefaultExecutionSpace,
			Kokkos::Schedule<Kokkos::Static>,
			Kokkos::IndexType<int32>
			>;

template<typename IncludeFunctionType>
void applyImplicitDragLift(
	real                            dt,
	deviceViewType1D<real>          mass,
	deviceViewType1D<real>          addedMass,
	deviceViewType1D<real>          dragCoeff,
	deviceViewType1D<realx3>        dragFluidVel,
	deviceViewType1D<realx3>        liftCrossVec,
	IncludeFunctionType             incld,
	deviceViewType1D<realx3>        velocity
)
{
	auto activeRange = incld.activeRange();
	if(incld.isAllActive())
	{
		Kokkos::parallel_for(
		"pFlow::sphereFluidParticlesKernels::applyImplicitDragLift",
		rpImplicitDrag(activeRange.first, activeRange.second),
		LAMBDA_HD(int32 i){
			const real    m   = mass[i] + addedMass[i];
			const real    A   = dt * dragCoeff[i] / m;
			const real    s   = 1 + A;
			const realx3  uf  = dragFluidVel[i];
			const realx3  B   = (dt / m) * liftCrossVec[i];

			const realx3  R   = velocity[i] + A * uf + cross(B, uf);
			const real    BR  = dot(B, R);
			const real    B2  = dot(B, B);
			const real    D   = s * s + B2;
			const realx3  BxR = cross(B, R);

			velocity[i] = (s * R + (BR / s) * B - BxR) / D;
		});
	}
	else
	{
		Kokkos::parallel_for(
		"pFlow::sphereFluidParticlesKernels::applyImplicitDragLift",
		rpImplicitDrag(activeRange.first, activeRange.second),
		LAMBDA_HD(int32 i){
			if(incld(i))
			{
				const real    m   = mass[i] + addedMass[i];
				const real    A   = dt * dragCoeff[i] / m;
				const real    s   = 1 + A;
				const realx3  uf  = dragFluidVel[i];
				const realx3  B   = (dt / m) * liftCrossVec[i];

				const realx3  R   = velocity[i] + A * uf + cross(B, uf);
				const real    BR  = dot(B, R);
				const real    B2  = dot(B, B);
				const real    D   = s * s + B2;
				const realx3  BxR = cross(B, R);

				velocity[i] = (s * R + (BR / s) * B - BxR) / D;
			}
		});
	}
	Kokkos::fence();
}

using rpAcceleration = Kokkos::RangePolicy<
			DefaultExecutionSpace,
			Kokkos::Schedule<Kokkos::Static>,
			Kokkos::IndexType<int32>
			>;

template<typename IncludeFunctionType>
void acceleration( 
	realx3		g,
	deviceViewType1D<real>  	mass,
	deviceViewType1D<real>  	addedMass,
	deviceViewType1D<real>  	displacedMass,
	deviceViewType1D<realx3>  	force,
	deviceViewType1D<realx3>  	fluidForce,
	deviceViewType1D<real>  	I,
	deviceViewType1D<real>  	addedInertia,
	deviceViewType1D<realx3>  	torque,
	deviceViewType1D<realx3>  	fluidTorque,
	IncludeFunctionType 		incld,
	deviceViewType1D<realx3> 	lAcc,
	deviceViewType1D<realx3> 	rAcc
	)
{

	auto activeRange = incld.activeRange();
	if(incld.isAllActive())
	{
		Kokkos::parallel_for(
		"pFlow::sphereParticlesKernels::acceleration",
		rpAcceleration(activeRange.first, activeRange.second),
		LAMBDA_HD(int32 i){
				const real me = mass[i] + addedMass[i];
				lAcc[i] = (force[i]+fluidForce[i])/me
						+ ((mass[i]-displacedMass[i])/me)*g;
				rAcc[i] = (torque[i]+fluidTorque[i])/(I[i]+addedInertia[i]);
		});
	}
	else
	{
		Kokkos::parallel_for(
		"pFlow::sphereParticlesKernels::acceleration",
		rpAcceleration(activeRange.first, activeRange.second),
		LAMBDA_HD(int32 i){
			if(incld(i))
			{
				const real me = mass[i] + addedMass[i];
				lAcc[i] = (force[i]+fluidForce[i])/me
						+ ((mass[i]-displacedMass[i])/me)*g;
				rAcc[i] = (torque[i]+fluidTorque[i])/(I[i]+addedInertia[i]);
			}
		});

	}
	
	Kokkos::fence();
}

}

namespace pFlow::sphereFluidParticlesImplicitKernels
{

template<typename IncludeFunctionType>
void applyImplicitRotDrag(
	real                            dt,
	deviceViewType1D<real>          I,
	deviceViewType1D<real>          addedInertia,
	deviceViewType1D<real>          rotDragCoeff,
	IncludeFunctionType             incld,
	deviceViewType1D<realx3>        rVelocity
)
{
	auto activeRange = incld.activeRange();
	if(incld.isAllActive())
	{
		Kokkos::parallel_for(
		"pFlow::sphereFluidParticlesKernels::applyImplicitRotDrag",
		rpImplicitDrag(activeRange.first, activeRange.second),
		LAMBDA_HD(int32 i){
			const real f = dt * rotDragCoeff[i] / (I[i] + addedInertia[i]);
			rVelocity[i] = rVelocity[i] / (1 + f);
		});
	}
	else
	{
		Kokkos::parallel_for(
		"pFlow::sphereFluidParticlesKernels::applyImplicitRotDrag",
		rpImplicitDrag(activeRange.first, activeRange.second),
		LAMBDA_HD(int32 i){
			if(incld(i))
			{
				const real f = dt * rotDragCoeff[i] / (I[i] + addedInertia[i]);
				rVelocity[i] = rVelocity[i] / (1 + f);
			}
		});
	}
	Kokkos::fence();
}

template<typename IncludeFunctionType>
void correctFusedAM4(
	real                     dt,
	real                     dampingVel,
	realx3                   g,
	deviceViewType1D<real>   mass,
	deviceViewType1D<real>   addedMass,
	deviceViewType1D<real>   displacedMass,
	deviceViewType1D<real>   I,
	deviceViewType1D<real>   addedInertia,
	deviceViewType1D<realx3> contactForce,
	deviceViewType1D<realx3> fluidForce,
	deviceViewType1D<realx3> contactTorque,
	deviceViewType1D<realx3> fluidTorque,
	deviceViewType1D<realx3> acceleration,
	deviceViewType1D<realx3> rAcceleration,
	deviceViewType1D<realx3> position,
	deviceViewType1D<realx3> velocity,
	deviceViewType1D<realx3> rVelocity,
	deviceViewType1D<realx3> pDy0, deviceViewType1D<realx3> pDy1,
	deviceViewType1D<realx3> pDy2, deviceViewType1D<realx3> pDy3,
	deviceViewType1D<realx3> vDy0, deviceViewType1D<realx3> vDy1,
	deviceViewType1D<realx3> vDy2, deviceViewType1D<realx3> vDy3,
	deviceViewType1D<realx3> wDy0, deviceViewType1D<realx3> wDy1,
	deviceViewType1D<realx3> wDy2, deviceViewType1D<realx3> wDy3,
	deviceViewType1D<real>   dragCoeff,
	deviceViewType1D<realx3> dragFluidVel,
	deviceViewType1D<realx3> liftCrossVec,
	deviceViewType1D<real>   rotDragCoeff,
	IncludeFunctionType      incld
)
{
	const real c9  = static_cast<real>( 9.0/24.0);
	const real c36 = static_cast<real>(36.0/24.0);
	const real c54 = static_cast<real>(54.0/24.0);
	auto activeRange = incld.activeRange();
	Kokkos::parallel_for(
		"pFlow::sphereFluidParticlesKernels::correctFusedAM4",
		rpImplicitDrag(activeRange.first, activeRange.second),
		LAMBDA_HD(int32 i){
			if( incld(i) )
			{
				const real   me   = mass[i] + addedMass[i];
				const realx3 lAcc = (contactForce[i]+fluidForce[i])/me
					+ ((mass[i]-displacedMass[i])/me)*g;
				const real   Ie   = I[i] + addedInertia[i];
				const realx3 rAcc = (contactTorque[i]+fluidTorque[i])/Ie;
				acceleration[i]  = lAcc;
				rAcceleration[i] = rAcc;

				{
					const realx3 dn = velocity[i];
					position[i] += dt*( c9*dn - c36*pDy0[i] + c54*pDy1[i]
					                    - c36*pDy2[i] + c9*pDy3[i] );
					pDy3[i]=pDy2[i]; pDy2[i]=pDy1[i]; pDy1[i]=pDy0[i]; pDy0[i]=dn;
				}
				{
					const realx3 dn = lAcc;
					velocity[i] += dampingVel*dt*( c9*dn - c36*vDy0[i] + c54*vDy1[i]
					                    - c36*vDy2[i] + c9*vDy3[i] );
					vDy3[i]=vDy2[i]; vDy2[i]=vDy1[i]; vDy1[i]=vDy0[i]; vDy0[i]=dn;
				}
				{
					const realx3 dn = rAcc;
					rVelocity[i] += dt*( c9*dn - c36*wDy0[i] + c54*wDy1[i]
					                    - c36*wDy2[i] + c9*wDy3[i] );
					wDy3[i]=wDy2[i]; wDy2[i]=wDy1[i]; wDy1[i]=wDy0[i]; wDy0[i]=dn;
				}
				{
					const real   m  = me;
					const real   A  = dt*dragCoeff[i]/m;
					const real   s  = 1 + A;
					const realx3 uf = dragFluidVel[i];
					const realx3 B  = (dt/m)*liftCrossVec[i];
					const realx3 R  = velocity[i] + A*uf + cross(B,uf);
					const real   BR = dot(B,R);
					const real   B2 = dot(B,B);
					const real   D  = s*s + B2;
					const realx3 BxR= cross(B,R);
					velocity[i] = (s*R + (BR/s)*B - BxR)/D;
				}
				{
					const real f = dt*rotDragCoeff[i]/Ie;
					rVelocity[i] = rVelocity[i]/(1+f);
				}
			}
		});
	Kokkos::fence();
}

}

namespace pFlow::sphereFluidParticlesImplicitKernels
{

template<typename IncludeFunctionType>
void setVelocityFromFluidForNew(
	int64                           lastMaxId,
	deviceViewType1D<uint32>        id,
	deviceViewType1D<realx3>        dragFluidVel,
	IncludeFunctionType             incld,
	deviceViewType1D<realx3>        velocity
)
{
	auto activeRange = incld.activeRange();
	Kokkos::parallel_for(
		"pFlow::sphereFluidParticlesKernels::setVelocityFromFluidForNew",
		rpImplicitDrag(activeRange.first, activeRange.second),
		LAMBDA_HD(int32 i){
			if( incld(i) && static_cast<int64>(id[i]) > lastMaxId )
			{
				velocity[i] = dragFluidVel[i];
			}
		});
	Kokkos::fence();
}

template<typename IncludeFunctionType>
int64 maxActiveId(
	deviceViewType1D<uint32>        id,
	IncludeFunctionType             incld
)
{
	auto activeRange = incld.activeRange();
	if(activeRange.second <= activeRange.first) return -1;

	int64 result = -1;
	Kokkos::parallel_reduce(
		"pFlow::sphereFluidParticlesKernels::maxActiveId",
		rpImplicitDrag(activeRange.first, activeRange.second),
		LAMBDA_HD(int32 i, int64& lmax){
			if( incld(i) )
			{
				const int64 thisId = static_cast<int64>(id[i]);
				if( thisId > lmax ) lmax = thisId;
			}
		},
		Kokkos::Max<int64>(result));
	Kokkos::fence();
	return result < 0 ? -1 : result;
}

}

#endif 
