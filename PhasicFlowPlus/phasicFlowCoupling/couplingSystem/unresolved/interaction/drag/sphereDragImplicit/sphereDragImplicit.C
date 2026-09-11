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

// Schiller & Naumann (1933) Z. Ver. Dtsch. Ing. 77, 318

template<typename DragClosureType>
pFlow::coupling::sphereDragImplicit<DragClosureType>::sphereDragImplicit
(
    const unresolvedCouplingSystem& uCS, 
    const porosity& 				prsty
)
:
	dragImplicit(uCS, prsty),
	dragClosure_(this->dict())
{}

template<typename DragClosureType>
void pFlow::coupling::sphereDragImplicit<DragClosureType>::calculateDragForceImplicit
(
    const fluidAveraging&           fluidVelocity,
    const solidAveraging&           parVelocity,
    const Plus::realx3ProcCMField&  parRotVelocity,
    const Plus::realProcCMField&    diameter,
    const distributionBase&         cellDistribution,    
    Plus::realx3ProcCMField&        particleForce,
    Plus::realProcCMField&          dragCoeff,
    Plus::realx3ProcCMField&        dragFluidVel
)
{
    setSuSpToZero();

    const auto* kernel = dynamic_cast<const particleGaussian*>(&cellDistribution);

    const auto& parCellInd =  this->parCellIndex();
    const auto& nu = this->mesh().template lookupObject<Foam::volScalarField>("nu");
    const auto& rho = this->mesh().template lookupObject<Foam::volScalarField>("rho");
    size_t numPar = parCellInd.size();
    const auto& alpha = this->alpha();
    
    const auto& Ufld =
        this->mesh().template lookupObject<Foam::volVectorField>("U");

    auto fluidVel = fluidVelocity.fieldSpan();
    auto solidVel = parVelocity.fieldSpan();

    auto& Su = this->Su();
    auto& Sp = this->Sp();

    auto pGradPtr = this->pressureGradient(rho);

    auto& quad = this->quadrature();

    if(this->faxen())
    {
        quad.resize(numPar);
        quad.resetCount();

        const auto& cm = this->Porosity().uCS().centerMass();

        #pragma omp parallel for schedule(dynamic)
        for(size_t i=0; i<numPar; ++i)
        {
            if(parCellInd[i] < 0) continue;
            quad.prepare
            (
                i,
                Foam::vector(cm[i].x(), cm[i].y(), cm[i].z()),
                0.5*diameter[i],
                parCellInd[i]
            );
        }
    }

    const auto& pGrad = this->sampleUndisturbed(pGradPtr());

    const auto& parCentre = this->Porosity().uCS().centerMass();
    const Foam::volVectorField& cellCentre = this->mesh().C();

    Foam::tmp<Foam::volScalarField> tYw;
    Foam::tmp<Foam::volVectorField> tNw;

    if(this->wallCorrection())
    {
        tYw = Foam::wallDist::New(this->mesh()).y();
        tNw = Foam::fvc::grad(tYw());
    }

    Foam::scalar wallFactorObs = 1.0;
    Foam::scalar wallReObs     = 0.0;
    Foam::scalar minHoverA     = Foam::GREAT;
    Foam::scalar wallBetaObs   = 1.0;
    Foam::scalar minHOverAEff  = Foam::GREAT;

    Foam::scalar betaObserved = 0.0;
    Foam::label  anyWidth     = 0;
    Foam::scalar faxenRatio   = 0.0;

    Plus::realProcCMField&   spPar = this->spScratch();
    Plus::realx3ProcCMField& upPar = this->upScratch();
    Plus::realx3ProcCMField& fbPar = this->fbScratch();
    std::fill(spPar.begin(), spPar.end(), 0.0);
    std::fill(upPar.begin(), upPar.end(), realx3(0.0));
    std::fill(fbPar.begin(), fbPar.end(), realx3(0.0));

    const bool cons = this->conservativeReaction();
    const bool pfr  = this->pointForceReaction();
    Plus::realx3ProcCMField& uePar = this->ueScratch();
    Plus::realx3ProcCMField& wxPar = this->wxScratch();
    Plus::realx3ProcCMField& ubPar = this->ubScratch();
    if(pfr)
    {
        std::fill(ubPar.begin(), ubPar.end(), realx3(0.0));
    }
    if(cons)
    {
        std::fill(uePar.begin(), uePar.end(), realx3(0.0));
        std::fill(wxPar.begin(), wxPar.end(), realx3(0.0));
        std::fill(ubPar.begin(), ubPar.end(), realx3(0.0));
    }

    this->bassetModel().beginCoupling(numPar);

    // pairInteraction
    if(this->pairModel().active())
    {
        if(!kernel)
        {
            fatalErrorInFunction
                << "pairInteraction needs distributionMethod particleGaussian "
                   "(the kernel width of each particle)." << endl;
            Plus::processor::abort(0);
        }

        this->pairModel().beginCoupling(numPar);

        for(size_t i=0; i<numPar; ++i)
        {
            const auto ci = parCellInd[i];
            if(ci < 0) continue;
            this->pairModel().publishWidth(i, kernel->kernelWidth(i), nu[ci]*rho[ci]);
        }

        this->pairModel().evaluate(this->faxen());
    }

    if(this->selfInducedLift()) this->resizeUndisturbed(numPar);

    // Gotoh (1990): rotation part of the local flow at each particle
    const Plus::procCMField<Foam::vector>* rotVort   = nullptr;
    const Plus::procCMField<Foam::vector>* rotStrain = nullptr;
    Foam::tmp<Foam::volVectorField> tRotVort;
    Foam::tmp<Foam::volVectorField> tRotStrain;
    if(this->rotationDrag())
    {
        tRotVort = Foam::fvc::curl(Ufld);
        tRotStrain = Foam::sqrt(2.0)*Foam::mag(Foam::symm(Foam::fvc::grad(Ufld)))
                   * Foam::dimensionedVector("ex", Foam::dimless, Foam::vector(1, 0, 0));
        rotVort   = &this->sampleRotationVorticity(tRotVort());
        rotStrain = &this->sampleRotationStrain(tRotStrain());
    }
    Foam::scalar rotFactorObs = 1.0;
    Foam::scalar rotTaObs     = 0.0;
    Foam::scalar rotPsiObs    = 1.0;
    Foam::scalar rotLatObs    = 0.0;

    // Faxen: Maxey & Riley (1983)
    if(this->faxen()) quad.bind(Ufld);

    #pragma omp parallel for schedule (dynamic) \
        reduction(max:betaObserved) reduction(max:anyWidth) \
        reduction(max:faxenRatio) reduction(max:wallFactorObs) \
        reduction(max:wallReObs) \
        reduction(min:minHoverA) reduction(min:wallBetaObs) \
        reduction(min:minHOverAEff) \
        reduction(max:rotFactorObs) reduction(max:rotTaObs) \
        reduction(min:rotPsiObs) reduction(max:rotLatObs)
    for(size_t parIndx=0; parIndx<numPar; parIndx++)
    {
        auto cellIndx = parCellInd[parIndx];

        if(cellIndx < 0 ) continue;

        auto rhoi = rho[cellIndx];
        auto mui = nu[cellIndx]* rhoi;
        Foam::scalar ef = alpha[cellIndx];
        auto dp = diameter[parIndx];
        auto vp =  Foam::constant::mathematical::pi/6 * Foam::pow(dp,3.0);

        Foam::vector up{solidVel[parIndx].x(), solidVel[parIndx].y(),solidVel[parIndx].z()};

        Foam::vector uEff = fluidVel[parIndx];
        
        if(this->faxen())
        {
            const Foam::vector uS = quad.surfaceAverage(parIndx, Ufld);
            faxenRatio = Foam::max
            (
                faxenRatio,
                Foam::mag(uS - uEff)/Foam::max(Foam::mag(uS), SMALL)
            );
            uEff = uS;
        }

        Foam::vector ur = uEff-up;
        
        Foam::scalar urMag = Foam::mag(ur);

        // Gotoh (1990): Omega = (|omega| - sigma)/2, Ta = Omega a^2/nu
        Foam::scalar omegaRot = 0.0;
        if(rotVort)
        {
            const Foam::scalar wMag  = Foam::mag((*rotVort)[parIndx]);
            const Foam::scalar sigma = (*rotStrain)[parIndx].x();
            omegaRot = 0.5*Foam::max(wMag - sigma, Foam::scalar(0));
        }
        const Foam::scalar nuP   = mui/rhoi;
        const Foam::scalar taRot = omegaRot*0.25*dp*dp/nuP;
        Foam::scalar fRot  = rotationDragFactorStokes(taRot);
        if(this->rotationDragGeostrophic())
        {
            // Stewartson (1953): 0.49 x 2 Omega rho V U over 6 pi mu a U,
            // weighted by Omega/(Omega + nu/a^2) = Ta/(1 + Ta)
            fRot += taRot/(1 + taRot)*(0.49*4.0/9.0)*taRot;
        }
        // drag factor max(f(Re), g), or 1 + sqrt((f - 1)^2 + (g - 1)^2) with
        // rotationDragCombination quadrature, applied to the closure's f(Re, eps)
        auto rotFactor = [&](const Foam::scalar ReX)
        {
            if(omegaRot <= 0) return Foam::scalar(1);
            const Foam::scalar f =
                Foam::max(dragClosure_.dimlessDrag(ReX, 1.0), Foam::SMALL);
            const Foam::scalar total = this->rotationDragQuadrature()
                ? 1 + Foam::sqrt(Foam::sqr(f - 1) + Foam::sqr(fRot - 1))
                : Foam::max(f, fRot);
            return Foam::max(Foam::scalar(1), total/f);
        };

        Foam::scalar hOverA   = Foam::GREAT;
        Foam::vector nw       = Foam::vector::zero;
        bool         haveWall = false;

        if(this->wallCorrection())
        {
            const Foam::scalar a = 0.5*dp;

            nw = tNw()[cellIndx];

            const Foam::scalar nMag = Foam::mag(nw);

            if(nMag > SMALL)
            {
                nw /= nMag;
                haveWall = true;
            }

            const Foam::vector xp
            (
                parCentre[parIndx].x(), parCentre[parIndx].y(), parCentre[parIndx].z()
            );
            Foam::scalar h = tYw()[cellIndx];
            if(haveWall) h += (nw & (xp - cellCentre[cellIndx]));

            hOverA = Foam::max(h, Foam::scalar(0))/Foam::max(a, SMALL);

            minHoverA = Foam::min(minHoverA, hOverA);

            hOverA = Foam::max(hOverA, this->minWallRadii());
        }

        const Foam::scalar sigmaK =
            (this->selfInduced_ && kernel) ? kernel->kernelWidth(parIndx) : 0.0;

        Foam::scalar beta = 0.0;

        if(sigmaK > 0)
        {
            anyWidth = 1;

            ef += selfInducedAlphaDeficit(selfInducedSigmaHat(dp, sigmaK));
            ef  = Foam::min(ef, Foam::scalar(1));

            const Foam::scalar Re0 = ef * rhoi * urMag * dp / mui;

            beta = selfInducedZetaU
            (
                this->selfInducedSampleFactor(), dp, sigmaK,
                Re0, dragClosure_.dimlessDrag(Re0, ef)
            );
            // rotation: the force per unit slip rises by rotFactor, and the far
            // field of the point force is screened by the stronger of the
            // inertial (Oseen) and rotational effects
            if(omegaRot > 0)
            {
                const Foam::scalar psiO = psiOseen(Re0*sigmaK/Foam::max(dp, SMALL));
                const Foam::scalar psiR = psiRotation(omegaRot, sigmaK, nuP);
                if(psiO > SMALL)
                {
                    const Foam::scalar psiRatio = Foam::min(psiO, psiR)/psiO;
                    beta *= rotFactor(Re0)*psiRatio;
                    rotPsiObs = Foam::min(rotPsiObs, psiRatio);
                }
            }
            if(haveWall)
            {
                const Foam::scalar h    = 0.5*dp*hOverA;
                const Foam::scalar aEff = selfInducedBlobRadius(sigmaK);

                const Foam::scalar hOverAEff = Foam::max
                (
                    h/Foam::max(aEff, SMALL),
                    this->minBlobWallRadii()
                );

                const Foam::scalar cosSW =
                    urMag > SMALL ? ((ur & nw)/urMag) : 0.0;

                Foam::scalar gWall =
                    selfInducedWallFactor(hOverAEff, cosSW);

                // gapReynolds
                if(this->wallGapRe())
                {
                    const Foam::scalar ReGapS =
                        Re0*0.5*(hOverA - 1.0)/this->wallGapReScale();
                    gWall = 1.0 - (1.0 - gWall)/(1.0 + Foam::max(ReGapS, Foam::scalar(0)));
                }

                wallBetaObs  = Foam::min(wallBetaObs, gWall);
                minHOverAEff = Foam::min(minHOverAEff, h/Foam::max(aEff, SMALL));

                beta *= gWall;
            }

            beta = Foam::min(beta, this->maxBeta_);
        }

        const Foam::scalar amp = 1.0/(1.0 - beta);

        // rotationSelfInducedLateral: the lateral part of the self-induced
        // velocity, beta_perp (omega_hat x w), with beta_perp = beta 0.0517 T/(1 + 0.524 T)
        // (Gotoh 1990; Candelier, Mehlig & Magnaudet 2019, eq. 4.10), removed from the
        // sampled slip; amp then removes the in-line part as before
        if(this->rotationSelfInducedLateral() && beta > 0 && omegaRot > 0 && rotVort && sigmaK > 0)
        {
            const Foam::vector wv    = (*rotVort)[parIndx];
            const Foam::scalar wvMag = Foam::mag(wv);
            if(wvMag > SMALL)
            {
                const Foam::scalar aS2 = Foam::constant::mathematical::pi*sigmaK*sigmaK;
                const Foam::scalar Ts  = Foam::sqrt(omegaRot*aS2/Foam::max(nuP, SMALL));
                const Foam::scalar betaPerp = beta*0.0517*Ts/(1.0 + 0.524*Ts);
                ur   -= (betaPerp/(1.0 - beta))*((wv/wvMag) ^ ur);
                urMag = Foam::mag(ur);
                uEff  = up + ur;
                rotLatObs = Foam::max(rotLatObs, betaPerp);
            }
        }

        // Esmaily & Horwitz (2018) eq. (2)
        if(this->selfInducedLift() && beta > 0)
        {
            const Foam::vector xc
            (
                parCentre[parIndx].x(), parCentre[parIndx].y(), parCentre[parIndx].z()
            );

            Foam::vector fL;

            if(this->liftFeedbackFor(parIndx, xc, 0.5*dp, fL))
            {
                const Foam::scalar Re0 = ef * rhoi * urMag * dp / mui;
                const Foam::scalar sp0 = 3 * Foam::constant::mathematical::pi
                                       * mui * ef * dp
                                       * dragClosure_.dimlessDrag(Re0, ef)
                                       * rotFactor(Re0);

                uEff += (beta/Foam::max(sp0, SMALL))*fL;
                ur    = uEff - up;
                urMag = Foam::mag(ur);
            }
        }

        // pairInteraction
        if(this->pairModel().active())
        {
            const Foam::vector dU = this->pairModel().du(parIndx);
            ur    += dU/amp;
            urMag  = Foam::mag(ur);
            uEff   = up + ur;
        }

        betaObserved = Foam::max(betaObserved, beta);

        const Foam::scalar Re = ef * rhoi * urMag*amp * dp /mui;

        Foam::scalar sp = 3 * Foam::constant::mathematical::pi * mui * ef * dp
                        * dragClosure_.dimlessDrag(Re, ef);

        {
            const Foam::scalar gRot = rotFactor(Re);
            sp *= gRot;
            rotFactorObs = Foam::max(rotFactorObs, gRot);
            rotTaObs     = Foam::max(rotTaObs, taRot);
        }

        sp *= amp;

        Foam::vector wallExtra = Foam::vector::zero;

        if(haveWall)
        {
            const Foam::scalar a = 0.5*dp;

            const bool fRe = this->wallFiniteRe();

            // Vasseur & Cox (1977) Fig. 4
            const Foam::scalar ReGap = Re*0.5*(hOverA - 1.0)/this->wallGapReScale();
            const Foam::scalar lamN =
                this->wallGapRe() ? wallDragNormalGapRe(hOverA, ReGap)
              : fRe ? wallDragNormalRe(hOverA, Re) : wallDragNormal(hOverA);
            const Foam::scalar lamP =
                fRe ? wallTransForceRe(hOverA, Re) : wallTransForce(hOverA);
            // Zeng et al. (2009); Goldman, Cox & Brenner (1967)
            const Foam::scalar lamS =
                fRe ? lamP : wallShearForce(hOverA);
            const Foam::scalar lamR = wallRotForce(hOverA);

            if(fRe && hOverA <= 9.0) wallReObs = Foam::max(wallReObs, Re);

            wallFactorObs = Foam::max(wallFactorObs, lamN);

            {
                const Foam::vector uPar = uEff - (uEff & nw)*nw;

                const Foam::vector wp
                (
                    parRotVelocity[parIndx].x(),
                    parRotVelocity[parIndx].y(),
                    parRotVelocity[parIndx].z()
                );

                const Foam::vector wPar = wp - (wp & nw)*nw;

                if(this->wallNormalImplicit())
                {
                    wallExtra = sp*(lamS - lamP)*uPar
                              + sp*a*lamR*(wPar ^ nw);
                }
                else
                {
                    wallExtra = sp*(lamN - lamP)*(ur & nw)*nw
                              + sp*(lamS - lamP)*uPar
                              + sp*a*lamR*(wPar ^ nw);
                }
            }

            if(this->wallNormalImplicit())
            {
                const Foam::scalar c2 = urMag > SMALL
                    ? Foam::sqr((ur & nw)/urMag) : Foam::scalar(0);
                sp *= lamP + (lamN - lamP)*c2;
            }
            else
            {
                sp *= lamP;
            }
        }

        Foam::vector pf = -vp*pGrad[parIndx] + wallExtra;
        
        particleForce[parIndx] += realx3(pf.x(), pf.y(), pf.z());

        // pairInteraction
        if(this->pairModel().active())
        {
            this->pairModel().addForce(parIndx, sp*(uEff - up) + wallExtra);
        }
        
        dragCoeff[parIndx]    += static_cast<real>(sp);
        dragFluidVel[parIndx]  = realx3
        (
            uEff.x(),
            uEff.y(),
            uEff.z()
        );
        
        
        spPar[parIndx]  = static_cast<real>(sp);
        upPar[parIndx]  = realx3(up.x(), up.y(), up.z());

        if(cons)
        {
            uePar[parIndx] = realx3(uEff.x(), uEff.y(), uEff.z());
            wxPar[parIndx] = realx3(wallExtra.x(), wallExtra.y(), wallExtra.z());
        }

        this->bassetModel().publish
        (
            parIndx,
            uEff,
            amp,
            1.5*dp*dp*Foam::sqrt(Foam::constant::mathematical::pi*rhoi*mui),
            dp,
            mui/rhoi
        );

        if(this->selfInducedLift())
        {
            this->setUndisturbed(parIndx, up + amp*(uEff - up));
        }
    }

    if(this->faxen()) quad.unbind();

    if(this->bassetActive())
    {
        const auto& runTime = this->mesh().time();

        this->bassetModel().evaluate(runTime.value(), runTime.deltaTValue());

        #pragma omp parallel for schedule (dynamic)
        for(size_t parIndx=0; parIndx<numPar; parIndx++)
        {
            if(parCellInd[parIndx] < 0) continue;

            const Foam::vector fE = this->bassetModel().explicitForce(parIndx);

            particleForce[parIndx] += realx3(fE.x(), fE.y(), fE.z());
            fbPar[parIndx]          = realx3(fE.x(), fE.y(), fE.z());
        }
    }

    semiImplicitCoupling(this->Porosity().uCS()).globalParticleSum(spPar);
    semiImplicitCoupling(this->Porosity().uCS()).globalParticleSum(upPar);

    if(this->bassetActive())
    {
        semiImplicitCoupling(this->Porosity().uCS()).globalParticleSum(fbPar);
    }

    Foam::scalar corrOverDrag = 0.0;
    if(cons)
    {
        semiImplicitCoupling(this->Porosity().uCS()).globalParticleSum(uePar);
        semiImplicitCoupling(this->Porosity().uCS()).globalParticleSum(wxPar);
    }

    if(cons || pfr)
    {
        for(size_t parIndx=0; parIndx<numPar; parIndx++)
        {
            if(spPar[parIndx] == 0) continue;
            if(!cellDistribution.hasSupport(parIndx) && parCellInd[parIndx] < 0) continue;

            Foam::vector ub = Foam::vector::zero;
            cellDistribution.inverseDistributeValue
            (
                parIndx, parCellInd[parIndx], Ufld.internalField(), ub
            );
            ubPar[parIndx] = realx3(ub.x(), ub.y(), ub.z());
        }

        semiImplicitCoupling(this->Porosity().uCS()).globalParticleSum(ubPar);
    }

    #pragma omp parallel for schedule (dynamic)
    for(size_t parIndx=0; parIndx<numPar; parIndx++)
    {
        const Foam::scalar sp = spPar[parIndx];

        const Foam::vector fb
        (
            fbPar[parIndx].x(), fbPar[parIndx].y(), fbPar[parIndx].z()
        );

        if(sp == 0 && Foam::magSqr(fb) == 0) continue;
        if(!cellDistribution.hasSupport(parIndx) && parCellInd[parIndx] < 0) continue;

        const Foam::vector up
        (
            upPar[parIndx].x(), upPar[parIndx].y(), upPar[parIndx].z()
        );

        Foam::vector corr = Foam::vector::zero;

        // point-force reaction: kernel average of U here, Sp*U subtracted
        // below, so the fluid gets -K(x)*sp*(Ubar - u_p) instead of
        // -K(x)*sp*(U(x) - u_p)
        Foam::vector dip = Foam::vector::zero;
        if(pfr)
        {
            dip = sp*Foam::vector(ubPar[parIndx].x(), ubPar[parIndx].y(), ubPar[parIndx].z());
        }

        if(cons)
        {
            const Foam::vector ue(uePar[parIndx].x(), uePar[parIndx].y(), uePar[parIndx].z());
            const Foam::vector wx(wxPar[parIndx].x(), wxPar[parIndx].y(), wxPar[parIndx].z());
            const Foam::vector ub(ubPar[parIndx].x(), ubPar[parIndx].y(), ubPar[parIndx].z());

            // Capecelatro & Desjardins (2013) eq. (28)
            corr = sp*(ub - ue) - wx;

            const Foam::scalar dragMag = Foam::mag(sp*(ue - up)) + Foam::mag(wx);
            if(dragMag > SMALL)
            {
                #pragma omp critical
                corrOverDrag = Foam::max(corrOverDrag, Foam::mag(corr)/dragMag);
            }
        }

        cellDistribution.distributeValue_OMP
        (
            parIndx, parCellInd[parIndx], Su, -(sp*up) + fb - corr + dip
        );
        cellDistribution.distributeValue_OMP
        (
            parIndx, parCellInd[parIndx], Sp, sp
        );
    }

    const auto& Vcells = this->mesh().V();

    forAll(Vcells, i)
    {
      Su[i] /= Vcells[i];
      Sp[i] /= Vcells[i];
    }

    cellDistribution.smoothenField(Sp);
    cellDistribution.smoothenField(Su);

    if(pfr)
    {
        // deferred correction: cancels the implicit Sp*U at coupling convergence
        forAll(Vcells, i)
        {
            Su[i] -= Sp[i]*Ufld.internalField()[i];
        }
    }

    this->reportSelfInduction(betaObserved, anyWidth);
    this->reportRotationDrag(rotFactorObs, rotTaObs, rotPsiObs, rotLatObs);

    {
        Foam::scalar dRep = 0, sRep = 0;

        for(size_t i=0; i<numPar; ++i)
        {
            if(parCellInd[i] < 0) continue;
            dRep = diameter[i];
            sRep = kernel ? kernel->kernelWidth(i) : 0.0;
            break;
        }

        Foam::reduce(dRep, Foam::maxOp<Foam::scalar>());
        Foam::reduce(sRep, Foam::maxOp<Foam::scalar>());

        this->reportDisturbanceRegime(dRep, sRep);

        if(this->bassetActive())
        {
            Foam::scalar slipRep = 0, nuRep = 0;

            for(size_t i=0; i<numPar; ++i)
            {
                if(parCellInd[i] < 0) continue;
                const Foam::vector uf(fluidVel[i].x(), fluidVel[i].y(), fluidVel[i].z());
                const Foam::vector vp(solidVel[i].x(), solidVel[i].y(), solidVel[i].z());
                slipRep = Foam::mag(uf - vp);
                nuRep   = nu[parCellInd[i]];
                break;
            }

            Foam::reduce(slipRep, Foam::maxOp<Foam::scalar>());
            Foam::reduce(nuRep,   Foam::maxOp<Foam::scalar>());

            this->bassetModel().reportKernelRegime
            (
                slipRep, dRep, nuRep, this->mesh().time().deltaTValue()
            );
        }
    }
    this->reportFaxen(faxenRatio);
    this->reportWall(wallFactorObs, minHoverA);
    this->reportConservative(corrOverDrag);
    this->reportWallRe(wallReObs);
    this->reportWallSelfInduction(wallBetaObs, minHOverAEff);

    if(this->faxen())
    {
        Foam::label nInc = quad.countIncomplete();
        Foam::reduce(nInc, Foam::sumOp<Foam::label>());

        if(nInc > 0 && Plus::processor::isMaster())
        {
            WARNING
                << nInc << " Faxen quadrature evaluation(s) this step could not "
                << "reach the whole particle: some nodes lie across a cyclic "
                << "patch, so the average was taken over the rest. Nodes on a "
                << "neighbouring rank are evaluated there, nodes inside a wall "
                << "take the wall value; expect this only for a particle within "
                << "a radius of a cyclic patch."
                << END_WARNING;
        }
    }

    Sp.correctBoundaryConditions();
	Su.correctBoundaryConditions();
    
}
