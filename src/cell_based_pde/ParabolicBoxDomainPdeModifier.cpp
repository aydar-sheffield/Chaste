/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "ParabolicBoxDomainPdeModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "TetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "AveragedSourceParabolicPde.hpp"


template<unsigned DIM>
ParabolicBoxDomainPdeModifier<DIM>::ParabolicBoxDomainPdeModifier()
    : AbstractBoxDomainPdeModifier<DIM>()
{
}

template<unsigned DIM>
ParabolicBoxDomainPdeModifier<DIM>::ParabolicBoxDomainPdeModifier(ParabolicPdeAndBoundaryConditions<DIM>* pPdeAndBcs,
                                                                  ChasteCuboid<DIM> meshCuboid,
                                                                  double stepSize)
    : AbstractBoxDomainPdeModifier<DIM>(),
      mpPdeAndBcs(pPdeAndBcs)
{
    assert(DIM == 2);

    //Generate mesh. Note only need to do this ones as the mesh is fixed.
    this->GenerateFeMesh(meshCuboid, stepSize);
}

template<unsigned DIM>
ParabolicBoxDomainPdeModifier<DIM>::~ParabolicBoxDomainPdeModifier()
{
    // It we have used this modifier, then we will have created a solution vector
    if (this->mSolution)
    {
        PetscTools::Destroy(this->mSolution);
    }
}

template<unsigned DIM>
void ParabolicBoxDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Set up boundary conditions
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer();

    this->UpdateCellPdeElementMap(rCellPopulation);

    // When using a PDE mesh which doesnt coincide with the cells, we must set up the source terms before solving the PDE.
    // Pass in already updated CellPdeElementMap to speed up finding cells.
    mpPdeAndBcs->SetUpSourceTermsForAveragedSourcePde(this->mpFeMesh, &this->mCellPdeElementMap);

    // Use SimpleLinearParabolicSolver as Averaged Source PDE
    SimpleLinearParabolicSolver<DIM,DIM> solver(this->mpFeMesh, mpPdeAndBcs->GetPde(), p_bcc.get());

    ///\todo Investigate more than one PDE time step per spatial step (#2687)
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();
    double dt = p_simulation_time->GetTimeStep();
    solver.SetTimes(current_time,current_time + dt);
    solver.SetTimeStep(dt);

    // Use previous solution as the initial condition
    Vec previous_solution = this->mSolution;
    solver.SetInitialCondition(previous_solution);

    // Note that the linear solver creates a vector, so we have to keep a handle on the old one
    // in order to destroy it.
    this->mSolution = solver.Solve();
    PetscTools::Destroy(previous_solution);
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ParabolicBoxDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractBoxDomainPdeModifier<DIM>::SetupSolve(rCellPopulation,outputDirectory);

    // Temporarily cache the variable name until we create an AbstractPdeAndBcs object
    // and move mpPdeAndBcs to the abstract class. See #2767
    this->mCachedDependentVariableName = mpPdeAndBcs->rGetDependentVariableName();

    // Cache the output directory
    this->mOutputDirectory = outputDirectory;

    // Copy the cell data to mSolution (this is the initial condition)
    SetupInitialSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ParabolicBoxDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer()
{
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    if (mpPdeAndBcs->IsNeumannBoundaryCondition())
    {
        // Impose any Neumann boundary conditions
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
             elem_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
             ++elem_iter)
        {
            p_bcc->AddNeumannBoundaryCondition(*elem_iter, mpPdeAndBcs->GetBoundaryCondition());
        }
    }
    else
    {
        // Impose any Dirichlet boundary conditions
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
                         node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
                         ++node_iter)
        {
            p_bcc->AddDirichletBoundaryCondition(*node_iter, mpPdeAndBcs->GetBoundaryCondition());
        }
    }

    return p_bcc;
}

template<unsigned DIM>
void ParabolicBoxDomainPdeModifier<DIM>::SetupInitialSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Specify Homogeneous initial conditons based upon the values stored in CellData.
    // Note need all the CellDataValues to be the same.

    double initial_condition = rCellPopulation.Begin()->GetCellData()->GetItem(mpPdeAndBcs->rGetDependentVariableName());

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        double initial_condition_at_cell = cell_iter->GetCellData()->GetItem(mpPdeAndBcs->rGetDependentVariableName());
        assert(fabs(initial_condition_at_cell - initial_condition)<1e-12);
    }

    // initialise mSolution
    this->mSolution = PetscTools::CreateAndSetVec(this->mpFeMesh->GetNumNodes(), initial_condition);
}

template<unsigned DIM>
void ParabolicBoxDomainPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractBoxDomainPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ParabolicBoxDomainPdeModifier<1>;
template class ParabolicBoxDomainPdeModifier<2>;
template class ParabolicBoxDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ParabolicBoxDomainPdeModifier)
