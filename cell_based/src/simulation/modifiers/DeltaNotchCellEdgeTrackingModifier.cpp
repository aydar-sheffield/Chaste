/*

Copyright (c) 2005-2019, University of Oxford.
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

#include "DeltaNotchCellEdgeTrackingModifier.hpp"
#include "SrnCellModel.hpp"
#include "DeltaNotchSrnInteriorModel.hpp"
#include "DeltaNotchSrnEdgeModel.hpp"

template<unsigned DIM>
DeltaNotchCellEdgeTrackingModifier<DIM>::DeltaNotchCellEdgeTrackingModifier()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
DeltaNotchCellEdgeTrackingModifier<DIM>::~DeltaNotchCellEdgeTrackingModifier()
{
}

template<unsigned DIM>
void DeltaNotchCellEdgeTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Update the cell
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DeltaNotchCellEdgeTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DeltaNotchCellEdgeTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Recovers each cell's edge levels of notch and delta, and those of its neighbor's
    // Then saves them
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        auto p_cell_edge_model = static_cast<SrnCellModel*>(cell_iter->GetSrnModel());

        /* Cells edge data */
        std::vector<double> notch_vec;
        std::vector<double> delta_vec;
        for (unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            boost::shared_ptr<DeltaNotchSrnEdgeModel> p_model
            = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_edge_model->GetEdgeSrn(i));
            double this_delta = p_model->GetDelta();
            double this_notch = p_model->GetNotch();
            delta_vec.push_back(this_delta);
            notch_vec.push_back(this_notch);
        }
        // Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        cell_iter->GetCellEdgeData()->SetItem("edge notch", notch_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge delta", delta_vec);

        /* Cell interior data */
        //We're not using interior model so interiors are set to 0
        cell_iter->GetCellData()->SetItem("interior delta", 0);
        cell_iter->GetCellData()->SetItem("interior notch", 0);
    }

    //After the edge data is filled, fill the edge neighbour data
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
    {
        auto p_cell_edge_model = static_cast<SrnCellModel*>(cell_iter->GetSrnModel());
        const unsigned int n_cell_edges = p_cell_edge_model->GetNumEdgeSrn();
        std::vector<double> neigh_mean_delta(n_cell_edges);

        for (unsigned int i=0; i<n_cell_edges; ++i)
        {
            //Get neighbouring cell's values of delta on this
            auto elemNeighbours = rCellPopulation.GetNeighbouringEdgeIndices(*cell_iter, i);
            double mean_delta = 0;
            for (auto neighbourIndex: elemNeighbours)
            {
                auto neighbourCell = rCellPopulation.GetCellUsingLocationIndex(neighbourIndex.first);
                std::vector<double> neighbour_delta_vec = neighbourCell->GetCellEdgeData()->GetItem("edge delta");
                mean_delta += neighbour_delta_vec[neighbourIndex.second];
            }
            if (elemNeighbours.size()>0)
                mean_delta = mean_delta/elemNeighbours.size();
            neigh_mean_delta[i] = mean_delta;
        }
        cell_iter->GetCellEdgeData()->SetItem("neighbour delta", neigh_mean_delta);
    }

}

template<unsigned DIM>
void DeltaNotchCellEdgeTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class DeltaNotchCellEdgeTrackingModifier<1>;
template class DeltaNotchCellEdgeTrackingModifier<2>;
template class DeltaNotchCellEdgeTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeltaNotchCellEdgeTrackingModifier)
