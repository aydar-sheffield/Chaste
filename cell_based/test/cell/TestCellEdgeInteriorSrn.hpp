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

#ifndef TESTCELLEDGEINTERIORSRN_HPP_
#define TESTCELLEDGEINTERIORSRN_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellAgesWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "AlwaysDivideCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"

/*
 * The next header file defines a simple subcellular reaction network model that includes the functionality
 * for solving each cell's Delta/Notch signalling ODE system at each time step, using information about neighbouring
 * cells through the {{{CellEdgeData}}} class.
 */
#include "DeltaNotchSrnEdgeModel.hpp"
#include "SrnCellModel.hpp"
#include "DeltaNotchCellEdgeTrackingModifier.hpp"

#include "DeltaNotchSrnInteriorModel.hpp"
#include "DeltaNotchEdgeInteriorTrackingModifier.hpp"

/**
 * The tests below are designed for the pure edge SRN case, and the case with both edge and
 * interior SRN
 */
class TestCellEdgeInteriorSrn: public AbstractCellBasedTestSuite
{
public:
    /**
     * Tests for pure edge based SRN
     */
    void TestDeltaNotchEdgeSrnCorrectBehaviour()
        {
            TS_ASSERT_THROWS_NOTHING(DeltaNotchSrnEdgeModel srn_model);

            // Create cell edge SRN with four edges
            auto p_cell_edge_srn_model = new SrnCellModel();
            for (int i = 0; i < 4; i++)
            {
                boost::shared_ptr<DeltaNotchSrnEdgeModel> p_delta_notch_edge_srn_model(new DeltaNotchSrnEdgeModel());

                // Create a vector of initial conditions
                std::vector<double> starter_conditions;
                starter_conditions.push_back(0.5);
                starter_conditions.push_back(0.5);
                p_delta_notch_edge_srn_model->SetInitialConditions(starter_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_delta_notch_edge_srn_model);
            }

            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_cell_edge_srn_model, false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            std::vector<double> neigbour_delta = {1.0, 1.0, 1.0, 1.0};
            std::vector<double> interior_delta(4);
            std::vector<double> interior_notch(4);
            p_cell->GetCellEdgeData()->SetItem("neighbour delta", neigbour_delta);
            p_cell->GetCellData()->SetItem("interior delta", 0);
            p_cell->GetCellData()->SetItem("interior notch", 0);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            for (unsigned i = 0; i < p_cell_edge_srn_model->GetNumEdgeSrn(); i++)
            {
                auto p_delta_notch_edge_srn_model = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_edge_srn_model->GetEdgeSrn(i));

                TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(), 0.5, 1e-4);
                TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(), 0.5, 1e-4);
            }

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_cell_edge_srn_model->SimulateToCurrentTime();
            }

            // Test converged to steady state
            for (unsigned i = 0; i < p_cell_edge_srn_model->GetNumEdgeSrn(); i++)
            {
                auto p_delta_notch_edge_srn_model = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_edge_srn_model->GetEdgeSrn(i));

                TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(), 0.9900, 1e-4);
                TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(), 0.0101, 1e-4);
            }
        }

        void TestDeltaNotchEdgeSrnCreateCopy()
        {
            int numEdges = 4;

            auto p_cell_edge_srn_model = new SrnCellModel();
            for (int i = 0; i < numEdges; i++)
            {
                boost::shared_ptr<DeltaNotchSrnEdgeModel> p_delta_notch_edge_srn_model(new DeltaNotchSrnEdgeModel());

                // Set ODE system
                std::vector<double> state_variables;
                state_variables.push_back(2.0);
                state_variables.push_back(3.0);
                p_delta_notch_edge_srn_model->SetOdeSystem(new DeltaNotchEdgeOdeSystem(state_variables));
                p_delta_notch_edge_srn_model->SetInitialConditions(state_variables);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_delta_notch_edge_srn_model);
            }

            // Create a copy
            SrnCellModel* p_cell_edge_srn_model2 = static_cast<SrnCellModel*> (p_cell_edge_srn_model->CreateSrnModel());

            for (int i = 0; i < numEdges; i++)
            {
                auto p_delta_notch_edge_srn_model = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_edge_srn_model2->GetEdgeSrn(i));
                // Check correct initializations
                TS_ASSERT_EQUALS(p_delta_notch_edge_srn_model->GetNotch(), 2.0);
                TS_ASSERT_EQUALS(p_delta_notch_edge_srn_model->GetDelta(), 3.0);
            }

            // Destroy models
            delete p_cell_edge_srn_model;
            delete p_cell_edge_srn_model2;
        }

        void TestArchiveDeltaNotchSrnModel()
        {
            int numEdges = 4;

            OutputFileHandler handler("archive", false);
            std::string archive_filename = handler.GetOutputDirectoryFullPath() + "delta_notch_edge_srn.arch";

            // Create an output archive
            {
                SimulationTime* p_simulation_time = SimulationTime::Instance();
                p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

                UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();

                // As usual, we archive via a pointer to the most abstract class possible
                AbstractSrnModel* p_srn_model = new SrnCellModel();
                for (int i = 0; i < numEdges; i++)
                {
                    MAKE_PTR(DeltaNotchSrnEdgeModel, p_delta_notch_edge_srn_model);
                    static_cast<SrnCellModel *>(p_srn_model)->AddEdgeSrnModel(p_delta_notch_edge_srn_model);
                }

                MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
                MAKE_PTR(TransitCellProliferativeType, p_transit_type);

                // We must create a cell to be able to initialise the cell SRN model's ODE system
                CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model));
                p_cell->SetCellProliferativeType(p_transit_type);
                p_cell->GetCellEdgeData()->SetItem("neighbour delta", std::vector<double>{10.0, 10.0, 10.0, 10.0});
                p_cell->GetCellData()->SetItem("interior delta", 5.0);
                p_cell->GetCellData()->SetItem("interior notch", 1.0);
                p_cell->InitialiseCellCycleModel();
                p_cell->InitialiseSrnModel();
                p_cell->SetBirthTime(0.0);

                std::ofstream ofs(archive_filename.c_str());
                boost::archive::text_oarchive output_arch(ofs);

                // Read neighbour/interior Delta from CellEdgeData
                for (int i = 0; i < numEdges; i++)
                {
                    auto p_delta_notch_edge_model
                    = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(static_cast<SrnCellModel*>(p_srn_model)->GetEdgeSrn(i));
                    p_delta_notch_edge_model->UpdateDeltaNotch();
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetNeighbouringDelta(), 10.0, 1e-12);
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetInteriorDelta(), 5.0, 1e-12);
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetInteriorNotch(), 1.0, 1e-12);
                }

                output_arch << p_srn_model;

                // Note that here, deletion of the cell-cycle model and SRN is handled by the cell destructor
                SimulationTime::Destroy();
            }

            {
                // We must set SimulationTime::mStartTime here to avoid tripping an assertion
                SimulationTime::Instance()->SetStartTime(0.0);

                AbstractSrnModel* p_srn_model;

                std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
                boost::archive::text_iarchive input_arch(ifs);

                input_arch >> p_srn_model;

                for (int i = 0; i < numEdges; i++)
                {
                    auto p_delta_notch_edge_model
                    = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(static_cast<SrnCellModel*>(p_srn_model)->GetEdgeSrn(i));
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetNeighbouringDelta(), 10.0, 1e-12);
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetInteriorDelta(), 5.0, 1e-12);
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetInteriorNotch(), 1.0, 1e-12);
                }

                delete p_srn_model;
            }
        }

        void TestVertexMeshCellEdgeSrnDivision()
        {
            /* First we create a regular vertex mesh. */
            HoneycombVertexMeshGenerator generator(2, 2);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            std::vector<CellPtr> cells;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
            {
                /* Initalise cell cycle */
                auto p_cc_model = new AlwaysDivideCellCycleModel();
                p_cc_model->SetDimension(2);

                /* Initialise edge based SRN */
                auto p_element = p_mesh->GetElement(elem_index);
                auto p_cell_edge_srn_model = new SrnCellModel();

                /* Gets the edges of the element and create an SRN for each edge */
                for (unsigned i = 0; i < p_element->GetNumEdges() ; i ++)
                {
                    std::vector<double> initial_conditions;

                    /* Initial concentration of delta and notch is the same */
                    initial_conditions.push_back(2.0);
                    initial_conditions.push_back(2.0);

                    MAKE_PTR(DeltaNotchSrnEdgeModel, p_srn_model);
                    p_srn_model->SetInitialConditions(initial_conditions);
                    p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
                }

                CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
                p_cell->SetCellProliferativeType(p_diff_type);
                p_cell->SetBirthTime(0.0);
                p_cell->InitialiseCellCycleModel();
                p_cell->InitialiseSrnModel();
                cells.push_back(p_cell);
            }

            /* Create the cell population */
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            /* Create an edge tracking modifier */
            MAKE_PTR(DeltaNotchCellEdgeTrackingModifier<2>, p_modifier);

            p_modifier->SetupSolve(cell_population,"testVertexMeshCellEdgeSrnDivision");

            /* Divide the 0th cell*/
            {
                auto p_cell = cell_population.GetCellUsingLocationIndex(0);
                p_cell->ReadyToDivide();
                auto p_new_cell = p_cell->Divide();
                cell_population.AddCell(p_new_cell, p_cell);
                cell_population.Update(true);
            }

            /* We should now have 5 cells after the divide*/
            TS_ASSERT_EQUALS(cell_population.GetNumAllCells(), 5);

            /* Check the 0th cell */
            {
                auto p_cell = cell_population.GetCellUsingLocationIndex(0);
                auto p_cell_edge_srn = static_cast<SrnCellModel*>(p_cell->GetSrnModel());
                TS_ASSERT_EQUALS(p_cell_edge_srn->GetNumEdgeSrn(), 5);
            }

            /* Check the 4th cell */
            {
                auto p_cell = cell_population.GetCellUsingLocationIndex(4);
                auto p_cell_edge_srn = static_cast<SrnCellModel*>(p_cell->GetSrnModel());
                TS_ASSERT_EQUALS(p_cell_edge_srn->GetNumEdgeSrn(), 5);
            }
        }

        /**
         * Tests with both interior and edge SRNs
         */
        void TestDeltaNotchEdgeInteriorSrnCorrectBehaviour()
        {
            TS_ASSERT_THROWS_NOTHING(DeltaNotchSrnEdgeModel srn_model);

            // Create cell edge SRN with four edges
            auto p_cell_srn_model = new SrnCellModel();
            for (int i = 0; i < 4; i++)
            {
                boost::shared_ptr<DeltaNotchSrnEdgeModel> p_delta_notch_edge_srn_model(new DeltaNotchSrnEdgeModel());

                // Create a vector of initial conditions
                std::vector<double> starter_conditions;
                starter_conditions.push_back(0.5);
                starter_conditions.push_back(0.5);
                p_delta_notch_edge_srn_model->SetInitialConditions(starter_conditions);
                p_cell_srn_model->AddEdgeSrnModel(p_delta_notch_edge_srn_model);
            }
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_delta_notch_interior_srn_model(new DeltaNotchSrnInteriorModel());

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(1.0);
            starter_conditions.push_back(1.0);
            p_delta_notch_interior_srn_model->SetInitialConditions(starter_conditions);
            p_cell_srn_model->SetInteriorSrnModel(p_delta_notch_interior_srn_model);
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_cell_srn_model, false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            std::vector<double> neigbour_delta = {1.0, 1.0, 1.0, 1.0};
            std::vector<double> interior_delta = {1.0, 1.0, 1.0, 1.0};
            std::vector<double> interior_notch = {1.0, 1.0, 1.0, 1.0};

            p_cell->GetCellEdgeData()->SetItem("neighbour delta", neigbour_delta);
            p_cell->GetCellData()->SetItem("interior delta", 1.0);
            p_cell->GetCellData()->SetItem("interior notch", 1.0);
            p_cell->GetCellData()->SetItem("total neighbour edge delta", 4.0);
            p_cell->GetCellData()->SetItem("total edge notch", 2.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            for (unsigned i = 0; i < p_cell_srn_model->GetNumEdgeSrn(); i++)
            {
                auto p_delta_notch_edge_srn_model = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_srn_model->GetEdgeSrn(i));

                TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(), 0.5, 1e-4);
                TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(), 0.5, 1e-4);
            }
            {
                auto p_interior_srn = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_srn_model->GetInteriorSrn());
                TS_ASSERT_DELTA(p_interior_srn->GetNotch(), 1.0, 1e-4);
                TS_ASSERT_DELTA(p_interior_srn->GetDelta(), 1.0, 1e-4);
            }
            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_cell_srn_model->SimulateToCurrentTime();
            }

            // Test convergence to the steady state
            for (unsigned i = 0; i < p_cell_srn_model->GetNumEdgeSrn(); i++)
            {
                auto p_delta_notch_edge_srn_model = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_srn_model->GetEdgeSrn(i));

                TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(), 1.0900, 1e-4);
                TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(), 0.1083, 1e-4);
            }
            auto p_interior_srn
            = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_srn_model->GetInteriorSrn());
            TS_ASSERT_DELTA(p_interior_srn->GetNotch(), 0.9085, 1e-4);
            TS_ASSERT_DELTA(p_interior_srn->GetDelta(), 0.0022, 1e-4);
        }

        void TestDeltaNotchEdgeInteriorSrnCreateCopy()
        {
            int numEdges = 4;

            auto p_cell_srn_model = new SrnCellModel();
            for (int i = 0; i < numEdges; i++)
            {
                boost::shared_ptr<DeltaNotchSrnEdgeModel> p_delta_notch_edge_srn_model(new DeltaNotchSrnEdgeModel());

                // Set ODE system
                std::vector<double> state_variables;
                state_variables.push_back(2.0);
                state_variables.push_back(3.0);
                p_delta_notch_edge_srn_model->SetOdeSystem(new DeltaNotchEdgeOdeSystem(state_variables));
                p_delta_notch_edge_srn_model->SetInitialConditions(state_variables);
                p_cell_srn_model->AddEdgeSrnModel(p_delta_notch_edge_srn_model);
            }
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_delta_notch_interior_srn_model(new DeltaNotchSrnInteriorModel());
            std::vector<double> state_variables;
            state_variables.push_back(5.0);
            state_variables.push_back(7.0);
            p_delta_notch_interior_srn_model->SetOdeSystem(new DeltaNotchInteriorOdeSystem(state_variables));
            p_delta_notch_interior_srn_model->SetInitialConditions(state_variables);
            p_cell_srn_model->SetInteriorSrnModel(p_delta_notch_interior_srn_model);
            // Create a copy
            SrnCellModel* p_cell_srn_model2 = static_cast<SrnCellModel*> (p_cell_srn_model->CreateSrnModel());

            for (int i = 0; i < numEdges; i++)
            {
                auto p_delta_notch_edge_srn_model = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_srn_model2->GetEdgeSrn(i));
                // Check correct initializations
                TS_ASSERT_EQUALS(p_delta_notch_edge_srn_model->GetNotch(), 2.0);
                TS_ASSERT_EQUALS(p_delta_notch_edge_srn_model->GetDelta(), 3.0);
            }
            auto p_interior_srn_model2 = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_srn_model2->GetInteriorSrn());
            TS_ASSERT_EQUALS(p_interior_srn_model2->GetNotch(), 5.0);
            TS_ASSERT_EQUALS(p_interior_srn_model2->GetDelta(), 7.0);
            // Destroy models
            delete p_cell_srn_model;
            delete p_cell_srn_model2;
        }

        void TestArchiveDeltaNotchEdgeInteriorSrnModel()
        {
            int numEdges = 4;

            OutputFileHandler handler("archive", false);
            std::string archive_filename = handler.GetOutputDirectoryFullPath() + "delta_notch_edge_srn.arch";

            // Create an output archive
            {
                SimulationTime* p_simulation_time = SimulationTime::Instance();
                p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

                UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();

                // As usual, we archive via a pointer to the most abstract class possible
                AbstractSrnModel* p_srn_model = new SrnCellModel();
                for (int i = 0; i < numEdges; i++)
                {
                    MAKE_PTR(DeltaNotchSrnEdgeModel, p_delta_notch_edge_srn_model);
                    static_cast<SrnCellModel *>(p_srn_model)->AddEdgeSrnModel(p_delta_notch_edge_srn_model);
                }
                MAKE_PTR(DeltaNotchSrnInteriorModel, p_interior_srn_model);
                static_cast<SrnCellModel*>(p_srn_model)->SetInteriorSrnModel(p_interior_srn_model);
                MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
                MAKE_PTR(TransitCellProliferativeType, p_transit_type);

                // We must create a cell to be able to initialise the cell SRN model's ODE system
                CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model));
                p_cell->SetCellProliferativeType(p_transit_type);
                p_cell->GetCellEdgeData()->SetItem("neighbour delta", std::vector<double>{10.0, 10.0, 10.0, 10.0});
                p_cell->GetCellData()->SetItem("interior delta", 5.0);
                p_cell->GetCellData()->SetItem("interior notch", 1.0);
                p_cell->GetCellData()->SetItem("total neighbour edge delta", 40.0);
                p_cell->GetCellData()->SetItem("total edge notch", 4.0);
                p_cell->InitialiseCellCycleModel();
                p_cell->InitialiseSrnModel();
                p_cell->SetBirthTime(0.0);

                std::ofstream ofs(archive_filename.c_str());
                boost::archive::text_oarchive output_arch(ofs);

                // Read neighbour/interior Delta from CellEdgeData
                for (int i = 0; i < numEdges; i++)
                {
                    auto p_delta_notch_edge_model
                    = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(static_cast<SrnCellModel*>(p_srn_model)->GetEdgeSrn(i));
                    p_delta_notch_edge_model->UpdateDeltaNotch();
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetNeighbouringDelta(), 10.0, 1e-12);
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetInteriorDelta(), 5.0, 1e-12);
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetInteriorNotch(), 1.0, 1e-12);
                }
                auto p_interior_srn
                = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(static_cast<SrnCellModel*>(p_srn_model)->GetInteriorSrn());
                p_interior_srn->UpdateDeltaNotch();
                TS_ASSERT_DELTA(p_interior_srn ->GetTotalEdgeDelta(), 40.0, 1e-12);
                TS_ASSERT_DELTA(p_interior_srn ->GetTotalEdgeNotch(), 4.0, 1e-12);

                output_arch << p_srn_model;

                // Note that here, deletion of the cell-cycle model and SRN is handled by the cell destructor
                SimulationTime::Destroy();
            }

            {
                // We must set SimulationTime::mStartTime here to avoid tripping an assertion
                SimulationTime::Instance()->SetStartTime(0.0);

                AbstractSrnModel* p_srn_model;

                std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
                boost::archive::text_iarchive input_arch(ifs);

                input_arch >> p_srn_model;

                for (int i = 0; i < numEdges; i++)
                {
                    auto p_delta_notch_edge_model
                    = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(static_cast<SrnCellModel*>(p_srn_model)->GetEdgeSrn(i));
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetNeighbouringDelta(), 10.0, 1e-12);
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetInteriorDelta(), 5.0, 1e-12);
                    TS_ASSERT_DELTA(p_delta_notch_edge_model->GetInteriorNotch(), 1.0, 1e-12);
                }
                /*auto p_interior_srn
                = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(static_cast<SrnCellModel*>(p_srn_model)->GetInteriorSrn());
                TS_ASSERT_DELTA(p_interior_srn ->GetTotalEdgeDelta(), 40.0, 1e-12);
                TS_ASSERT_DELTA(p_interior_srn ->GetTotalEdgeNotch(), 4.0, 1e-12);*/
                delete p_srn_model;
            }
        }

        void TestVertexMeshCellEdgeInteriorSrnSwap()
        {
            /* First we create a regular vertex mesh. */
            HoneycombVertexMeshGenerator generator(2, 2);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            std::vector<CellPtr> cells;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
            {
                /* Initalise cell cycle */
                auto p_cc_model = new NoCellCycleModel();
                p_cc_model->SetDimension(2);

                /* Initialise edge based SRN */
                auto p_element = p_mesh->GetElement(elem_index);
                auto p_cell_srn_model = new SrnCellModel();

                /* Gets the edges of the element and create an SRN for each edge */
                for (unsigned i = 0; i < p_element->GetNumEdges() ; i ++)
                {
                    std::vector<double> initial_conditions;

                    /* Initial concentration of delta and notch is the same */
                    initial_conditions.push_back(1.0);
                    initial_conditions.push_back(1.0);

                    MAKE_PTR(DeltaNotchSrnEdgeModel, p_srn_model);
                    p_srn_model->SetInitialConditions(initial_conditions);
                    p_cell_srn_model->AddEdgeSrnModel(p_srn_model);
                }
                MAKE_PTR(DeltaNotchSrnInteriorModel, p_interior_model);
                p_interior_model->SetInitialConditions(std::vector<double>(2,1.0));
                p_cell_srn_model->SetInteriorSrnModel(p_interior_model);

                CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_srn_model));
                p_cell->SetCellProliferativeType(p_diff_type);
                p_cell->SetBirthTime(0.0);
                p_cell->InitialiseCellCycleModel();
                p_cell->InitialiseSrnModel();
                cells.push_back(p_cell);
            }

            /* Create the cell population */
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            /* Create an edge tracking modifier */
            MAKE_PTR(DeltaNotchCellEdgeTrackingModifier<2>, p_modifier);
            p_modifier->UpdateCellData(cell_population);
            /* Force T1 swap on the first shared edge of the first element
             * This means that edge 0 and 1 of the element after the swap
             * correspond to the edges 0 and 2 of the old element
             */
            auto elem = p_mesh->GetElement(0);
            const unsigned int old_n_edges = elem->GetNumEdges();

            const unsigned int edge_index_swapped = 1;
            auto edge = elem->GetEdge(edge_index_swapped);
            p_mesh->IdentifySwapType(edge->GetNode(0), edge->GetNode(1));

            p_mesh->RemoveDeletedNodes();
            const unsigned int new_n_edges = elem->GetNumEdges();
            TS_ASSERT_EQUALS(old_n_edges, new_n_edges+1);

            auto edge_operations = cell_population.GetCellEdgeChangeOperations();
            std::vector<std::vector<unsigned int> > pop_edge_status(4);
            std::vector<std::vector<long int> > pop_edge_mapping(4);
            for (auto operation:edge_operations)
            {
                const unsigned int element = operation->GetElementIndex();
                pop_edge_mapping[element] = operation->GetNewEdges()->GetEdgesMapping();
                pop_edge_status[element] = operation->GetNewEdges()->GetEdgesStatus();
            }
            cell_population.Update(false);
            //p_modifier->UpdateCellData(cell_population);

            //Check if the srn quantities are correct in cell 1. First check edges and then interior
            auto cell_0_srn
            = static_cast<SrnCellModel*>(cell_population.GetCellUsingLocationIndex(0)->GetSrnModel());
            const unsigned int n_edge_srns = cell_0_srn->GetEdges().size();
            TS_ASSERT_EQUALS(n_edge_srns,5);
            TS_ASSERT_EQUALS(n_edge_srns, p_mesh->GetElement(0)->GetNumEdges());
            for (unsigned int i=0; i<pop_edge_mapping[0].size(); ++i)
            {
                auto p_delta_notch_edge_srn_model
                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(cell_0_srn->GetEdgeSrn(i));
                if (pop_edge_status[0][i]!=0)
                {
                    TS_ASSERT_EQUALS(pop_edge_status[0][i],3);
                    //Quarters go for each adjacent edge...
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(),1.25,1e-10);
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(),1.25,1e-10);
                    assert(i==0||i==1);
                }
                else
                {
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(),1.0,1e-10);
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(),1.0,1e-10);
                }
            }

            //.. and a half to the interior
            auto p_interior_cell_0 = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(cell_0_srn->GetInteriorSrn());
            TS_ASSERT_DELTA(p_interior_cell_0->GetDelta(), 1.5, 1e-10);
            TS_ASSERT_DELTA(p_interior_cell_0->GetNotch(), 1.5, 1e-10);

            //Check if the srn quantities are correct in cell 2
            auto cell_1_srn
            = static_cast<SrnCellModel*>(cell_population.GetCellUsingLocationIndex(1)->GetSrnModel());
            const unsigned int n_edge_srns_1 = cell_1_srn->GetEdges().size();
            TS_ASSERT_EQUALS(n_edge_srns_1,5);
            TS_ASSERT_EQUALS(n_edge_srns_1, p_mesh->GetElement(1)->GetNumEdges());
            for (unsigned int i=0; i<pop_edge_mapping[1].size(); ++i)
            {
                auto p_delta_notch_edge_srn_model
                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(cell_1_srn->GetEdgeSrn(i));
                if (pop_edge_status[1][i]!=0)
                {
                    TS_ASSERT_EQUALS(pop_edge_status[1][i],3);
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(),1.25,1e-10);
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(),1.25,1e-10);
                    assert(i==3||i==4);//edges
                }
                else
                {
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(),1.0,1e-10);
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(),1.0,1e-10);
                }
            }
            //Check if the interior srn quantities have been updated properly
            auto p_interior_cell_1 = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(cell_1_srn->GetInteriorSrn());
            TS_ASSERT_DELTA(p_interior_cell_1->GetDelta(), 1.5, 1e-10);
            TS_ASSERT_DELTA(p_interior_cell_1->GetNotch(), 1.5, 1e-10);

            //Check if the srn quantities are correct in cell 3
            auto cell_2_srn
            = static_cast<SrnCellModel*>(cell_population.GetCellUsingLocationIndex(2)->GetSrnModel());
            const unsigned int n_edge_srns_2 = cell_2_srn->GetEdges().size();
            TS_ASSERT_EQUALS(n_edge_srns_2,7);
            TS_ASSERT_EQUALS(n_edge_srns_2, p_mesh->GetElement(2)->GetNumEdges());
            for (unsigned int i=0; i<pop_edge_mapping[2].size(); ++i)
            {
                auto p_delta_notch_edge_srn_model
                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(cell_2_srn->GetEdgeSrn(i));
                if (pop_edge_status[2][i]!=0)
                {
                    TS_ASSERT_EQUALS(pop_edge_status[2][i],2);
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(),0.0,1e-10);
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(),0.0,1e-10);
                    assert(i==0);
                }
                else
                {
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(),1.0,1e-10);
                    TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(),1.0,1e-10);
                }
            }

            auto p_interior_cell_2 = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(cell_2_srn->GetInteriorSrn());
            TS_ASSERT_DELTA(p_interior_cell_2->GetDelta(), 1.0, 1e-10);
            TS_ASSERT_DELTA(p_interior_cell_2->GetNotch(), 1.0, 1e-10);

            //Check if the srn quantities are correct in cell 4
            auto cell_3_srn
            = static_cast<SrnCellModel*>(cell_population.GetCellUsingLocationIndex(3)->GetSrnModel());
            const unsigned int n_edge_srns_3 = cell_3_srn->GetEdges().size();
            TS_ASSERT_EQUALS(n_edge_srns_3,6);
            TS_ASSERT_EQUALS(n_edge_srns_3, p_mesh->GetElement(3)->GetNumEdges());
            //The fourth cell srns are unaffected
            TS_ASSERT(pop_edge_mapping[3].empty());
            TS_ASSERT(pop_edge_status[3].empty());

            auto p_interior_cell_3 = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(cell_3_srn->GetInteriorSrn());
            TS_ASSERT_DELTA(p_interior_cell_3->GetDelta(), 1.0, 1e-10);
            TS_ASSERT_DELTA(p_interior_cell_3->GetNotch(), 1.0, 1e-10);
        }

        void TestVertexMeshCellEdgeInteriorSrnDivision()
        {
            /* First we create a regular vertex mesh. */
            HoneycombVertexMeshGenerator generator(2, 2);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            std::vector<CellPtr> cells;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
            {
                /* Initalise cell cycle */
                auto p_cc_model = new AlwaysDivideCellCycleModel();
                p_cc_model->SetDimension(2);

                /* Initialise edge based SRN */
                auto p_element = p_mesh->GetElement(elem_index);
                auto p_cell_edge_srn_model = new SrnCellModel();

                /* Gets the edges of the element and create an SRN for each edge */
                for (unsigned i = 0; i < p_element->GetNumEdges() ; i ++)
                {
                    std::vector<double> initial_conditions;

                    /* Initial concentration of delta and notch is the same */
                    initial_conditions.push_back(2.0);
                    initial_conditions.push_back(2.0);

                    MAKE_PTR(DeltaNotchSrnEdgeModel, p_srn_model);
                    p_srn_model->SetInitialConditions(initial_conditions);
                    p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
                }
                MAKE_PTR(DeltaNotchSrnInteriorModel, p_interior_srn);
                p_interior_srn->SetInitialConditions(std::vector<double>(2, 4.0));
                p_cell_edge_srn_model->SetInteriorSrnModel(p_interior_srn);

                CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
                p_cell->SetCellProliferativeType(p_diff_type);
                p_cell->SetBirthTime(0.0);
                p_cell->InitialiseCellCycleModel();
                p_cell->InitialiseSrnModel();
                cells.push_back(p_cell);
            }

            /* Create the cell population */
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            /* Create an edge tracking modifier */
            MAKE_PTR(DeltaNotchEdgeInteriorTrackingModifier<2>, p_modifier);
            p_modifier->SetupSolve(cell_population,"testVertexMeshCellEdgeInteriorSrnDivision");

            /* Divide the 0th cell*/
            {
                auto p_cell = cell_population.GetCellUsingLocationIndex(0);
                p_cell->ReadyToDivide();
                auto p_new_cell = p_cell->Divide();
                cell_population.AddCell(p_new_cell, p_cell);
                cell_population.Update(true);
            }

            /* We should now have 5 cells after the divide*/
            TS_ASSERT_EQUALS(cell_population.GetNumAllCells(), 5);

            /* The 0th and 4th index cells should not have only 5 edges and concentration split */

            /* Check the 0th cell */
            {
                auto p_cell = cell_population.GetCellUsingLocationIndex(0);
                auto p_cell_edge_srn = static_cast<SrnCellModel*>(p_cell->GetSrnModel());
                TS_ASSERT(p_cell_edge_srn->GetInteriorSrn() != nullptr);
                TS_ASSERT_EQUALS(p_cell_edge_srn->GetNumEdgeSrn(), 5);
            }
            /* Check the 4th cell */
            {
                auto p_cell = cell_population.GetCellUsingLocationIndex(4);
                auto p_cell_edge_srn = static_cast<SrnCellModel*>(p_cell->GetSrnModel());
                TS_ASSERT(p_cell_edge_srn->GetInteriorSrn() != nullptr);
                TS_ASSERT_EQUALS(p_cell_edge_srn->GetNumEdgeSrn(), 5);
            }
            /* Check interior SRN variable split in half*/
            {
                auto p_cell_1 = cell_population.GetCellUsingLocationIndex(0);
                auto p_cell_2 = cell_population.GetCellUsingLocationIndex(4);
                auto p_cell_srn_1 = static_cast<SrnCellModel*>(p_cell_1->GetSrnModel());
                auto p_cell_srn_2 = static_cast<SrnCellModel*>(p_cell_2->GetSrnModel());
                TS_ASSERT(p_cell_srn_1->GetInteriorSrn() != nullptr);
                TS_ASSERT(p_cell_srn_2->GetInteriorSrn() != nullptr);
                auto interior_srn_1
                = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(static_cast<SrnCellModel*>(p_cell_srn_1)->GetInteriorSrn());
                auto interior_srn_2
                = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(static_cast<SrnCellModel*>(p_cell_srn_2)->GetInteriorSrn());
                const double delta_1 = interior_srn_1->GetDelta();
                const double delta_2 = interior_srn_2->GetDelta();
                const double notch_1 = interior_srn_1->GetNotch();
                const double notch_2 = interior_srn_2->GetNotch();
                TS_ASSERT_DELTA(delta_1, delta_2,1e-12);
                TS_ASSERT_DELTA(notch_1, notch_2,1e-12);
                TS_ASSERT_DELTA(delta_1,2.0,1e-12);
                TS_ASSERT_DELTA(notch_1,2.0,1e-12);
            }
        }


};

#endif //TESTCELLEDGEINTERIORSRN_HPP_