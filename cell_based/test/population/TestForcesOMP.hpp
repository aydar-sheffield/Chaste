//
// Created by twin on 09/09/19.
//

#ifndef TESTFORCESOMP_HPP_
#define TESTFORCESOMP_HPP_

#include <chrono>
#include <cxxtest/TestSuite.h>
#include <AbstractCellBasedTestSuite.hpp>

#include "CheckpointArchiveTypes.hpp"



#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

#include "FarhadifarForce.hpp"


#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"

#include "PetscSetupAndFinalize.hpp"


class TestForcesOMP : public AbstractCellBasedTestSuite
{
public:

    void TestFarhadifarForceOMP()
    {
        #ifndef _OPENMP
        TS_FAIL("OpenMP not used");
        return;
        #endif

        std::vector<c_vector<double, 2> > applied_forces(6);
        {
            /**
             * Here we test that the forces are applied correctly to individual nodes.
             * We apply the force to something like this:
             *  . ____ . ____ .
             *  |      |      |
             *  |      |      |
             *  . ____ . ____ .
             */
            std::vector<Node<2>*> nodes;
            // the boolean says wether the node is a boundary node or not
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
            nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
            nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
            nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
            nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

            // make two square elements out of these nodes
            std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
            unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
            unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

            for (unsigned i=0; i<4; i++)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }

            std::vector<VertexElement<2,2>*> vertex_elements;
            vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
            vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

            // Make a vertex mesh
            MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

            TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
            TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

            // Get a cell population
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
            VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

            // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                cell_iter->SetBirthTime(-5.0);
            }

            MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
            p_growth_modifier->UpdateTargetAreas(cell_population);

            // Now let's make a FarhadifarForce and apply it to the population
            FarhadifarForce<2> force;

            force.AddForceContribution(cell_population);

            for (unsigned int i=0; i<applied_forces.size(); ++i)
            {
                applied_forces[i] = cell_population.rGetMesh().GetNode(i)->rGetAppliedForce();
            }

            c_vector<double, 2> applied_force_0;
            c_vector<double, 2> applied_force_1;
            applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
            applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

            // If this is a Farhadifar force, this will be the force at the vertices
            TS_ASSERT_DELTA(applied_force_0[0], 3.44, 1e-10);
            TS_ASSERT_DELTA(applied_force_0[1], 3.44, 1e-10);
            TS_ASSERT_DELTA(applied_force_1[0], 0.0, 1e-10);
            TS_ASSERT_DELTA(applied_force_1[1], 6.76, 1e-10);
        }
        std::vector<c_vector<double, 2> > applied_forces_OMP(6);
        {
            /**
             * Here we test that the forces are applied correctly to individual nodes.
             * We apply the force to something like this:
             *  . ____ . ____ .
             *  |      |      |
             *  |      |      |
             *  . ____ . ____ .
             */
            std::vector<Node<2>*> nodes;
            // the boolean says wether the node is a boundary node or not
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
            nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
            nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
            nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
            nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

            // make two square elements out of these nodes
            std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
            unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
            unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

            for (unsigned i=0; i<4; i++)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }

            std::vector<VertexElement<2,2>*> vertex_elements;
            vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
            vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

            // Make a vertex mesh
            MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

            TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
            TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

            // Get a cell population
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
            VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

            // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                cell_iter->SetBirthTime(-5.0);
            }

            MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
            p_growth_modifier->UpdateTargetAreas(cell_population);

            // Now let's make a FarhadifarForce and apply it to the population
            FarhadifarForce<2> force;

            force.AddForceContribution(cell_population);

            for (unsigned int i=0; i<applied_forces.size(); ++i)
            {
                applied_forces_OMP[i] = cell_population.rGetMesh().GetNode(i)->rGetAppliedForce();
            }
            c_vector<double, 2> applied_force_0;
            c_vector<double, 2> applied_force_1;
            applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
            applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

            // If this is a Farhadifar force, this will be the force at the vertices
            TS_ASSERT_DELTA(applied_force_0[0], 3.44, 1e-10);
            TS_ASSERT_DELTA(applied_force_0[1], 3.44, 1e-10);
            TS_ASSERT_DELTA(applied_force_1[0], 0.0, 1e-10);
            TS_ASSERT_DELTA(applied_force_1[1], 6.76, 1e-10);
        }
        for (unsigned int i=0; i<applied_forces.size(); ++i)
        {
            TS_ASSERT_DELTA(applied_forces[i][0], applied_forces_OMP[i][0], 1e-12);
            TS_ASSERT_DELTA(applied_forces[i][1], applied_forces_OMP[i][1], 1e-12);
        }
    }

    void TestFarhadifarForceOMPPopulationSpeedUp()
    {
        #ifndef _OPENMP
        TS_FAIL("OpenMP not used");
        return;
        #endif

        const unsigned int n_x_cells = 100;
        const unsigned int n_y_cells = 100;
        std::vector<c_vector<double,2> > applied_force;
        std::vector<c_vector<double,2> > applied_force_OMP;
        double wall_time, cpu_time, wall_time_OMP, cpu_time_OMP;
        unsigned int n_population_nodes;

        std::vector<c_vector<double, 2> > applied_forces;
        std::vector<c_vector<double, 2> > applied_forces_OMP;

        {

            HoneycombVertexMeshGenerator generator(n_x_cells, n_y_cells);
            MutableVertexMesh<2,2>* vertex_mesh = generator.GetMesh();
            // Get a cell population
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            cells_generator.GenerateBasic(cells, vertex_mesh->GetNumElements(), std::vector<unsigned>());
            VertexBasedCellPopulation<2> cell_population(*vertex_mesh, cells);
            n_population_nodes = vertex_mesh->GetNumNodes();
            // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                cell_iter->SetBirthTime(-5.0);
            }

            MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
            p_growth_modifier->UpdateTargetAreas(cell_population);

            // Now let's make a FarhadifarForce and apply it to the population
            FarhadifarForce<2> force;
            std::clock_t c_start = std::clock();
            auto t_start = std::chrono::high_resolution_clock::now();
            force.AddForceContribution(cell_population);
            std::clock_t c_end = std::clock();
            auto t_end = std::chrono::high_resolution_clock::now();
            for (unsigned int i=0; i<n_population_nodes; ++i)
            {
                applied_forces.push_back(cell_population.rGetMesh().GetNode(i)->rGetAppliedForce());
            }

            cpu_time = 1000.0*(c_end-c_start)/CLOCKS_PER_SEC;
            wall_time = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
        }

        {
            HoneycombVertexMeshGenerator generator(n_x_cells, n_y_cells);
            MutableVertexMesh<2,2>* vertex_mesh = generator.GetMesh();
            // Get a cell population
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            cells_generator.GenerateBasic(cells, vertex_mesh->GetNumElements(), std::vector<unsigned>());
            VertexBasedCellPopulation<2> cell_population(*vertex_mesh, cells);
            n_population_nodes = vertex_mesh->GetNumNodes();
            // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                cell_iter->SetBirthTime(-5.0);
            }

            MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
            p_growth_modifier->UpdateTargetAreas(cell_population);

            // Now let's make a FarhadifarForce and apply it to the population
            FarhadifarForce<2> force;
            std::clock_t c_start = std::clock();
            auto t_start = std::chrono::high_resolution_clock::now();
            force.AddForceContribution(cell_population);
            std::clock_t c_end = std::clock();
            auto t_end = std::chrono::high_resolution_clock::now();
            for (unsigned int i=0; i<n_population_nodes; ++i)
            {
                applied_forces_OMP.push_back(cell_population.rGetMesh().GetNode(i)->rGetAppliedForce());
            }

            cpu_time_OMP = 1000.0*(c_end-c_start)/CLOCKS_PER_SEC;
            wall_time_OMP = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
        }
        std::cout<<"CPU time, ms"<<std::endl;
        std::cout<<"In serial: "<<cpu_time<<" ; Using OMP: "<<cpu_time_OMP<<std::endl;
        std::cout<<"Wall time, ms"<<std::endl;
        std::cout<<"In serial: "<<wall_time<<" ; Using OMP: "<<wall_time_OMP<<std::endl;
        for (unsigned int i=0; i<applied_forces_OMP.size(); ++i)
        {
            TS_ASSERT_DELTA(applied_forces[i][0], applied_forces_OMP[i][0], 1e-12);
            TS_ASSERT_DELTA(applied_forces[i][1], applied_forces_OMP[i][1], 1e-12);
        }
    }

    void TestFarhadifarForceOMPMethods()
    {
        #ifndef _OPENMP
        TS_FAIL("OpenMP not used");
        return;
        #endif

        // This is the same test as for other vertex based forces. It comprises a sanity check that forces point in the right direction.
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        std::vector<double> angles = std::vector<double>(num_nodes);

        for (unsigned i=0; i<num_nodes; i++)
        {
            angles[i] = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, true, cos(angles[i]), sin(angles[i])));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Set up the cell
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(-1.0);
        cells.push_back(p_cell);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.InitialiseCells();

        // Create a force system
        FarhadifarForce<2> force;

        // Test get/set methods
        TS_ASSERT_DELTA(force.GetAreaElasticityParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(force.GetPerimeterContractilityParameter(), 0.04, 1e-12);
        TS_ASSERT_DELTA(force.GetLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetBoundaryLineTensionParameter(), 0.12, 1e-12);

        force.SetAreaElasticityParameter(5.8);
        force.SetPerimeterContractilityParameter(17.9);
        force.SetLineTensionParameter(0.5);
        force.SetBoundaryLineTensionParameter(0.6);

        TS_ASSERT_DELTA(force.GetAreaElasticityParameter(), 5.8, 1e-12);
        TS_ASSERT_DELTA(force.GetPerimeterContractilityParameter(), 17.9, 1e-12);
        TS_ASSERT_DELTA(force.GetLineTensionParameter(), 0.5, 1e-12);
        TS_ASSERT_DELTA(force.GetBoundaryLineTensionParameter(), 0.6, 1e-12);

        force.SetAreaElasticityParameter(1.0);
        force.SetPerimeterContractilityParameter(0.04);
        force.SetLineTensionParameter(0.12);
        force.SetBoundaryLineTensionParameter(0.12);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Currently, the Farhadifar force only works if used together with a target area growth modifier
        // This tests that a meaningful error appears if we don't use a growth modifier
        TS_ASSERT_THROWS_THIS(force.AddForceContribution(cell_population),
                              "You need to add an AbstractTargetAreaModifier to the simulation in order to use a FarhadifarForce");

        // create our modifier, which sets the target areas for the cell population

        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // The force on each node should be radially inward, with the same magnitude for all nodes
        double force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.25, 2);

        // Set the cell to be necrotic
        cell_population.GetCellUsingLocationIndex(0)->StartApoptosis();

        // Reset force vector
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        //#2488 workaround
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // The force on each node should not yet be affected by setting the cell to be apoptotic
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Increment time
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_DELTA(cell_population.GetCellUsingLocationIndex(0)->GetTimeUntilDeath(), 0.125, 1e-6);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        //#2488 workaround
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // Now the forces should be affected
        double apoptotic_force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());
        TS_ASSERT_LESS_THAN(force_magnitude, apoptotic_force_magnitude);
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), apoptotic_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -apoptotic_force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -apoptotic_force_magnitude*sin(angles[i]), 1e-4);
        }
    }

    void TestFarhadifarForceOMPTerms()
    {
        #ifndef _OPENMP
        TS_FAIL("OpenMP not used");
        return;
        #endif

        /**
         * Here we test that the forces are applied correctly to individual nodes.
         * We apply the force to something like this:
         *  . ____ . ____ .
         *  |      |      |
         *  |      |      |
         *  . ____ . ____ .
         */
        std::vector<Node<2>*> nodes;
        // the boolean says wether the node is a boundary node or not
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

        // make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->SetBirthTime(-5.0);
        }

        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        // Now let's make a FarhadifarForce and apply it to the population
        FarhadifarForce<2> force;

        force.AddForceContribution(cell_population);

        c_vector<double, 2> applied_force_0;
        applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1;
        applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        // If this is a Farhadifar force, this will be the force at the vertices
        TS_ASSERT_DELTA(applied_force_0[0], 3.44, 1e-10);
        TS_ASSERT_DELTA(applied_force_0[1], 3.44, 1e-10);
        TS_ASSERT_DELTA(applied_force_1[0], 0.0, 1e-10);
        TS_ASSERT_DELTA(applied_force_1[1], 6.76, 1e-10);
    }

    void TestFarhadifarForceOMPInSimulation()
    {
        #ifndef _OPENMP
        TS_FAIL("OpenMP not used");
        return;
        #endif

        /**
         * This is the same test as above, just that now we don't check that the applied forces are calculated correctly,
         * but rather that in a simulation the displacement of vertices is as we expect.
         *
         * This is the mesh:
         *  . ____ . ____ .
         *  |      |      |
         *  |      |      |
         *  . ____ . ____ .
         */
        std::vector<Node<2>*> nodes;
        // the boolean says wether the node is a boundary node or not
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

        // make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->SetBirthTime(-5.0);
        }

        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        // Now let's make a FarhadifarForce and add it to the simulation.
        MAKE_PTR(FarhadifarForce<2>, p_force);

        // We need to reset the cell rearrangement threshold - vertex movements are kept below that threshold
        cell_population.rGetMesh().SetCellRearrangementThreshold(0.5);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFarhadifarForce");
        simulator.SetEndTime(0.01);
        simulator.SetDt(0.01);
        simulator.AddForce(p_force);

        simulator.Solve();

        c_vector<double, 2> applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        // New Location = Old Location + (Dt * applied force), since viscosity should be one
        c_vector<double, 2> expected_new_node_location_0;
        expected_new_node_location_0[0] = 0.0+0.01*3.44;
        expected_new_node_location_0[1] = 0.0+0.01*3.44;
        c_vector<double, 2> expected_new_node_location_1;
        expected_new_node_location_1[0] = 2.0 + 0.01*0.0;
        expected_new_node_location_1[1] = 0.0 + 0.01*6.76;

        // If this is a Farhadifar force, this will be the location of the first two vertices.
        TS_ASSERT_DELTA(expected_new_node_location_0[0], (cell_population.rGetMesh().GetNode(0)->rGetLocation())[0], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_0[1], (cell_population.rGetMesh().GetNode(0)->rGetLocation())[1], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_1[0], (cell_population.rGetMesh().GetNode(1)->rGetLocation())[0], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_1[1], (cell_population.rGetMesh().GetNode(1)->rGetLocation())[1], 1e-10);

    }

    void TestFarhadifarForceOMPArchiving()
    {
        #ifndef _OPENMP
        TS_FAIL("OpenMP not used");
        return;
        #endif

        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "FarhadifarForce.arch";

        {
            FarhadifarForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetAreaElasticityParameter(5.8);
            force.SetPerimeterContractilityParameter(17.9);
            force.SetLineTensionParameter(0.5);
            force.SetBoundaryLineTensionParameter(0.6);

            // Serialize via pointer to most abstract class possible
            AbstractForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractForce<2>* p_abstract_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_abstract_force;

            FarhadifarForce<2>* p_farhadifar_force = static_cast<FarhadifarForce<2>*>(p_abstract_force);

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(p_farhadifar_force->GetAreaElasticityParameter(), 5.8, 1e-12);
            TS_ASSERT_DELTA(p_farhadifar_force->GetPerimeterContractilityParameter(), 17.9, 1e-12);
            TS_ASSERT_DELTA(p_farhadifar_force->GetLineTensionParameter(), 0.5, 1e-12);
            TS_ASSERT_DELTA(p_farhadifar_force->GetBoundaryLineTensionParameter(), 0.6, 1e-12);

            // Tidy up
            delete p_abstract_force;
        }
    }

};

#endif /*TESTFORCESOMP_HPP_*/
