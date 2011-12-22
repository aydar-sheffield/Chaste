/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and visualize vertex-based simulations on periodic meshes with different cell-cycle models =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize vertex-based simulations.
 * This mechanical model was originally proposed by T. Nagai and H. Honda ("A dynamic cell model for
 * the formation of epithelial tissues", Philosophical Magazine Part B 81:699-719).
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The remaining header files define classes that will be used in the cell population
 * simulation test. We have encountered some of these header files in previous cell-based
 * Chaste tutorials. */
#include "CellsGenerator.hpp"
#include "CryptCellsGenerator.hpp"
#include "WntConcentration.hpp"
#include "SloughingCellKiller.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/*
 * The next three header files define two different types of cell-cycle model,
 * one with fixed cell-cycle times and one
 * where the cell-cycle time depends on the Wnt concentration.
 */
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable mesh. */
#include "HoneycombVertexMeshGenerator.hpp"
/* The next header file defines a helper class for generating a periodic vertex mesh. */
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
/* The next header file defines the class that simulates the evolution of a crypt {{{CellPopulation}}}
 * for a vertex mesh. */
#include "CryptSimulation2d.hpp"
/* The next header file defines a vertex-based {{{CellPopulation}}} class.*/
#include "VertexBasedCellPopulation.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population, subject to each vertex.
 */
#include "NagaiHondaForce.hpp"

/* Next, we define the test class. */
class TestRunningVertexBasedCryptSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:

    /* EMPTYLINE
     *
     *
     * == Test 1 - create a vertex-based crypt simulation ==
     *
     * EMPTYLINE
     *
     * The first test generates a crypt, in which we use a cylindrical vertex mesh,
     * give each cell a fixed cell-cycle model, and enforce sloughing at the top of
     * the crypt.
     */
    void TestVertexBasedCrypt() throw(Exception)
    {
        /* Create a cylindrical mesh, and get the cell location indices. To enforce
         * periodicity at the left and right hand sides of the mesh, we use a subclass
         * called {{{Cylindrical2dMesh}}}, which has extra methods for maintaining
         * periodicity.
         */
          CylindricalHoneycombVertexMeshGenerator generator(6, 9);
         Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we the `CryptCellsGenerator` helper class, which is templated over the type
         * of cell model required (here {{{FixedDurationGenerationBasedCellCycleModel}}})
         * and the dimension. We create an empty vector of cells and pass this into the
         * method along with the mesh. The third argument 'true' indicates that the cells
         * should be assigned random birth times, to avoid synchronous division. The
         * {{{cells}}} vector is populated once the method {{{Generate}}} is
         * called.
         * The last four arguments represent the height below which cells belong to generations 0,
         * 1, 2, 3 and 4, respectively.
         */
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 1.0, 2.0, 3.0, 4.0);

        /* Create a cell population, as before. */
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        /* Create a simulator as before (except setting a different output directory). */
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("VertexCrypt");
        simulator.SetEndTime(0.1);

        /* Before running the simulation, we add a one or more force laws, which determine the mechanics of
         * the cell population.  For this test, we use a {{{NagaiHondaForce}}}.
         */
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* Before running the simulation, we add a cell killer. This object
         * dictates conditions under which cells die. For this test, we use
         * a {{{SloughingCellKiller}}}, which kills cells above a certain height.
         */
        double crypt_length = 6.0;
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }

    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/VertexCrypt/results_from_time_0}}}.
     * You may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
     * java executable.
     *
     * EMPTYLINE
     *
     * When we visualize the results, we should see three colours of cells: a row of blue stem cells, 3 rows of yellow transit
     * cells, and 5 rows of pink differentiated cells. Cells above 6.0 will be sloughed off immediately.
     *
     * EMPTYLINE
     *
     * == Test 2 - create a vertex-based crypt simulation with a simple wnt dependent cell-cycle model ==
     *
     * EMPTYLINE
     *
     * The next test generates a crypt, in which we use a cylindrical vertex mesh, and
     * impose a linearly decreasing concentration gradient of Wnt. Cells detect the level of Wnt
     * at their centre and those that are in a region of sufficient Wnt are defined to be transit cells,
     * whilst those above this Wnt threshold are defined to be differentiated. The cell cycle length of
     * transit cells is then assigned randomly from a uniform distribution.
     */
    void TestVertexBasedCryptWithSimpleWntCellCycleModel() throw(Exception)
    {
        /* Create a cylindrical mesh, and get the cell location indices, as before. */
        CylindricalHoneycombVertexMeshGenerator generator(6, 9);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        /* Create a {{{std::vector}}} of {{{CellPtr}}}s.
         * Generate cells, which are assigned a {{{SimpleWntCellCycleModel}}} using
         * the {{{CryptCellsGenerator}}}. The final boolean argument 'true' indicates
         * to assign randomly chosen birth times.
         */
        std::vector<CellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        /* Create a cell population, as before. */
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        /* Define the crypt length; this will be used for sloughing and calculating the Wnt gradient. */
        double crypt_length = 6.0;

        /* Set up a {{{WntConcentration}}} object, as in UserTutorials/RunningMeshBasedCryptSimulations. */
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        /* Create a simulator as before, and add a force law and sloughing cell killer to it. */
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("VertexCryptWithSimpleWntCellCycleModel");
        simulator.SetEndTime(0.1);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        /* Here we impose a boundary condition at the base: that cells
        * at the bottom of the crypt are repelled if they move past 0.*/
        simulator.UseJiggledBottomCells();

        /* Run the simulation, by calling {{{Solve()}}}. */
        simulator.Solve();
    }
    /*
    * EMPTYLINE
    *
    * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
    * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/VertexCryptWithSimpleWntCellCycleModel/results_from_time_0}}}.
    * You may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
    * java executable.
    *
    * EMPTYLINE
    *
    * When we visualize the results, we should see two colours of cells: yellow transit
    * cells and pink differentiated cells. Cells above 6.0 will be sloughed off immediately.
    */
};

#endif /* TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_ */
