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

#ifndef TESTRUSHLARSEN_HPP_
#define TESTRUSHLARSEN_HPP_

#include <cxxtest/TestSuite.h>

#include "LuoRudy1991.hpp"
#include "LuoRudy1991Opt.hpp"
#include "AbstractRushLarsenCardiacCell.hpp" // Needed for chaste_libs=0 build

#include "ZeroStimulus.hpp"
#include "SimpleStimulus.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"

#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "HeartConfig.hpp"
#include "CellMLToSharedLibraryConverter.hpp"
#include "DynamicCellModelLoader.hpp"
#include "Warnings.hpp"

/**
 * This test is based on code contributed by Megan Lewis, University of Saskatchewan.
 */
class TestRushLarsen : public CxxTest::TestSuite
{
    AbstractCardiacCell* mpRushLarsenCell;
    AbstractCardiacCell* mpRushLarsenCellOpt;

    void GenerateCells() throw (Exception)
    {
        // Do the conversions preserving generated sources
        CellMLToSharedLibraryConverter converter(true);
        std::string dirname = "TestRushLarsen";
        std::string model = "LuoRudy1991";
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus());

        std::vector<std::string> args;
        args.push_back("--rush-larsen");

        { // No opt
            // Copy CellML file into output dir
            OutputFileHandler handler(dirname + "/normal");
            FileFinder cellml_file("heart/src/odes/cellml/" + model + ".cellml", RelativeTo::ChasteSourceRoot);
            FileFinder copied_file = handler.CopyFileTo(cellml_file);

            // Create options file & convert
            converter.CreateOptionsFile(handler, model, args);
            DynamicCellModelLoader* p_loader = converter.Convert(copied_file);
            mpRushLarsenCell = dynamic_cast<AbstractCardiacCell*>(p_loader->CreateCell(p_solver, p_stimulus));
        }

        { // Optimised
            // Copy CellML file into output dir
            OutputFileHandler handler(dirname + "/opt");
            FileFinder cellml_file("heart/src/odes/cellml/" + model + ".cellml", RelativeTo::ChasteSourceRoot);
            FileFinder copied_file = handler.CopyFileTo(cellml_file);

            // Create options file & convert
            args.push_back("--opt");
            FileFinder conf_file("heart/src/odes/cellml/" + model + "-conf.xml", RelativeTo::ChasteSourceRoot);
            args.push_back("--conf=" + conf_file.GetAbsolutePath());
            converter.CreateOptionsFile(handler, model, args);
            DynamicCellModelLoader* p_loader = converter.Convert(copied_file);
            mpRushLarsenCellOpt = dynamic_cast<AbstractCardiacCell*>(p_loader->CreateCell(p_solver, p_stimulus));
        }
    }

public:
    void TestLuoRudyRushLarsenMethod()
    {
        HeartConfig::Instance()->SetOdeTimeStep(0.01);
        GenerateCells();

        // Check the models really use Rush-Larsen
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);

        // Normal Luo-Rudy for comparison
        boost::shared_ptr<EulerIvpOdeSolver> p_euler_solver(new EulerIvpOdeSolver());
        CellLuoRudy1991FromCellML reference_model(p_euler_solver, mpRushLarsenCell->GetStimulusFunction());
        CellLuoRudy1991FromCellMLOpt reference_model_opt(p_euler_solver, mpRushLarsenCell->GetStimulusFunction());

        // Check GetIIonic is identical
        TS_ASSERT_DELTA(mpRushLarsenCell->GetIIonic(), reference_model.GetIIonic(), 1e-12);
        TS_ASSERT_DELTA(mpRushLarsenCellOpt->GetIIonic(), reference_model_opt.GetIIonic(), 1e-12);


        // Test non-stimulated cell (using ComputeExceptVoltage)
        mpRushLarsenCell->ComputeExceptVoltage(0.0, 1.0);
        reference_model.ComputeExceptVoltage(0.0, 1.0);

        TS_ASSERT_EQUALS(mpRushLarsenCell->GetNumberOfStateVariables(),
                         reference_model.GetNumberOfStateVariables());
        for (unsigned i=0; i<reference_model.GetNumberOfStateVariables(); i++)
        {
            TS_ASSERT_DELTA(mpRushLarsenCell->rGetStateVariables()[i],
                            reference_model.rGetStateVariables()[i], 1e-6);
        }


        // Test stimulated cell (using Compute)
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(-25.5, 1.99, 0.0));
        mpRushLarsenCell->ResetToInitialConditions();
        mpRushLarsenCell->SetStimulusFunction(p_stimulus);
        OdeSolution solutions_RL = mpRushLarsenCell->Compute(0.0, 1.0, 0.01);
        TS_ASSERT_EQUALS(solutions_RL.GetNumberOfTimeSteps(), 100u);

        reference_model.ResetToInitialConditions();
        reference_model.SetStimulusFunction(p_stimulus);
        OdeSolution solutions_ref = reference_model.Compute(0.0, 1.0, 0.01);

        for (unsigned i=0; i<reference_model.GetNumberOfStateVariables(); i++)
        {
            TS_ASSERT_DELTA(solutions_RL.rGetSolutions().back()[i],
                            solutions_ref.rGetSolutions().back()[i], 1e-2);
        }


        // Test optimised cell (using SolveAndUpdateState)
        mpRushLarsenCellOpt->ResetToInitialConditions();
        mpRushLarsenCellOpt->SetStimulusFunction(p_stimulus);
        mpRushLarsenCellOpt->SolveAndUpdateState(0.0, 1.0);

        reference_model_opt.ResetToInitialConditions();
        reference_model_opt.SetStimulusFunction(p_stimulus);
        reference_model_opt.SolveAndUpdateState(0.0, 1.0);

        for (unsigned i=0; i<reference_model.GetNumberOfStateVariables(); i++)
        {
            TS_ASSERT_DELTA(mpRushLarsenCellOpt->rGetStateVariables()[i],
                            reference_model_opt.rGetStateVariables()[i], 1e-2);
        }


        // Test order of convergence (first-order method)

        // Get two solutions with halved stepsize
        mpRushLarsenCell->ResetToInitialConditions();
        mpRushLarsenCell->SetTimestep(0.01);
        OdeSolution solutions_RL_stimulated_cell_order1 = mpRushLarsenCell->Compute(0.0, 1.0);

        mpRushLarsenCell->ResetToInitialConditions();
        mpRushLarsenCell->SetTimestep(0.005);
        OdeSolution solutions_RL_stimulated_cell_order2 = mpRushLarsenCell->Compute(0.0, 1.0);

        // Get accurate solution to be used as reference solution
        reference_model.ResetToInitialConditions();
        reference_model.SetTimestep(1e-6);
        OdeSolution ref_solution = reference_model.Compute(0.0, 1.0);

        for (unsigned i=0; i<reference_model.GetNumberOfStateVariables(); i++)
        {
            double error1 = solutions_RL_stimulated_cell_order1.rGetSolutions().back()[i] - ref_solution.rGetSolutions().back()[i];
            double error2 = solutions_RL_stimulated_cell_order2.rGetSolutions().back()[i] - ref_solution.rGetSolutions().back()[i];
            TS_ASSERT_DELTA(error1/error2, 2, 0.1);
        }

        // Free memory
        delete mpRushLarsenCell;
        delete mpRushLarsenCellOpt;
    }

};

#endif // TESTRUSHLARSEN_HPP_
