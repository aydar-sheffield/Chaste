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
#ifndef TESTHEVERTEXMESHWRITER_HPP_
#define TESTHEVERTEXMESHWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "HEVertexMeshWriter.hpp"
#include "VtkMeshWriter.hpp"
#include "FileComparison.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 // Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

class TestVertexMeshWriter : public CxxTest::TestSuite
{
public:

    void TestHEVertexMeshWriterIn2d()
    {
        // Create a 2D mesh comprising seven nodes and two elements
        std::vector<HENode<2>*> nodes_2d;
        nodes_2d.push_back(new HENode<2>(0, false, 0.0, 0.0));
        nodes_2d.push_back(new HENode<2>(1, false, 1.0, 0.0));
        nodes_2d.push_back(new HENode<2>(2, false, 1.5, 1.0));
        nodes_2d.push_back(new HENode<2>(3, false, 1.0, 2.0));
        nodes_2d.push_back(new HENode<2>(4, false, 0.0, 1.0));
        nodes_2d.push_back(new HENode<2>(5, false, 2.0, 0.0));

        std::vector<std::vector<HENode<2>*> > nodes_elements_2d(2);
        for (unsigned i=0; i<5; i++)
        {
            nodes_elements_2d[0].push_back(nodes_2d[i]);
        }
        nodes_elements_2d[1].push_back(nodes_2d[1]);
        nodes_elements_2d[1].push_back(nodes_2d[5]);
        nodes_elements_2d[1].push_back(nodes_2d[2]);

        std::vector<HEElement<2>*> elements_2d;
        elements_2d.push_back(new HEElement<2>(0, nodes_elements_2d[0]));
        elements_2d.push_back(new HEElement<2>(1, nodes_elements_2d[1]));

        // Make a vertex mesh
        HEVertexMesh<2> basic_vertex_mesh(nodes_2d, elements_2d);

        // Create a vertex mesh writer
        HEVertexMeshWriter<2> vertex_mesh_writer("TestHEVertexMeshWriterIn2d", "vertex_mesh_2d");

        // TODO: adapt below to HEVertexMesh
        /*// Test files are written correctly
        vertex_mesh_writer.WriteFilesUsingMesh(basic_vertex_mesh);

        OutputFileHandler handler("TestVertexMeshWriterIn2d", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "vertex_mesh_2d.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "vertex_mesh_2d.cell";

        FileComparison comparer1(results_file1,"mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d.node");
        TS_ASSERT(comparer1.CompareFiles());

        FileComparison comparer2(results_file2,"mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d.cell");
        TS_ASSERT(comparer2.CompareFiles());*/

#ifdef CHASTE_VTK
        /**
         * Test below is general enough to check arbitrary meshes.
         */
        //The order of stored data is illustrated below:
        // [_____|_][____|_]
        //  ^^^^  ^
        // edge   cell interior
        const unsigned int num_elements = basic_vertex_mesh.GetNumElements();
        std::vector<unsigned int> cell_offset_dist(num_elements,0);
        unsigned int num_edges= 0;
        for (unsigned int i=1; i<num_elements; ++i)
        {
            cell_offset_dist[i] = cell_offset_dist[i-1] + basic_vertex_mesh.GetElement(i-1)->GetNumNodes()+1;
            num_edges+=basic_vertex_mesh.GetElement(i)->GetNumNodes();
        }
        num_edges += basic_vertex_mesh.GetElement(0)->GetNumNodes();
        const unsigned int total_num_values = basic_vertex_mesh.GetNumElements()+num_edges;

        //Cell data are all element ids
        std::vector<double> cell_ids(total_num_values);
        typename HEVertexMesh<2>::HEElementIterator elem_end =basic_vertex_mesh.GetElementIteratorEnd();
        for (typename HEVertexMesh<2>::HEElementIterator elem = basic_vertex_mesh.GetElementIteratorBegin();
                elem != elem_end; ++elem)
        {
            const unsigned int elem_index = elem->GetIndex();
            const unsigned int elem_num_edges = elem->GetNumNodes();
            for (unsigned int e=0; e<elem_num_edges; ++e)
                cell_ids[cell_offset_dist[elem_index]+e] = elem_index;
            cell_ids[cell_offset_dist[elem_index]+elem_num_edges] = elem_index;
        }

        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);

        vertex_mesh_writer.WriteVtkUsingMesh(basic_vertex_mesh);

        {
            OutputFileHandler handler("TestHEVertexMeshWriterIn2d", false);
            ///\todo #1076.  We need a way to test the contents of the VTK file
            std::string results_file3 = handler.GetOutputDirectoryFullPath() + "vertex_mesh_2d.vtu";
            FileFinder vtk_file(results_file3, RelativeTo::Absolute);
            TS_ASSERT(vtk_file.Exists());
        }
#else
        std::cout << "This test ran, but did not test VTK-dependent functions as VTK visualization is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestMeshWriterWithDeletedNode()
    {
        /*// Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 30u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 9u);


         * Delete element 0. This element contains 3 nodes that are
         * not contained in any other element and so will be marked
         * as deleted.

        mesh.DeleteElementPriorToReMesh(0);

        // Write mesh to file
        VertexMeshWriter<2,2> mesh_writer("TestMeshWriterWithDeletedNode", "vertex_mesh");
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        // Read mesh back in from file
        std::string output_dir = mesh_writer.GetOutputDirectory();
        VertexMeshReader<2,2> mesh_reader2(output_dir + "vertex_mesh");

        // We should have one less element and three less nodes
        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 27u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 8u);*/
    }

    void TestReadingAndWritingElementAttributes()
    {
        /*// Read in a mesh with element attributes
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshReader2d/vertex_mesh_with_element_attributes");
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        // Construct the mesh
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetUnsignedAttribute(), 97u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetUnsignedAttribute(), 152u);

        // Write the mesh to file
        // Nested scope so the reader is destroyed before we try writing to the folder again
        {
            VertexMeshWriter<2,2> mesh_writer("TestReadingAndWritingElementAttributes", "vertex_mesh_with_element_attributes");
            mesh_writer.WriteFilesUsingMesh(mesh);

            // Now read in the mesh that was written
            OutputFileHandler handler("TestReadingAndWritingElementAttributes", false);
            VertexMeshReader<2,2> mesh_reader2(handler.GetOutputDirectoryFullPath() + "vertex_mesh_with_element_attributes");
            TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 2u);
            TS_ASSERT_EQUALS(mesh_reader2.GetNumElementAttributes(), 1u);

            // Construct the mesh again
            VertexMesh<2,2> mesh2;
            mesh2.ConstructFromMeshReader(mesh_reader);
            TS_ASSERT_EQUALS(mesh2.GetElement(0)->GetUnsignedAttribute(), 97u);
            TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetUnsignedAttribute(), 152u);
        }

        // For coverage, repeat this test for a vertex mesh whose elements have faces
        VertexMeshReader<3,3> mesh_reader3d("mesh/test/data/TestVertexMeshWriter/vertex_mesh_3d_with_faces_and_attributes");
        TS_ASSERT_EQUALS(mesh_reader3d.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh_reader3d.GetNumElementAttributes(), 1u);

        // Construct the mesh
        VertexMesh<3,3> mesh3d;
        mesh3d.ConstructFromMeshReader(mesh_reader3d);
        TS_ASSERT_EQUALS(mesh3d.GetElement(0)->GetUnsignedAttribute(), 49u);

        // Write the mesh to file
        VertexMeshWriter<3,3> mesh_writer3d("TestReadingAndWritingElementAttributes", "vertex_mesh_3d_with_faces_and_attributes");
        mesh_writer3d.WriteFilesUsingMesh(mesh3d);

        // Now read in the mesh that was written
        OutputFileHandler handler3d("TestReadingAndWritingElementAttributes", false);
        VertexMeshReader<3,3> mesh_reader3d2(handler3d.GetOutputDirectoryFullPath() + "vertex_mesh_3d_with_faces_and_attributes");
        TS_ASSERT_EQUALS(mesh_reader3d2.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh_reader3d2.GetNumElementAttributes(), 1u);

        // Construct the mesh again
        VertexMesh<3,3> mesh3d2;
        mesh3d2.ConstructFromMeshReader(mesh_reader3d2);
        TS_ASSERT_EQUALS(mesh3d2.GetElement(0)->GetUnsignedAttribute(), 49u);*/
    }

    void TestWriteFilesUsingMeshReader()
    {
        /*// Create a VertexMeshReader and use it to write mesh files
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshReader2d/vertex_mesh_with_element_attributes");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        VertexMeshWriter<2,2> mesh_writer("TestWriteFilesUsingMeshReader", "vertex_mesh");
        mesh_writer.WriteFilesUsingMeshReader(mesh_reader);

        // Now read in the mesh that was written
        OutputFileHandler handler("TestWriteFilesUsingMeshReader", false);
        VertexMeshReader<2,2> mesh_reader2(handler.GetOutputDirectoryFullPath() + "vertex_mesh");
        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElementAttributes(), 1u);*/
    }
};

#endif /*TESTHEVERTEXMESHWRITER_HPP_*/
