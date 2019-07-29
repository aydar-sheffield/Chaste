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

#include "TrapEdgeVertexMeshWriter.hpp"


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::TrapEdgeVertexMeshWriter(const std::string &rDirectory,
                                                                           const std::string &rBaseName,
                                                                           const bool clearOutputDir)
        : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir),
        mpMesh(nullptr)
{

#ifdef CHASTE_VTK
    // Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
#endif //CHASTE_VTK

}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::~TrapEdgeVertexMeshWriter()
{

}


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteVtkUsingMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                                                                         std::string stamp)
{
#ifdef CHASTE_VTK
    assert(SPACE_DIM==3 || SPACE_DIM == 2);    // LCOV_EXCL_LINE

    // Create VTK mesh
    MakeVtkMesh(rMesh);

    // Now write VTK mesh to file
    assert(mpVtkUnstructedMesh->CheckAttributes() == 0);
    vtkXMLUnstructuredGridWriter* p_writer = vtkXMLUnstructuredGridWriter::New();
#if VTK_MAJOR_VERSION >= 6
    p_writer->SetInputData(mpVtkUnstructedMesh);
#else
    p_writer->SetInput(mpVtkUnstructedMesh);
#endif
    // Uninitialised stuff arises (see #1079), but you can remove valgrind problems by removing compression:
    // **** REMOVE WITH CAUTION *****
    p_writer->SetCompressor(nullptr);
    // **** REMOVE WITH CAUTION *****

    std::string vtk_file_name = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName;
    if (stamp != "")
    {
        vtk_file_name += "_" + stamp;
    }
    vtk_file_name += ".vtu";

    p_writer->SetFileName(vtk_file_name.c_str());
    //p_writer->PrintSelf(std::cout, vtkIndent());
    p_writer->Write();
    p_writer->Delete(); // Reference counted
#endif //CHASTE_VTK

}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::MakeVtkMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM> &rMesh)
{
#ifdef CHASTE_VTK
    // Make the Vtk mesh
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");
    /*
     * Populating points. First, outer points of elements
     */
    const unsigned n_vertices = rMesh.GetNumNodes();
    for (unsigned i=0; i<n_vertices; ++i)
    {
        c_vector<double, SPACE_DIM> position = rMesh.GetNode(i)->rGetLocation();
        if (SPACE_DIM==2)
        {
            p_pts->InsertPoint(i, position[0], position[1], 0.0);
        }
        else
        {
            p_pts->InsertPoint(i, position[0], position[1], position[2]);
        }
    }
    /*
     * Populating inner points.
     * [_________________][_____][_________]....[_____]
     *      ^^^^^^^^^^^    ^^^^^  ^^^^^^^^        ^^^^
     *  Outer points       Cell_1  Cell_2        Cell_{n_elements}
     *                              Inner Points
     * cell_offset_dist stores the distance from the beginning of p_pts array for each element
     * Note that the number of inner points equals to the number of nodes of each element
     */
    const unsigned n_elements = rMesh.GetNumElements();
    std::vector<unsigned> cell_offset_dist(rMesh.GetNumElements());
    cell_offset_dist[0] = n_vertices;
    for (unsigned i=1; i<n_elements; ++i)
        cell_offset_dist[i] = cell_offset_dist[i-1]+rMesh.GetElement(i-1)->GetNumNodes();
    //Coefficient 0<=alpha<=1 represents how thin the trapezoid is
    const double alpha = 0.8;
    typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem_end =rMesh.GetElementIteratorEnd();
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem = rMesh.GetElementIteratorBegin();
            elem != elem_end; ++elem)
    {
        const unsigned n_elem_nodes = elem->GetNumNodes();
        const c_vector<double, SPACE_DIM> elem_centroid = rMesh.GetCentroidOfElement(elem->GetIndex());
        for (unsigned elem_node_num = 0; elem_node_num < n_elem_nodes; elem_node_num++)
        {
            const c_vector<double, SPACE_DIM> node_position = elem->GetNode(elem_node_num)->rGetLocation();
            const double new_x = (node_position[0]-elem_centroid[0])*alpha + elem_centroid[0];
            const double new_y = (node_position[1]-elem_centroid[1])*alpha + elem_centroid[1];
            if (SPACE_DIM == 2)
            {
                p_pts->InsertPoint(cell_offset_dist[elem->GetIndex()]+elem_node_num,
                                   new_x, new_y, 0.0);
            }
            else
            {
                const double new_z = (node_position[2]-elem_centroid[2])*alpha + elem_centroid[2];
                p_pts->InsertPoint(cell_offset_dist[elem->GetIndex()]+elem_node_num,
                                   new_x, new_y, new_z);
            }
        }
    }
    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); // Reference counted
    /*Only 2D is fully supported.
     *For 3D, if "edge" should be synonymous with face
     *then, in principle, below should work
     */
    unsigned total_n_edges = 0;
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem = rMesh.GetElementIteratorBegin();
            elem != elem_end; ++elem)
    {
        if (SPACE_DIM==2)
        {
            //First do the trapezoids for each edge
            for (unsigned edge_index = 0; edge_index <elem->GetNumEdges(); ++edge_index)
            {
                vtkCell* p_cell;
                p_cell = vtkQuad::New();
                const unsigned n_trap_nodes = p_cell->GetNumberOfEdges(); //4 in 2D, 8 in 3D
                assert(n_trap_nodes == 4);
                vtkIdList* p_cell_id_list = p_cell->GetPointIds();
                p_cell_id_list->SetNumberOfIds(n_trap_nodes);
                auto p_edge = elem->GetEdge(edge_index);
                assert(p_edge->GetNumNodes()==2);
                //See the diagram above for storing pattern
                std::array<unsigned, 2> base_ids{p_edge->GetNode(0)->GetIndex(),
                    p_edge->GetNode(1)->GetIndex()};
                std::array<unsigned, 2> top_ids{elem->GetNodeLocalIndex(base_ids[0])
                    +cell_offset_dist[elem->GetIndex()],
                    elem->GetNodeLocalIndex(base_ids[1])+
                    cell_offset_dist[elem->GetIndex()]};
                //Assuming counter-clockwise ordering
                p_cell_id_list->SetId(0, base_ids[0]);
                p_cell_id_list->SetId(1, base_ids[1]);
                p_cell_id_list->SetId(2, top_ids[1]);
                p_cell_id_list->SetId(3, top_ids[0]);
                mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
                p_cell->Delete(); // Reference counted
                total_n_edges++;
            }
            //Now do the internal cell
            vtkCell* p_cell;
            p_cell = vtkPolygon::New();
            const unsigned n_elem_nodes = elem->GetNumNodes();
            vtkIdList* p_cell_id_list = p_cell->GetPointIds();
            p_cell_id_list->SetNumberOfIds(n_elem_nodes);
            for (unsigned j=0; j<n_elem_nodes; ++j)
            {
                p_cell_id_list->SetId(j,cell_offset_dist[elem->GetIndex()]+j);
            }

            mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);

            p_cell->Delete(); // Reference counted
        }
        else
        {
            /*
             * \todo For each face ...
             */
            vtkCell* p_cell;
            p_cell = vtkHexahedron::New();
            const unsigned n_trap_nodes = p_cell->GetNumberOfFaces();
            assert(n_trap_nodes == 6);
            p_cell->Delete();
        }
    }
    //For 2D case. For 3D, we should sum the total number of faces + n_elements
    assert(total_n_edges+n_elements==mpVtkUnstructedMesh->GetNumberOfCells());
#endif //CHASTE_VTK
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void
TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_scalars);
    p_scalars->Delete(); // Reference counted
#endif //CHASTE_VTK

}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    //Blank as we're only using the class for VTK at the moment
}

///////// Explicit instantiation///////

template class TrapEdgeVertexMeshWriter<1,1>;
template class TrapEdgeVertexMeshWriter<1,2>;
template class TrapEdgeVertexMeshWriter<1,3>;
template class TrapEdgeVertexMeshWriter<2,2>;
template class TrapEdgeVertexMeshWriter<2,3>;
template class TrapEdgeVertexMeshWriter<3,3>;
