/*
 * HEVertexMeshWriter.cpp
 *
 *  Created on: 13 Mar 2020
 *      Author: aydar
 */

#include "HEVertexMeshWriter.hpp"

template<unsigned int SPACE_DIM>
HEVertexMeshWriter<SPACE_DIM>::HEVertexMeshWriter(const std::string &rDirectory,
                                                  const std::string &rBaseName,
                                                  const bool clearOutputDir)
        : AbstractMeshWriter<SPACE_DIM,SPACE_DIM>(rDirectory, rBaseName, clearOutputDir),
        mpMesh(nullptr)
{

#ifdef CHASTE_VTK
    // Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
#endif //CHASTE_VTK

}

template<unsigned int SPACE_DIM>
HEVertexMeshWriter<SPACE_DIM>::~HEVertexMeshWriter()
{

}


template<unsigned int SPACE_DIM>
void HEVertexMeshWriter<SPACE_DIM>::WriteVtkUsingMesh(HEVertexMesh<SPACE_DIM> &rMesh,
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

template<unsigned int SPACE_DIM>
void HEVertexMeshWriter<SPACE_DIM>::MakeVtkMesh(HEVertexMesh<SPACE_DIM> &rMesh)
{

#ifdef CHASTE_VTK
    // Make the Vtk mesh
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");

    //Populating points. First, outer points of elements
    const unsigned n_vertices = rMesh.GetNumNodes();
    for (unsigned i=0; i<n_vertices; ++i)
    {
        c_vector<double, SPACE_DIM> position;
        position = rMesh.GetNode(i)->rGetLocation();
        if (SPACE_DIM==2)
        {
            p_pts->InsertPoint(i, position[0], position[1], 0.0);
        }
        else
        {
            p_pts->InsertPoint(i, position[0], position[1], position[2]);
        }
    }

     /* Populating inner points.
     * [_________________][_____][_________]....[_____]
     *      ^^^^^^^^^^^    ^^^^^  ^^^^^^^^        ^^^^
     *  Outer points       Cell_1  Cell_2        Cell_{n_elements}
     *                              Inner Points
     * cell_offset_dist stores the distance from the beginning of p_pts array for each element
     * Note that the number of inner points equals to the number of nodes of each element
     * */

    const unsigned n_elements = rMesh.GetNumElements();
    std::vector<unsigned> cell_offset_dist(rMesh.GetNumElements());
    cell_offset_dist[0] = n_vertices;
    for (unsigned i=1; i<n_elements; ++i)
        cell_offset_dist[i] = cell_offset_dist[i-1]+rMesh.GetElement(i-1)->GetNumNodes();
    //Coefficient 0<=alpha<=1 represents how thin the trapezoid is
    const double alpha = 0.8;
    typename HEVertexMesh<SPACE_DIM>::HEElementIterator elem_end =rMesh.GetElementIteratorEnd();
    for (typename HEVertexMesh<SPACE_DIM>::HEElementIterator elem = rMesh.GetElementIteratorBegin();
            elem != elem_end; ++elem)
    {
        const unsigned int n_elem_nodes = elem->GetNumNodes();
        c_vector<double, SPACE_DIM> elem_centroid;
        elem_centroid = rMesh.GetCentroidOfElement(elem->GetIndex());
        //Iterate over element nodes
        HalfEdge<SPACE_DIM>* edge = elem->GetHalfEdge();
        HalfEdge<SPACE_DIM>* next_edge = edge;
        unsigned int elem_node_num = 1;
        do
        {
            c_vector<double, SPACE_DIM> node_position;
            node_position = next_edge->GetTargetVertex()->rGetLocation();
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
            next_edge = next_edge->GetNextHalfEdge();
            //Node traversal starts at the FIRST local node
            elem_node_num = (elem_node_num+1)%n_elem_nodes;
        }while(next_edge != edge);

    }
    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); // Reference counted
    /* Only 2D is fully supported.
     * For 3D, if "edge" should be synonymous with face
     * then, in principle, below should work
     * */
    unsigned total_n_edges = 0;
    for (typename HEVertexMesh<SPACE_DIM>::HEElementIterator elem = rMesh.GetElementIteratorBegin();
            elem != elem_end; ++elem)
    {
        if (SPACE_DIM==2)
        {
            HalfEdge<SPACE_DIM>* edge = elem->GetHalfEdge();
            HalfEdge<SPACE_DIM>* next_edge = edge;
            //First do the trapezoids for each edge
            unsigned int local_index = 0;
            unsigned int next_index = 0;
            const unsigned int num_edges = elem->GetNumNodes();
            do
            {
                vtkCell* p_cell;
                p_cell = vtkQuad::New();
                const unsigned n_trapez_nodes = p_cell->GetNumberOfEdges(); //4 in 2D, 8 in 3D
                assert(n_trapez_nodes == 4);

                vtkIdList* p_cell_id_list = p_cell->GetPointIds();
                p_cell_id_list->SetNumberOfIds(n_trapez_nodes);

                //See the diagram above for storing pattern
                std::array<unsigned, 2> base_ids{next_edge->GetTwinHalfEdge()->GetTargetVertex()->GetIndex(),
                    next_edge->GetTargetVertex()->GetIndex()};
                next_index = (local_index+1)%num_edges;
                std::array<unsigned, 2> top_ids{local_index+cell_offset_dist[elem->GetIndex()],
                    next_index+cell_offset_dist[elem->GetIndex()]};
                //Assuming counter-clockwise ordering
                p_cell_id_list->SetId(0, base_ids[0]);
                p_cell_id_list->SetId(1, base_ids[1]);
                p_cell_id_list->SetId(2, top_ids[1]);
                p_cell_id_list->SetId(3, top_ids[0]);
                mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
                p_cell->Delete(); // Reference counted

                total_n_edges++;
                next_edge = next_edge->GetNextHalfEdge();
                local_index++;
            }while(next_edge!=edge);

            //Now do the internal cell
            vtkCell* p_cell;
            p_cell = vtkPolygon::New();
            vtkIdList* p_cell_id_list = p_cell->GetPointIds();
            p_cell_id_list->SetNumberOfIds(num_edges);
            for (unsigned j=0; j<num_edges; ++j)
            {
                p_cell_id_list->SetId(j,cell_offset_dist[elem->GetIndex()]+j);
            }

            mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);

            p_cell->Delete(); // Reference counted
        }
        else
        {

             /* \todo For each face ...*/

            vtkCell* p_cell;
            p_cell = vtkHexahedron::New();
            assert(p_cell->GetNumberOfFaces() == 6);
            p_cell->Delete();
        }
    }
    //For 2D case. For 3D, we should sum the total number of faces + n_elements
    assert(total_n_edges+n_elements==mpVtkUnstructedMesh->GetNumberOfCells());
#endif //CHASTE_VTK

}

template<unsigned int SPACE_DIM>
void HEVertexMeshWriter<SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
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

template<unsigned SPACE_DIM>
void HEVertexMeshWriter<SPACE_DIM>::AddPointData(std::string dataName, std::vector<double> dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkPointData* p_point_data = mpVtkUnstructedMesh->GetPointData();
    p_point_data->AddArray(p_scalars);
    p_scalars->Delete(); // Reference counted
#endif //CHASTE_VTK
}

template<unsigned int SPACE_DIM>
void HEVertexMeshWriter<SPACE_DIM>::WriteFiles()
{
    //Blank as we're only using the class for VTK at the moment
}

///////// Explicit instantiation///////

template class HEVertexMeshWriter<1>;
template class HEVertexMeshWriter<2>;
template class HEVertexMeshWriter<3>;

