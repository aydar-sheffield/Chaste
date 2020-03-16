/*
 * HEVertexMeshWriter.hpp
 *
 *  Created on: 13 Mar 2020
 *      Author: aydar
 */

#ifndef HEVERTEXMESHWRITER_HPP_
#define HEVERTEXMESHWRITER_HPP_

// Forward declaration prevents circular include chain
template<unsigned SPACE_DIM>
class HEVertexMesh;

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkConvexPointSet.h>
#include <vtkQuad.h>
#include <vtkPolygon.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDataCompressor.h>
#endif //CHASTE_VTK

#include "HEVertexMesh.hpp"
#include "AbstractMeshWriter.hpp"


// Forward declaration prevents circular include chain
template<unsigned SPACE_DIM>
class HEVertexMesh;

template<unsigned SPACE_DIM>
struct MeshWriterIterators;

/**
 * A mesh writer class for vertex-based meshes with edge representation. Each edge is associated with
 * trapezoid such that the base is the edge and the sides are parallel to the line joining the
 * cell's centroid and the nodes.
 *     __             __
 *   /    \         /    \
 *  /      \       /      \
 * |        |     |\      /|
 * |        |     | |    | |
 * |        | ==> |/      \|
 *  \      /       \      /
 *   \ __ /         \ __ /
 */
template<unsigned SPACE_DIM>
class HEVertexMeshWriter : public AbstractMeshWriter<SPACE_DIM, SPACE_DIM>
{
private:

    /**
     * If writing from a mesh object, the mesh to write to disk.
     * Otherwise NULL.
     */
    HEVertexMesh<SPACE_DIM>* mpMesh;

#ifdef CHASTE_VTK
    //Requires  "sudo aptitude install libvtk5-dev" or similar
    ///\todo Merge into VtkMeshWriter (#1076)
    vtkUnstructuredGrid* mpVtkUnstructedMesh;
#endif //CHASTE_VTK


public:

    /**
     * Constructor.
     *
     * @param rDirectory reference to the output directory, relative to where Chaste output is stored
     * @param rBaseName reference to the base name for results files
     * @param clearOutputDir whether to clear the output directory prior to writing files
     */
    HEVertexMeshWriter(const std::string& rDirectory,
                       const std::string& rBaseName,
                       const bool clearOutputDir=true);

    /**
     * Write VTK file using a mesh.
     *
     * @param rMesh reference to the vertex-based mesh
     * @param stamp is an optional stamp (like a time-stamp) to put into the name of the file
     */
    void WriteVtkUsingMesh(HEVertexMesh<SPACE_DIM>& rMesh, std::string stamp="");

    /**
     * Populate mpVtkUnstructedMesh using a vertex-based mesh.
     * Called by WriteVtkUsingMesh().
     *
     * @param rMesh reference to the vertex-based mesh
     */
    void MakeVtkMesh(HEVertexMesh<SPACE_DIM>& rMesh);

    /**
     * Add data to a future VTK file.
     *
     * @param dataName a tag to go into the VTK file
     * @param dataPayload a pay-load of length (number of elements)
     */
    void AddCellData(std::string dataName, std::vector<double> dataPayload);

    /**
     * Destructor.
     */
    ~HEVertexMeshWriter();

    void WriteFiles() override;
};

#endif /* MHEVERTEXMESHWRITER_HPP_ */
