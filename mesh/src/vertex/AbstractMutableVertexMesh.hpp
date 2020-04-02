/*
 * AbstractMutableVertexMesh.hpp
 *
 *  Created on: 30 Mar 2020
 *      Author: aydar
 */

#ifndef ABSTRACTMUTABLEVERTEXMESH_HPP_
#define ABSTRACTMUTABLEVERTEXMESH_HPP_
#include "VertexElementMap.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <vector>

/**
 * This is a base class for mutable meshes. Derived classes can utilise distinct data structures
 * for mesh representation. Hence, remeshing and associated methods can be implemented in derived classses.
 */
template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
class AbstractMutableVertexMesh
{
protected:
    /** Indices of nodes that have been deleted. These indices can be reused when adding new elements/nodes. */
    std::vector<unsigned> mDeletedNodeIndices;

    /** Indices of elements that have been deleted. These indices can be reused when adding new elements. */
    std::vector<unsigned> mDeletedElementIndices;

    /** The minimum distance apart that two nodes in the mesh can be without causing element rearrangement. */
    double mCellRearrangementThreshold;

    /**
     * The ratio between the minimum distance apart that two nodes in the mesh can be without causing element
     * rearrangement and their separation after remeshing.
     */
    double mCellRearrangementRatio;

    /** The area threshold at which T2 swaps occur in an apoptotic, triangular cell/element. */
    double mT2Threshold;

    /** The probability that, instead of a T1 swap, the relevant nodes merge to form a protorosette */
    double mProtorosetteFormationProbability;

    /** The probability that, in a given timestep, a protorosette node resolves into two rank-3 nodes */
    double mProtorosetteResolutionProbabilityPerTimestep;

    /** The probability that, in a given timestep, a rosette node resolves into two lower-rank nodes */
    double mRosetteResolutionProbabilityPerTimestep;

    /** Whether to check for edges intersections (true) or not (false). */
    bool mCheckForInternalIntersections;

    /**
     * Distance for T3 swap checking. At each time step we check for each boundary node whether
     * it intersects with any boundary elements (cells) whose centroids lie within this distance
     * to the node. Note that T3 swaps may not be resolved correctly if this distance is chosen
     * too small, while large values for mDistanceForT3SwapChecking may slow down the simulation.
     */
    double mDistanceForT3SwapChecking;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialize the mesh.
     *
     * Note that if you are calling this method (from subclasses) you should archive your
     * member variables FIRST. So that this method can call a ReMesh
     * (to convert from TrianglesMeshReader input format into your native format).
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // NOTE - Subclasses must archive their member variables BEFORE calling this method.
        archive & mDeletedNodeIndices;
        archive & mDeletedElementIndices;
        archive & mCellRearrangementThreshold;
        archive & mCellRearrangementRatio;
        archive & mT2Threshold;
        archive & mProtorosetteFormationProbability;
        archive & mProtorosetteResolutionProbabilityPerTimestep;
        archive & mRosetteResolutionProbabilityPerTimestep;
        archive & mCheckForInternalIntersections;
        archive & mDistanceForT3SwapChecking;
    }

public:

    /**
     * Default constructor.
     *
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     * @param cellRearrangementRatio ratio between the minimum threshold distance for element
     *                                rearrangement node separation after remeshing (defaults to 1.5)
     * @param protorosetteFormationProbability the probability of a protorosette formation event happening instead of
     *                                a T1 swap
     * @param protorosetteResolutionProbabilityPerTimestep the probability that, in a given timestep, a protorosette
     *                                will resolve (similar to the completion of a T1 swap)
     * @param rosetteResolutionProbabilityPerTimestep the probability that, in a given timestep, a rosette will
     *                                resolve (reduce the number of cells sharing a common vertex by 1)
     */
    AbstractMutableVertexMesh(double cellRearrangementThreshold=0.01,
                              double t2Threshold=0.001,
                              double cellRearrangementRatio=1.5,
                              double protorosetteFormationProbability=0.0,
                              double protorosetteResolutionProbabilityPerTimestep=0.0,
                              double rosetteResolutionProbabilityPerTimestep=0.0);

    virtual ~AbstractMutableVertexMesh();

    /**
     * Set method for mCellRearrangementThreshold.
     *
     * @param cellRearrangementThreshold
     */
    void SetCellRearrangementThreshold(double cellRearrangementThreshold);

    /**
     * Set method for mT2Threshold.
     *
     * @param t2Threshold
     */
    void SetT2Threshold(double t2Threshold);

    /**
     * Set method for mCellRearrangementRatio.
     *
     * @param cellRearrangementRatio
     */
    void SetCellRearrangementRatio(double cellRearrangementRatio);

    /**
     * Set method for mProtoRosetteFormationProbability.
     *
     * @param protorosetteFormationProbability the new value of mProtoRosetteFormationProbability
     */
    void SetProtorosetteFormationProbability(double protorosetteFormationProbability);

    /**
     * Set method for mProtoRosetteResolutionProbabilityPerTimestep.
     *
     * @param protorosetteResolutionProbabilityPerTimestep the new value of mProtoRosetteResolutionProbabilityPerTimestep
     */
    void SetProtorosetteResolutionProbabilityPerTimestep(double protorosetteResolutionProbabilityPerTimestep);

    /**
     * Set method for mRosetteResolutionProbabilityPerTimestep.
     *
     * @param rosetteResolutionProbabilityPerTimestep the new value of mRosetteResolutionProbabilityPerTimestep
     */
    void SetRosetteResolutionProbabilityPerTimestep(double rosetteResolutionProbabilityPerTimestep);

    /**
     * Set method for mCheckForInternalIntersections.
     *
     * @param checkForInternalIntersections
     */
    void SetCheckForInternalIntersections(bool checkForInternalIntersections);

    /**
     * @return mCellRearrangementThreshold
     */
    double GetCellRearrangementThreshold() const;

    /**
     * @return mT2Threshold
     */
    double GetT2Threshold() const;

    /**
     * @return mCellRearrangementRatio
     */
    double GetCellRearrangementRatio() const;

    /**
     * Get method for mProtoRosetteFormationProbability.
     *
     * @return mProtoRosetteFormationProbability
     */
    double GetProtorosetteFormationProbability() const;

    /**
     * Get method for mProtoRosetteResolutionProbabilityPerTimestep.
     *
     * @return mProtoRosetteResolutionProbabilityPerTimestep
     */
    double GetProtorosetteResolutionProbabilityPerTimestep() const;

    /**
     * Get method for mRosetteResolutionProbabilityPerTimestep.
     *
     * @return mRosetteResolutionProbabilityPerTimestep
     */
    double GetRosetteResolutionProbabilityPerTimestep() const;

    /**
     * Set distance for T3 swap checking. At each time step we check for each boundary node whether
     * it intersects with any boundary elements (cells) whose centroids lie within this distance
     * to the node. Note that T3 swaps may not be resolved correctly if this distance is chosen
     * too small, while large values for mDistanceForT3SwapChecking may slow down the simulation.
     *
     * @param distanceForT3SwapChecking
     */
    void SetDistanceForT3SwapChecking( double distanceForT3SwapChecking );

    /**
     * Get Distance for T3 swap checking.
     *
     * @return mDistanceForT3SwapChecking
     */
    double GetDistanceForT3SwapChecking() const;

    /**
     * @return mCheckForInternalIntersections, either to check for edges intersections or not.
     */
    bool GetCheckForInternalIntersections() const;

    /**
     * Update the state of the mesh by implementing any local remeshing operations (node merging,
     * or T1, T2 or T3 swaps) that are required, and store any changes in element indices using
     * the given VertexElementMap.
     *
     * This method calls several other methods, in particular CheckForT2Swaps(), CheckForSwapsFromShortEdges()
     * and CheckForIntersections().
     *
     * @param rElementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     */
    virtual void ReMesh(VertexElementMap& rElementMap);

};
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AbstractMutableVertexMesh)

#endif /* ABSTRACTMUTABLEVERTEXMESH_HPP_ */
