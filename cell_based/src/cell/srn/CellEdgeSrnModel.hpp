//
// Created by twin on 22/01/19.
//

#ifndef CHASTE_CELLEGESRNMODEL_H
#define CHASTE_CELLEGESRNMODEL_H

#include <vector>
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "AbstractSrnModel.hpp"
#include "CellCycleModelOdeHandler.hpp"
#include "SimulationTime.hpp"


class CellEdgeSrnModel : public AbstractSrnModel {

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the SRN model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSrnModel>(*this);
        archive & boost::serialization::base_object<CellCycleModelOdeHandler>(*this);
//        archive & mInitialConditions;
//        archive & mStateSize;
    }

    std::vector<boost::shared_ptr<AbstractSrnModel>> mEdgeSrnModels;
    using abstractsrnmodel_t = std::vector<boost::shared_ptr<AbstractSrnModel>>;

protected:


    CellEdgeSrnModel(const CellEdgeSrnModel &rModel) : AbstractSrnModel(rModel) {

        //Makes a copy of all SRN models inside the system
        for(auto srnModel: rModel.mEdgeSrnModels){
            this->mEdgeSrnModels.push_back(boost::shared_ptr<AbstractSrnModel>(srnModel->CreateSrnModel()));
        }
    }

public:

    /* Makes the class iterable which returns the individual edge srn models */
    using iterator = abstractsrnmodel_t::iterator;
    using const_iterator = abstractsrnmodel_t::const_iterator;
    iterator begin() { return mEdgeSrnModels.begin(); }
    iterator end() { return mEdgeSrnModels.end(); }
    const_iterator begin() const { return mEdgeSrnModels.begin(); }
    const_iterator end() const { return mEdgeSrnModels.end(); }
    const_iterator cbegin() const { return mEdgeSrnModels.cbegin(); }
    const_iterator cend() const { return mEdgeSrnModels.cend(); }


    CellEdgeSrnModel(){}

    ~CellEdgeSrnModel(){}

    void Initialise() override
    {
        for(auto edgeModel: mEdgeSrnModels)
        {
            edgeModel->Initialise();
        }
    }

    void SimulateToCurrentTime() override
    {
        for(auto srnModel: mEdgeSrnModels){
            srnModel->SimulateToCurrentTime();
        }

    }

    AbstractSrnModel* CreateSrnModel() override
    {
        return new CellEdgeSrnModel(*this);
    }

    void AddEdgeSrn(std::vector<boost::shared_ptr<AbstractSrnModel>> edgeSrn)
    {
        mEdgeSrnModels = edgeSrn;
    }

    void AddEdgeSrn(boost::shared_ptr<AbstractSrnModel> edgeSrn)
    {
        edgeSrn->SetEdgeLocalIndex(mEdgeSrnModels.size());
        mEdgeSrnModels.push_back(edgeSrn);

    }

    void InsertEdgeSrn(unsigned index, boost::shared_ptr<AbstractSrnModel> edgeSrn)
    {
        mEdgeSrnModels.insert(mEdgeSrnModels.begin() + index, edgeSrn);
    }

    boost::shared_ptr<AbstractSrnModel> RemoveEdgeSrn(unsigned index)
    {
        auto edgeSrn = mEdgeSrnModels[index];
        mEdgeSrnModels.erase(mEdgeSrnModels.begin() + index);
        return edgeSrn;
    }

    unsigned GetNumEdgeSrn()
    {
        return mEdgeSrnModels.size();
    }

    boost::shared_ptr<AbstractSrnModel> GetEdgeSrn(unsigned index)
    {
        assert(index < mEdgeSrnModels.size());
        return mEdgeSrnModels[index];
    }

    void SetCell(CellPtr pCell) override {
        AbstractSrnModel::SetCell(pCell);
        //Makes a copy of all SRN models inside the system
        for(auto srnModel: mEdgeSrnModels){
            srnModel->SetCell(pCell);
        }
    }

};


#endif //CHASTE_CELLEGESRNMODEL_H
