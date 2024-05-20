//! @file LmrRate.h
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LMRRATE_H
#define CT_LMRRATE_H
#include "cantera/kinetics/Arrhenius.h"
#include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/ChebyshevRate.h"
#include "cantera/kinetics/PlogRate.h"
// #include "cantera/base/Array.h"
// #include "cantera/base/FactoryBase.h"
// #include "cantera/kinetics/Reaction.h"

namespace Cantera{

struct LmrData : public ReactionData{
    // LmrData() = default;
    LmrData();

    void update(double T, double P) override {
        ReactionData::update(T);
        pressure = P;
        logP = std::log(P);
    }

    bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;
    void perturbPressure(double deltaP);
    void restore() override;

    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        moleFractions.resize(nSpecies, NAN);
        ready = true;
    }

    void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }
    double pressure = NAN; //!< pressure
    double logP = 0.0; //!< logarithm of pressure
    bool ready = false; //!< boolean indicating whether vectors are accessible
    vector<double> moleFractions;
    int mfNumber; 
    vector<string> allSpecies_; //list of all yaml species (not just those for which LMRR data exists)  
    
// protected:
    double m_pressure_buf = -1.0; //!< buffered pressure
};



class LmrRate final : public ReactionRate
{
public:
    LmrRate() = default;//! Default constructor.
    explicit LmrRate(const std::multimap<double, ArrheniusRate>& rates);
    LmrRate(const AnyMap& node, const UnitStack& rate_units={});
    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<LmrRate, LmrData>>();
    }
    const string type() const override { return "LMR_R"; } //! Identifier of reaction rate type

    // map<string, pair<const AnyMap&,const UnitStack&>> colliderInfo;
    map<string, pair<AnyMap,UnitStack>> colliderInfo;
    
    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;
    void getParameters(AnyMap& rateNode, const Units& rate_units) const;
    void getParameters(AnyMap& rateNode) const override {
        return getParameters(rateNode, Units(0));
    }

    // void syncPlogData(const LmrData& shared_data, PlogData plog_data);
    // void syncTroeData(const LmrData& shared_data, FalloffData troe_data);
    // void syncChebData(const LmrData& shared_data, ChebyshevData cheb_data);

    // void syncPlogData(const LmrData& shared_data, PlogData);
    double evalPlogRate(const LmrData& shared_data, double eig0val);
    double evalTroeRate(const LmrData& shared_data, double eig0val);
    double evalChebyshevRate(const LmrData& shared_data, double eig0val);

    double evalFromStruct(const LmrData& shared_data);
    void validate(const string& equation, const Kinetics& kin) override; //removed from cpp, but re-insert later

    double eig0_mix=0.0;
    double eig0_M;
    AnyMap colliders_i;
    UnitStack rate_units_i;




protected:
    double logP_ = -1000;
};


}
#endif