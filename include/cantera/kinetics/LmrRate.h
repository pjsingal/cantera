//! @file LmrRate.h
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LMRRATE_H
#define CT_LMRRATE_H
#include "cantera/kinetics/Arrhenius.h"
#include "cantera/base/Array.h"
#include "cantera/base/FactoryBase.h"
// #include "cantera/kinetics/Reaction.h"

namespace Cantera{

// struct LmrTroeData : public ReactionData
// {
//     //Need to actively assign values
//     double T;
//     double logT;
//     double recipT;
//     bool ready;
//     vector<double> conc_3b; //!< vector of effective third-body concentrations
//     AnyMap node;
//     UnitStack rate_units;
// };

// struct LmrTroeRate : public ReactionRate
// {
//     //Need to actively assign values
//     vector<double> m_work; //!< Work vector
//     double m_rc_low = NAN; //!< Evaluated reaction rate in the low-pressure limit
//     double m_rc_high = NAN; //!< Evaluated reaction rate in the high-pressure limit
//     size_t m_rate_index; //assign it as m_rate_index=npos;

//     //Don't need to assign values
//     bool m_chemicallyActivated = false;
    
//     bool m_negativeA_ok;
//     ArrheniusRate m_lowRate; //!< The reaction rate in the low-pressure limit
//     ArrheniusRate m_highRate; //!< The reaction rate in the high-pressure limit

//     //Functions needed for variable definitions
//     //! Get reaction rate in the low-pressure limit
//     ArrheniusRate& lowRate() {
//         return m_lowRate;
//     }
//     //! Get reaction rate in the high-pressure limit
//     ArrheniusRate& highRate() {
//         return m_highRate;
//     }
// };

// struct LmrTroe : public
// {
//     LmrTroeData data;

// };

struct LmrData : public ReactionData{
    // LmrData() = default;
    LmrData();
    
    //recent adds
    void update(double T, double P) override {
        ReactionData::update(T);
        pressure = P;
        logP = std::log(P);
    }
    //end of recent adds

    bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;
    void perturbPressure(double deltaP);
    void restore() override;

    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        conc_3b.resize(nReactions, NAN); //Troe param
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
    vector<double> conc_3b; //Troe param
    
protected:
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

    // map<string, AnyMap> plogList_;
    // map<string, AnyMap> troeList_;
    // map<string, AnyMap> chebList_;

    map<string, PlogRate*> plogObjects;
    map<string, TroeRate*> troeObjects;
    map<string, ChebyshevRate*> chebObjects;

    // new TroeRate(colliders[i],rate_units); //Borrowed this syntax from reg() in rxn rate factory 
    
    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;
    void getParameters(AnyMap& rateNode, const Units& rate_units) const;
    void getParameters(AnyMap& rateNode) const override {
        return getParameters(rateNode, Units(0));
    }
    double evalFromStruct(const LmrData& shared_data);
    void validate(const string& equation, const Kinetics& kin) override; //removed from cpp, but re-insert later
    UnitStack rate_units_;
    map<string,ArrheniusRate> eig0_;
    // map<string, map<double, pair<size_t, size_t>>> pressures_;
    // map<string, vector<ArrheniusRate>> rates_;
    // double logPeff_;
    // double eig0_mix_ = 0.0;
    // double log_eig0_mix_ = 0.0;
    // double k_LMR_;

    // //Troe member variables
    // map<string,ArrheniusRate> k0_;
    // map<string,ArrheniusRate> kinf_;
    // map<string,double> m_a;//! parameter a in the 4-parameter Troe falloff function. Dimensionless
    // map<string,double> m_rt3;//! parameter 1/T_3 in the 4-parameter Troe falloff function. [K^-1]
    // map<string,double> m_rt1;//! parameter 1/T_1 in the 4-parameter Troe falloff function. [K^-1]
    // map<string,double> m_t2;//! parameter T_2 in the 4-parameter Troe falloff function. [K]

    // //Chebyshev member variables
    // map<string,double> Tmin_, Tmax_; //!< valid temperature range [K]
    // map<string,double> Pmin_, Pmax_; //!< valid pressure range [Pa]
    // map<string,double> TrNum_, TrDen_; //!< terms appearing in the reduced temperature
    // map<string,double> PrNum_, PrDen_; //!< terms appearing in the reduced pressure
    // map<string,Array2D> m_coeffs_; //!<< coefficient array
    // map<string,vector<double>> dotProd_; //!< dot product of coeffs with the reduced pressure polynomial


    //! Create an object using the object construction function corresponding to
    //! "name" and the provided constructor arguments

    // template <class T, typename ... Args>

    // T* create(const string& name, Args... args) {
    //     return m_creators.at(canonicalize(name))(args...);
    // }

    // //! Register a new object construction function
    // void reg(const string& name, function<T*(Args...)> f) {
    //     m_creators[name] = f;
    // }

    //     // Troe falloff evaluator
    // reg("Troe", [](const AnyMap& node, const UnitStack& rate_units) {
    //     return new TroeRate(node, rate_units);
    // });

protected:
    std::unordered_map<string, function<T*(Args...)>> m_creators;
    double logP_ = -1000;
    // double logP1_ = 1000;
    // double logP2_ = -1000;
    // size_t ilow1_, ilow2_, ihigh1_, ihigh2_;
    // double rDeltaP_ = -1.0; //!< reciprocal of (logP2 - logP1)
};


}
#endif