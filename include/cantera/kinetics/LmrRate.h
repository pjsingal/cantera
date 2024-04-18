//! @file LmrRate.h
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LMRRATE_H
#define CT_LMRRATE_H
#include "cantera/kinetics/Arrhenius.h"
#include "cantera/base/Array.h"
// #include "cantera/kinetics/Reaction.h"

namespace Cantera{
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
   

    void setParametersPLOG(const AnyMap& node, const UnitStack& rate_units, int i, vector<AnyMap> colliders, string species);
    void setParametersTroe(const AnyMap& node, const UnitStack& rate_units, int i, vector<AnyMap> colliders, string species);
    void setParametersChebyshev(const AnyMap& node, const UnitStack& rate_units, int i, vector<AnyMap> colliders, string species);
    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    void getParametersPLOG(AnyMap& rateNode, const Units& rate_units) const;
    void getParametersTroe(AnyMap& rateNode, const Units& rate_units) const;
    void getParametersChebyshev(AnyMap& rateNode, const Units& rate_units) const;
    void getParameters(AnyMap& rateNode, const Units& rate_units) const;
    void getParameters(AnyMap& rateNode) const override {
        return getParameters(rateNode, Units(0));
    }

    double evalFromStructPLOG(const LmrData& shared_data);
    double evalFromStructTroe(const LmrData& shared_data);
    double evalFromStructChebyshev(const LmrData& shared_data);
    double evalFromStruct(const LmrData& shared_data);

    void validatePLOG(const string& equation, const Kinetics& kin);
    void validateTroe(const string& equation, const Kinetics& kin);
    void validateChebyshev(const string& equation, const Kinetics& kin);
    void validate(const string& equation, const Kinetics& kin) override;

    //Global member variables
    map<string,ArrheniusRate> eig0_;
    string fit_;
    vector<string> fitlist_; //stores the fit type (Troe, PLOG, or Chebyshev) for each reaction in sequential order
    string s_;
    map<string,ArrheniusRate> eig0_extra_; // for colliders where no fit data is specified
    

    // PLOG member variables
    map<string, map<double, pair<size_t, size_t>>> pressures_;
    map<string, vector<ArrheniusRate>> rates_;
    double logPeff_;
    double eig0_mix_ = 0.0;
    double log_eig0_mix_ = 0.0;
    double k_LMR_;

    //Troe member variables
    map<string,ArrheniusRate> k0_;
    map<string,ArrheniusRate> kinf_;
    map<string,double> m_a;//! parameter a in the 4-parameter Troe falloff function. Dimensionless
    map<string,double> m_rt3;//! parameter 1/T_3 in the 4-parameter Troe falloff function. [K^-1]
    map<string,double> m_rt1;//! parameter 1/T_1 in the 4-parameter Troe falloff function. [K^-1]
    map<string,double> m_t2;//! parameter T_2 in the 4-parameter Troe falloff function. [K]

    //Chebyshev member variables
    map<string,double> Tmin_, Tmax_; //!< valid temperature range [K]
    map<string,double> Pmin_, Pmax_; //!< valid pressure range [Pa]
    map<string,double> TrNum_, TrDen_; //!< terms appearing in the reduced temperature
    map<string,double> PrNum_, PrDen_; //!< terms appearing in the reduced pressure
    map<string,Array2D> m_coeffs_; //!<< coefficient array
    map<string,vector<double>> dotProd_; //!< dot product of coeffs with the reduced pressure polynomial

protected:
    double logP_ = -1000;
    double logP1_ = 1000;
    double logP2_ = -1000;
    size_t ilow1_, ilow2_, ihigh1_, ihigh2_;
    double rDeltaP_ = -1.0; //!< reciprocal of (logP2 - logP1)
};


}
#endif