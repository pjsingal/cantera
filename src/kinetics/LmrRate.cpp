//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/ChebyshevRate.h"
#include "cantera/kinetics/PlogRate.h"
// #include "cantera/kinetics/ReactionRateFactory.h"
// #include "cantera/base/FactoryBase.h"
// #include "cantera/kinetics/Reaction.h"
// #include <sstream>

namespace Cantera{

LmrData::LmrData(){ //THIS METHOD WAS ADAPTED SOMEWHAT BLINDLY FROM FALLOFF.CPP, PLEASE VERIFY IF CORRECT
    moleFractions.resize(1, NAN);
}

bool LmrData::update(const ThermoPhase& phase, const Kinetics& kin){
    double T = phase.temperature();
    double P = phase.pressure(); //find out what units this is in
    int X = phase.stateMFNumber();

    //Get the list of all species in yaml (not just the ones for which LMRR data exists)
    if (allSpecies_.empty()){
        allSpecies_ = phase.speciesNames();
    }
    if (moleFractions.empty()){
        moleFractions.resize(allSpecies_.size());
    }
    if (P != pressure || T != temperature || X != mfNumber) {
        ReactionData::update(T);
        pressure = P;
        logP = std::log(P);
        mfNumber=X;
        phase.getMoleFractions(moleFractions.data());
        return true;
    }
    return false;
}

void LmrData::perturbPressure(double deltaP){
    if (m_pressure_buf > 0.) {
        throw CanteraError("LmrData::perturbPressure",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_pressure_buf = pressure;
    update(temperature,pressure*(1. + deltaP));
}

void LmrData::restore(){
    ReactionData::restore();
    // only restore if there is a valid buffered value
    if (m_pressure_buf < 0.) {
        return;
    }
    update(temperature,m_pressure_buf);
    m_pressure_buf = -1.;
}

// Methods of class LmrRate
LmrRate::LmrRate(const AnyMap& node, const UnitStack& rate_units){
    setParameters(node, rate_units);
}

void LmrRate::setParameters(const AnyMap& node, const UnitStack& rate_units){
    ReactionRate::setParameters(node, rate_units);
    if(node.hasKey("collider-list")){
        auto& colliders = node["collider-list"].asVector<AnyMap>();
        for (int i = 0; i < colliders.size(); i++){
            if (colliders[i].hasKey("collider") && colliders[i].hasKey("eig0")) {
                colliderInfo.insert({colliders[i]["collider"].as<std::string>(), std::make_pair(colliders[i],rate_units)});
            } else {
                throw InputFileError("LmrRate::setParameters", m_input,"An eig0 value must be provided for all explicitly declared colliders in LMRR yaml entry.");
            }
        }
    } else {
        throw InputFileError("LmrRate::setParameters", m_input,"Yaml input for LMR-R does not follow the necessary structure.");
    }
}

void LmrRate::validate(const string& equation, const Kinetics& kin)
{
    if (!valid()) {
        throw InputFileError("LmrRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
    }
}

void LmrRate::syncPlogData(const LmrData& shared_data, PlogData& plog_data){
    plog_data.logP = shared_data.logP;
    plog_data.logT = shared_data.logT;
    plog_data.pressure = shared_data.pressure;
    plog_data.recipT = shared_data.recipT;
    plog_data.temperature = shared_data.temperature;
}

void LmrRate::syncTroeData(const LmrData& shared_data, FalloffData& troe_data){
    troe_data.conc_3b = shared_data.moleFractions;
    troe_data.logT = shared_data.logT;
    // troe_data.molar_density = shared_data.pressure; //
    troe_data.ready = shared_data.ready;
    troe_data.recipT = shared_data.recipT;
    troe_data.temperature = shared_data.temperature;
}

void LmrRate::syncChebData(const LmrData& shared_data, ChebyshevData& cheb_data){
    cheb_data.log10P=shared_data.logP;
    cheb_data.logT=shared_data.logT;
    cheb_data.pressure=shared_data.pressure;
    cheb_data.recipT=shared_data.recipT;
    cheb_data.temperature=shared_data.temperature;
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    //STEP 0: CALL PARAMS FROM SHARED DATA STRUCT
    vector<double> X = shared_data.moleFractions;
    double logT = shared_data.logT;
    double recipT = shared_data.recipT;
    logP_ = shared_data.logP;
    vector<string> yamlSpecies = shared_data.allSpecies_;

    //STEP 1: GET EIG0_M
    double eig0_M;
    map<string, pair<const AnyMap&,const UnitStack&>>::iterator it1 = colliderInfo.find("M");
        if (it1 != colliderInfo.end()){ 
            const AnyMap& colliders_i = it1->second.first;
            const UnitStack& rate_units_i = it1->second.second;
            eig0_M=ArrheniusRate(AnyValue(colliders_i["eig0"]), colliders_i.units(), rate_units_i).evalRate(logT, recipT);
    }

    //STEP 2: GET EIG0_MIX
    double eig0_mix=0.0;
    for (size_t i=0; i<yamlSpecies.size(); i++){ //testing each species listed at the top of yaml file
        map<string, pair<const AnyMap&,const UnitStack&>>::iterator it2 = colliderInfo.find(yamlSpecies[i]);
        //Case 2.1: yaml species has at least an eig0 value provided for LMRR 
        if (it2 != colliderInfo.end()) {
            const AnyMap& colliders_i = it2->second.first;
            const UnitStack& rate_units_i = it2->second.second;
            eig0_mix += X[i]*ArrheniusRate(AnyValue(colliders_i["eig0"]), colliders_i.units(), rate_units_i).evalRate(logT, recipT);
        }
        //Case 2.2: yaml species has no data provided for LMRR (treat as "M")
        else {
            eig0_mix += X[i]*eig0_M;
        }
    }


    //STEP 3: DEFINE DATA OBJECTS FOR PLOG, TROE, AND CHEBYSHEV
    
    PlogData plog_data;
    FalloffData troe_data;
    ChebyshevData cheb_data;
    syncPlogData(shared_data, plog_data);
    syncTroeData(shared_data, troe_data);
    syncChebData(shared_data, cheb_data);
    
    //STEP 4: GET K FOR THE GENERIC COLLIDER 'M'
    double k_M;
    map<string, pair<const AnyMap&,const UnitStack&>>::iterator it3 = colliderInfo.find("M");
    if (it3 != colliderInfo.end()){ 
        const AnyMap& colliders_i = it3->second.first;
        const UnitStack& rate_units_i = it3->second.second;
        //Case 4.1: "M" is of PLOG type
        if(colliders_i.hasKey("rate-constants")){ 
            
            plog_data.logP = logP_+log(eig0_mix)-log(eig0_M); //CRUCIAL STEP: replace logP with log of the effective pressure w.r.t. eig0_M
            
            //calling PlogRate(colliders_i,rate_units_i); is basically just calling setParameters within PlogRate
            // newobj.
            
            k_M = PlogRate(colliders_i,rate_units_i).evalFromStruct(plog_data);
            // k_M = 0.005;
        } 
        //Case 4.2: "M" is of Troe type
        else if(colliders_i.hasKey("Troe")){ 
            k_M = TroeRate(colliders_i,rate_units_i).evalFromStruct(troe_data);
            // k_M = 0.005;
        } 
        //Case 4.3: "M" is of Chebyshev type
        else if(colliders_i.hasKey("data")&&colliders_i.hasKey("pressure-range")&&colliders_i.hasKey("temperature-range")){ 
            k_M = ChebyshevRate(colliders_i,rate_units_i).evalFromStruct(cheb_data);
            // k_M = 0.005;
        } 
        //Case 4.4: insufficient data has been provided for "M"
        else {
            throw InputFileError("LmrRate::evalFromStruct", m_input,"Not enough data provided for species 'M'.");
        }
    }

    //STEP 5: GET K AND EIG0 FOR ALL SPECIES LISTED AT TOP OF YAML FILE AND THEN COMPUTE K_LMR
    double k_LMR=0.0;
    for (size_t i=0; i<yamlSpecies.size(); i++){ //testing each species listed at the top of yaml file
        double k;
        double eig0;
        map<string, pair<const AnyMap&,const UnitStack&>>::iterator it4 = colliderInfo.find(yamlSpecies[i]);
        //Case 5.1: yaml species has at least an eig0 value provided for LMRR
        if (it4 != colliderInfo.end()){ 
            const AnyMap& colliders_i = it4->second.first;
            const UnitStack& rate_units_i = it4->second.second;
            eig0 = ArrheniusRate(AnyValue(colliders_i["eig0"]), colliders_i.units(), rate_units_i).evalRate(logT, recipT);
            //Case 5.1.1: yaml species has LMRR data provided in PLOG format
            if(colliders_i.hasKey("rate-constants")){
                plog_data.logP = logP_+log(eig0_mix)-log(eig0); //CRUCIAL STEP: replace logP with log of the effective pressure w.r.t. eig0
                k = PlogRate(colliders_i,rate_units_i).evalFromStruct(plog_data);
                // k=0.005;
            } 
            //Case 5.1.2: yaml species has LMRR data provided in Troe format
            else if(colliders_i.hasKey("Troe")){ //"M" is of Troe type
                k = TroeRate(colliders_i,rate_units_i).evalFromStruct(troe_data);
                // k=0.005;
            } 
            //Case 5.1.3: yaml species has LMRR data provided in Chebyshev format
            else if(colliders_i.hasKey("data")&&colliders_i.hasKey("pressure-range")&&colliders_i.hasKey("temperature-range")){ //"M" is of Chebyshev type
                k = ChebyshevRate(colliders_i,rate_units_i).evalFromStruct(cheb_data);
                // k=0.005;
            } 
            //Case 5.1.4: yaml species has an eig0 but no additional LMRR data, so treat its rate as same as "M"
            else {
                k = k_M;
            }
        }
        //Case 5.2: yaml species has no LMRR data, so treat its rate and eig0 as same as "M"
        else{
            k=k_M;
            eig0=eig0_M;
        }
        double Xtilde = eig0*X[i]/eig0_mix; 
        k_LMR += k*Xtilde;
    }
    return k_LMR;
}

void LmrRate::getParameters(AnyMap& rateNode, const Units& rate_units) const{ //STILL NEED TO ACCOUNT FOR TROE AND PLOG DATA ENTRY
    vector<AnyMap> topLevelList;
    for (const auto& entry : colliderInfo) {
        string name = entry.first;
        const AnyMap& colliders_i = entry.second.first;
        const UnitStack& rate_units_i = entry.second.second;
        AnyMap colliderNode;
        if(colliders_i.hasKey("rate-constants")){
            colliderNode["collider"]=name;
            colliderNode["eig0"]=colliders_i["eig0"];
            colliderNode["rate-constants"]=colliders_i["rate-constants"];
        } 
        else if(colliders_i.hasKey("Troe")){
            colliderNode["collider"]=name;
            colliderNode["eig0"]=colliders_i["eig0"];
            colliderNode["low-P-rate-constant"]=colliders_i["low-P-rate-constant"];
            colliderNode["high-P-rate-constant"]=colliders_i["high-P-rate-constant"];
            colliderNode["Troe"]=colliders_i["Troe"];
        } 
        else if(colliders_i.hasKey("data")&&colliders_i.hasKey("pressure-range")&&colliders_i.hasKey("temperature-range")){ //"M" is of Chebyshev type
            colliderNode["collider"]=name;
            colliderNode["eig0"]=colliders_i["eig0"];
            colliderNode["temperature-range"]=colliders_i["temperature-range"];
            colliderNode["pressure-range"]=colliders_i["pressure-range"];
            colliderNode["data"]=colliders_i["data"];
        } 
        else {
            colliderNode["collider"]=name;
            colliderNode["eig0"]=colliders_i["eig0"];
        }
        topLevelList.push_back(std::move(colliderNode));
    }
    rateNode["collider-list"]=std::move(topLevelList);
}

}

