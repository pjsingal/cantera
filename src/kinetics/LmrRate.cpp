//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/ChebyshevRate.h"
#include "cantera/kinetics/PlogRate.h"
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
    rate_units_=rate_units;
    ReactionRate::setParameters(node, rate_units);
    if(node.hasKey("fit") && node.hasKey("collider-list")){
        auto& colliders = node["collider-list"].asVector<AnyMap>();
        for (int i = 0; i < colliders.size(); i++){
            if (colliders[i].hasKey("collider") && colliders[i].hasKey("eig0")) {
                string species = colliders[i]["collider"].as<std::string>();
                eig0_.insert({species , ArrheniusRate(AnyValue(colliders[i]["eig0"]), node.units(), rate_units)});
                if(colliders[i].hasKey("rate-constants")){ //PLOG type            
                    plogList_.insert({species, colliders[i]});
                } else if(colliders[i].hasKey("Troe")){ //Troe type
                    troeList_.insert({species, colliders[i]});
                } else if(colliders[i].hasKey("data")&&colliders[i].hasKey("pressure-range")&&colliders[i].hasKey("temperature-range")){ //Chebyshev type
                    chebList_.insert({species, colliders[i]});
                }
            } else {
                throw InputFileError("LmrRate::setParameters", m_input,"An eig0 value must be provided for all colliders");
            }
        }
    } else {
        throw InputFileError("LmrRate::setParameters", m_input,"LMR-R data is incorrectly entered in yaml.");
    }
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    double k_LMR_=0.0;
    double k_M;
    double eig0;
    double eig0_mix_ = 0.0;
    double eig0_M=eig0_["M"].evalRate(shared_data.logT, shared_data.recipT);
    logP_=shared_data.logP; 
    double eig0_mix=0;
    vector<double> moleFractions = shared_data.moleFractions;
    for (size_t i=0; i<shared_data.allSpecies_.size(); i++){ //testing each species listed at the top of yaml file
        double Xi = moleFractions[i];
        std::map<string, ArrheniusRate>::iterator it = eig0_.find(shared_data.allSpecies_[i]);
        if (it != eig0_.end()) {//key found, i.e. species has at least an eig0 value explicitly specified in the LMR-R reaction in yaml  
            eig0_mix += Xi*eig0_[shared_data.allSpecies_[i]].evalRate(shared_data.logT, shared_data.recipT);
        }
        else {//no eig0 specified for this species, so use eig0 of species M as default
            eig0_mix += Xi*eig0_M;
        }
    }
    double log_eig0_mix = std::log(eig0_mix);
    double k;
    for (size_t i=0; i<shared_data.allSpecies_.size(); i++){ //testing each species listed at the top of yaml file
        string flag;
        double Xi = moleFractions[i];
        std::map<string, ArrheniusRate>::iterator it1 = eig0_.find(shared_data.allSpecies_[i]);
        std::map<string, AnyMap>::iterator it2 = plogList_.find(shared_data.allSpecies_[i]);
        eig0=eig0_[shared_data.allSpecies_[i]].evalRate(shared_data.logT, shared_data.recipT); 
        //Case 1: collider specified in PLOG format
        if (it2 != plogList_.end()){ 
            PlogData dataObj;
            PlogRate rateObj;
            AnyMap colliderNode = plogList_[shared_data.allSpecies_[i]];
            double logPeff_=logP_+log_eig0_mix-log(eig0);
            dataObj.logP=logPeff_; //use the effective pressure instead of P_abs for PLOG calculation
            rateObj.setParameters(colliderNode,rate_units_);
            k = rateObj.evalFromStruct(dataObj);
            if (shared_data.allSpecies_[i]=="M"){
                k_M = k;
            }
        }
        //Case 2: collider specified in Troe format
        else if (std::find(troeList_.begin(), troeList_.end(), shared_data.allSpecies_[i]) != troeList_.end()){ 
            FalloffData dataObj;
            FalloffRate rateObj;
            AnyMap colliderNode = troeList_[shared_data.allSpecies_[i]];
            rateObj.setParameters(colliderNode,rate_units_);
            k = rateObj.evalFromStruct(dataObj);
            if (shared_data.allSpecies_[i]=="M"){
                k_M = k;
            }
        }
        //Case 3: collider specified in Chebyshev format
        else if (std::find(chebList_.begin(), chebList_.end(), shared_data.allSpecies_[i]) != chebList_.end()){ 
            ChebyshevData dataObj;
            ChebyshevRate rateObj;
            AnyMap colliderNode = chebList_[shared_data.allSpecies_[i]];
            rateObj.setParameters(colliderNode,rate_units_);
            k = rateObj.evalFromStruct(dataObj);
            if (shared_data.allSpecies_[i]=="M"){
                k_M = k;
            }
        }
        //Case 4: only eig0 specified for collider, so treat it as same as "M" but with different eig0
        else if (it1 != eig0_.end()){ 
            k=k_M;
        }
        //Case 5: collider not specified in yaml, so treat it the same as "M"
        else { 
            k=k_M;
            eig0=eig0_["M"].evalRate(shared_data.logT, shared_data.recipT);
        }
        double Xtilde = eig0*Xi/eig0_mix; 
        k_LMR_ += k*Xtilde;
    }
    return k_LMR_;
}

// void LmrRate::getParameters(AnyMap& rateNode, const Units& rate_units) const{ //STILL NEED TO ACCOUNT FOR TROE AND PLOG DATA ENTRY
//     //Create copies of existing variables
//     map<string, map<double, pair<size_t, size_t>>> pressures__ = pressures_;
//     map<string, vector<ArrheniusRate>> rates__ = rates_;
//     map<string,ArrheniusRate> eig0__ = eig0_;
//     vector<string> colliderList;
//     for (const auto& entry : eig0__) {
//         colliderList.push_back(entry.first);
//     }
//     vector<AnyMap> topLevelList;
//     for (size_t i=0; i<colliderList.size(); i++) {
//         AnyMap speciesNode; //will be filled with all LMR data for a single collider (species)
//         if (!valid()) { //WHERE TO PUT THIS CONDITIONAL STATEMENT?
//             return;
//         }
//         string& s = colliderList[i];
//         //1) Save name of species to "name"
//         speciesNode["name"]=s;
//         //2) Save single set of arrhenius params for eig0 to "low-P-rate-constant"
//         AnyMap tempNode;
//         eig0__[s].getRateParameters(tempNode);
//         if (!tempNode.empty()){
//             speciesNode["eig0"]=std::move(tempNode);
//         }
//         tempNode.clear();
//         //3) Save list of rate constant params to "rate-constants"
//         std::multimap<double, ArrheniusRate> rateMap;
//         for (auto iter = ++pressures__[s].begin(); iter->first < 1000; ++iter) {
//             for (size_t i = iter->second.first; i < iter->second.second; i++) {
//                 rateMap.insert({std::exp(iter->first), rates__[s][i]});
//             }
//         }
//         vector<AnyMap> rateList;
//         for (const auto& [pressure, rate] : rateMap) {
//             tempNode["P"].setQuantity(pressure, "atm");
//             rate.getRateParameters(tempNode);
//             rateList.push_back(std::move(tempNode));
//             tempNode.clear();
//         }    
//         speciesNode["rate-constants"] = std::move(rateList);
//         topLevelList.push_back(std::move(speciesNode));
//     }
//     rateNode["collider-list"]=std::move(topLevelList);
// }

}

