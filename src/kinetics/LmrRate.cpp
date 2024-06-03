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
    double P = phase.pressure();
    int X = phase.stateMFNumber();
    if (allSpecies.empty()){
        allSpecies = phase.speciesNames(); //Get the list of all species in yaml (not just the ones for which LMRR data exists)
    }
    if (moleFractions.empty()){
        moleFractions.resize(allSpecies.size());
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
    if (m_pressure_buf < 0.) {
        return;
    }
    update(temperature,m_pressure_buf);
    m_pressure_buf = -1.;
}

LmrRate::LmrRate(const AnyMap& node, const UnitStack& rate_units){
    setParameters(node, rate_units);
}

void LmrRate::setParameters(const AnyMap& node, const UnitStack& rate_units){
    ReactionRate::setParameters(node, rate_units);
    rate_units_=rate_units;
    if(node.hasKey("collider-list")){
        auto& colliders = node["collider-list"].asVector<AnyMap>();
        for (int i = 0; i < colliders.size(); i++){
            if (colliders[i].hasKey("collider") && colliders[i].hasKey("eig0")) {
                colliderInfo.insert({colliders[i]["collider"].as<std::string>(), colliders[i]});
            } else {
                throw InputFileError("LmrRate::setParameters", m_input,"An eig0 value must be provided for all explicitly declared colliders in LMRR yaml entry.");
            }
        }
    } else {
        throw InputFileError("LmrRate::setParameters", m_input,"Yaml input for LMR-R does not follow the necessary structure.");
    }
}

void LmrRate::validate(const string& equation, const Kinetics& kin){}

double LmrRate::evalPlogRate(const LmrData& shared_data, double eig0val, map<string,AnyMap>::iterator it){
    PlogData plog_data;
    PlogRate plog_rate;
    plog_data.logP = shared_data.logP+log(eig0_mix)-log(eig0val); //replaces logP with log of the effective pressure w.r.t. eig0_M
    plog_data.logT = shared_data.logT;
    plog_data.pressure = shared_data.pressure;
    plog_data.recipT = shared_data.recipT;
    plog_data.temperature = shared_data.temperature;
    plog_rate.setParameters(it->second,rate_units_); //it->second refers to the yaml data for the ith collider
    plog_rate.updateFromStruct(plog_data);
    return plog_rate.evalFromStruct(plog_data);
}

double LmrRate::evalTroeRate(const LmrData& shared_data, double eig0val, map<string,AnyMap>::iterator it){
    FalloffData troe_data;
    TroeRate troe_rate;
    troe_data.conc_3b = shared_data.moleFractions;
    troe_data.logT = shared_data.logT;
    // troe_data.molar_density = shared_data.pressure; //
    troe_data.ready = shared_data.ready;
    troe_data.recipT = shared_data.recipT;
    troe_data.temperature = shared_data.temperature;
    // colliders_i["type"]; 
    troe_rate.setParameters(it->second,rate_units_); //it->second refers to the yaml data for the ith collider
    return troe_rate.evalFromStruct(troe_data);
}

double LmrRate::evalChebyshevRate(const LmrData& shared_data, double eig0val, map<string,AnyMap>::iterator it){
    ChebyshevData cheb_data;
    ChebyshevRate cheb_rate;
    cheb_data.log10P=shared_data.logP;
    cheb_data.logT=shared_data.logT;
    cheb_data.pressure=shared_data.pressure;
    cheb_data.recipT=shared_data.recipT;
    cheb_data.temperature=shared_data.temperature;
    cheb_rate.setParameters(it->second,rate_units_); //it->second refers to the yaml data for the ith collider
    cheb_rate.updateFromStruct(cheb_data);
    return cheb_rate.evalFromStruct(cheb_data);
}

double LmrRate::geteig0(const LmrData& shared_data, map<string,AnyMap>::iterator it){
    // ArrheniusRate eig0_i_ = ArrheniusRate(AnyValue(colliders[i]["low-P-rate-constant"]), node.units(), rate_units_);
    // double eig0 = eig0_[s].evalRate(log(T[i]), 1.0/T[i]);
    ArrheniusRate eig = ArrheniusRate(AnyValue(it->second["eig0"]), it->second.units(), rate_units_);
    return eig.evalRate(shared_data.logT, shared_data.recipT);
}

double LmrRate::geteig0mix(const LmrData& shared_data){
    double eig0mix=0.0;
    for (size_t i=0; i<shared_data.allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        auto it = colliderInfo.find(shared_data.allSpecies[i]);
        if (it != colliderInfo.end()) { //yaml species has at least an eig0 value provided for LMRR 
            eig0mix += shared_data.moleFractions[i]*geteig0(shared_data,it);
        }
        else if (it == colliderInfo.end()) { //yaml species has no data provided for LMRR (treat as "M")
            eig0mix += shared_data.moleFractions[i]*eig0_M;
        }
        else{
            throw InputFileError("LmrRate::geteig0mix", m_input,"Cannot compute eig0mix due to invalid LMRR yaml input.");
        }
    }
    writeMsg("eig0mix = ",eig0mix);
    return eig0mix;
}

vector<double> LmrRate::get_eig0M_kM(const LmrData& shared_data){
    double kM;
    double eig0M;
    auto it = colliderInfo.find("M");
    if (it != colliderInfo.end() && it->second.hasKey("rate-constants")){ 
        eig0M=ArrheniusRate(AnyValue(it->second["eig0"]), it->second.units(), rate_units_).evalRate(shared_data.logT, shared_data.recipT);
        kM = evalPlogRate(shared_data, eig0M,it);
        writeMsg("eig0_M = ",eig0M);
        writeMsg("kM_plog = ",kM);
    } 
    else if(it != colliderInfo.end() && it->second.hasKey("Troe")){ 
        eig0M=ArrheniusRate(AnyValue(it->second["eig0"]), it->second.units(), rate_units_).evalRate(shared_data.logT, shared_data.recipT);
        kM = evalTroeRate(shared_data, eig0M,it);
        writeMsg("eig0_M = ",eig0M);
        writeMsg("kM_troe = ",kM);
        // colliders_i["type"]="LMR_R"; //this dummy key needs to be defined because falloff.cpp requires it
    }
    else if(it != colliderInfo.end() && it->second.hasKey("pressure-range")){
        eig0M=ArrheniusRate(AnyValue(it->second["eig0"]), it->second.units(), rate_units_).evalRate(shared_data.logT, shared_data.recipT);
        kM = evalChebyshevRate(shared_data, eig0M,it);
        writeMsg("eig0_M = ",eig0M);
        writeMsg("kM_cheb = ",kM);
    }
    else {
        throw InputFileError("LmrRate::getkM", m_input,"Not enough data provided for species 'M'.");
    }
    return {eig0M,kM};
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    // auto it = colliderInfo.find("M");
    // if (it != colliderInfo.end()){ 
    //     eig0_M = geteig0(shared_data,it);
    // }
    vector<double> Mvals = get_eig0M_kM(shared_data);
    eig0_M = Mvals[0];
    k_M = Mvals[1];

    eig0_mix = geteig0mix(shared_data);
    writeMsg("k_M_overall = ",k_M);
    k_LMR=0.0;
    for (size_t i=0; i<shared_data.allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        auto it = colliderInfo.find(shared_data.allSpecies[i]);
        // if (it != colliderInfo.end() && it->second.hasKey("rate-constants")){ 
        // double k_i;
        if (it != colliderInfo.end() && it->second.hasKey("rate-constants")){ 
            writelog(it->first);writelog("\n"); //speciesID
            eig0_i = geteig0(shared_data,it);
            k_i = evalPlogRate(shared_data, eig0_i,it);
            writeMsg("eig0_i_plog = ",eig0_i); writeMsg("k_i_plog = ",k_i);
        } 
        else if(it != colliderInfo.end() && it->second.hasKey("Troe")){ 
            writelog(it->first);writelog("\n"); //speciesID
            // colliders_i["type"]="LMR_R"; //this dummy key needs to be defined because falloff.cpp requires it
            eig0_i = geteig0(shared_data,it);
            k_i = evalTroeRate(shared_data, eig0_i,it);
            writeMsg("eig0_i_troe = ",eig0_i); writeMsg("k_i_troe = ",k_i);
        }
        else if(it != colliderInfo.end() && it->second.hasKey("pressure-range")){ 
            writelog(it->first);writelog("\n"); //speciesID
            eig0_i = geteig0(shared_data,it);
            k_i = evalChebyshevRate(shared_data, eig0_i,it);
            writeMsg("eig0_i_cheb = ",eig0_i); writeMsg("k_i_cheb = ",k_i);
        }
        else if(it != colliderInfo.end() && !(it->second.hasKey("pressure-range")) && !(it->second.hasKey("Troe")) && !(it->second.hasKey("rate-constants"))){ //yaml species has an eig0 but no additional LMRR data, so treat its rate as same as "M"
            writelog(it->first);writelog("\n"); //speciesID
            eig0_i = geteig0(shared_data,it);
            k_i=k_M;
            writeMsg("eig0_i_eigOnly = ",eig0_i); writeMsg("k_i_eigOnly = ",k_i);
        }
        else if(it == colliderInfo.end()){ //yaml species has no LMRR data, so treat its rate and eig0 as same as "M"
            k_i=k_M;
            eig0_i=eig0_M;
            writeMsg("eig0_i_noLMR = ",eig0_i); writeMsg("k_i_noLMR = ",k_i);
        }
        else{
            throw InputFileError("LmrRate::evalFromStruct", m_input,"LMRR reaction has invalid yaml input.");
        }
        writeMsg("k_i_final = ",k_i);
        k_LMR += k_i*eig0_i*shared_data.moleFractions[i]/eig0_mix; //Note: Xtilde = eig0_i*shared_data.moleFractions[i]/eig0_mix;
    }
    writeMsg("k_LMR = ",k_LMR);
    return k_LMR;
}

void LmrRate::getParameters(AnyMap& rateNode, const Units& rate_units) const{ //STILL NEED TO ACCOUNT FOR TROE AND PLOG DATA ENTRY
    vector<AnyMap> topLevelList;
    for (const auto& entry : colliderInfo) {
        string name = entry.first;
        auto colliders_i = entry.second;
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

