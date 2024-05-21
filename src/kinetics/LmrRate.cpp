//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
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
    // writelog("T = {}\n", T);
    // writelog("P = {}\n", P);
    // writelog("X = {}\n", X);

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
    //writelog("m_pressure_buf = {}\n", m_pressure_buf);
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
    // Reaction reactionInstance;
    // size_t numReactants = reactionInstance.reactants.size();
    // UnitStack rate_units_ = rate_units;
    // if (numReactants>0) {
    //     rate_units_.join(1);
    // }
    // writelog("numReactants = {}\n", numReactants);

    setParameters(node, rate_units);
}

// void BlowersMaselRate::setContext(const Reaction& rxn, const Kinetics& kin)
// {
//     m_stoich_coeffs.clear();
//     for (const auto& [name, stoich] : rxn.reactants) {
//         m_stoich_coeffs.emplace_back(kin.kineticsSpeciesIndex(name), -stoich);
//     }
//     for (const auto& [name, stoich] : rxn.products) {
//         m_stoich_coeffs.emplace_back(kin.kineticsSpeciesIndex(name), stoich);
//     }
// }

void LmrRate::setParameters(const AnyMap& node, const UnitStack& rate_units){
    std::string rxn = node["equation"].as<std::string>();
    vector<string> reactantList;
    string component = "";
    for (int i = 0; i < rxn.length(); i++) {
        if (!(rxn[i]=='(') && !(rxn[i]==')') && !(rxn[i]=='+') && !(rxn[i]==' ') && !(rxn[i]=='<') && !(rxn[i]=='=') && !(rxn[i]=='>')){
            component+=rxn[i];
        } else if ((i>0 && rxn[i]==' ' && !(rxn[i-1]==')') && !(rxn[i-1]=='+'))  || rxn[i]==')') {
            reactantList.push_back(component);
            // writelog("component = {}\n", component);
            component = "";
        } else if (rxn[i]=='<'){
            break;
        }
    } 
    UnitStack rate_units_ = rate_units;
    if (reactantList.size()>2){
        rate_units_.join(1);
    }
    // writelog("numReactants = {}\n", reactantList.size());
    ReactionRate::setParameters(node, rate_units_);
    if (node.hasKey("collider-list")) {
        // writelog_direct("TRUE1\n");
        auto& colliders = node["collider-list"].asVector<AnyMap>();
        for (int i = 0; i < colliders.size(); i++){
            if (colliders[i].hasKey("collider") && colliders[i].hasKey("low-P-rate-constant")) {
                if(colliders[i].hasKey("rate-constants")){
                    string species_i_ = colliders[i]["collider"].as<std::string>();
                    // writelog("species_i_ = {}\n", species_i_);
                    ArrheniusRate eig0_i_ = ArrheniusRate(AnyValue(colliders[i]["low-P-rate-constant"]), node.units(), rate_units_);       
                    map<double, pair<size_t, size_t>> pressures_i_;
                    vector<ArrheniusRate> rates_i_; 
                    std::multimap<double, ArrheniusRate> multi_rates;
                    auto& rates = colliders[i]["rate-constants"].asVector<AnyMap>();
                    for (const auto& rate : rates){
                        multi_rates.insert({rate.convert("P","Pa"),ArrheniusRate(AnyValue(rate), node.units(), rate_units_)});
                    }
                    rates_i_.reserve(multi_rates.size());
                    m_valid = !multi_rates.empty(); //if rates object empty, m_valid==FALSE. if rates is not empty, m_valid==TRUE
                    size_t j = 0;
                    for (const auto& [pressure, rate] : multi_rates) { 
                        double logp = std::log(pressure);
                        if (pressures_i_.empty() || pressures_i_.rbegin()->first != logp) {
                            pressures_i_[logp] = {j, j+1};
                        } else {
                            pressures_i_[logp].second = j+1;
                        }
                        j++;
                        rates_i_.push_back(rate); 
                    }
                    if (!m_valid) { //runs if multi_rates is empty
                        // ensure that reaction rate can be evaluated (but returns NaN)
                        rates_i_.reserve(1);
                        pressures_i_[std::log(OneBar)] = {0, 0};
                        rates_i_.push_back(ArrheniusRate());
                    }
                    // Duplicate the first and last groups to handle P < P_0 and P > P_N
                    pressures_i_.insert({-1000.0, pressures_i_.begin()->second});
                    pressures_i_.insert({1000.0, pressures_i_.rbegin()->second});

                    rates_.insert({species_i_, rates_i_});
                    pressures_.insert({species_i_, pressures_i_});
                    eig0_.insert({species_i_, eig0_i_});
                }
                else{
                    string species_i_extra_ = colliders[i]["collider"].as<std::string>();
                    ArrheniusRate eig0_i_extra_ = ArrheniusRate(AnyValue(colliders[i]["low-P-rate-constant"]), node.units(), rate_units_);
                    eig0_extra_.insert({species_i_extra_, eig0_i_extra_});
                }
            }
        }
    }
}

void LmrRate::validate(const string& equation, const Kinetics& kin){
    //Validate the LMRR input data for each species
    if (!valid()) {
        throw InputFileError("LmrRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
    }
    fmt::memory_buffer err_reactions1; //for k-related errors
    fmt::memory_buffer err_reactions2; //for eig0-related errors
    fmt::memory_buffer err_reactions3; //for eig0_extra-related errors
    double T[] = {300.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
    // LmrData data;
    // Iterate through the outer map (string to inner map)
    // writelog_direct("TRUE0");
    for (const auto& outer_pair : pressures_) { //iterating through species
        // writelog_direct("TRUE1");
        const std::string& s = outer_pair.first; //s refers only to the species for which LMR data is provided in yaml (e.g. 'H2O', 'M')
        const std::map<double, std::pair<size_t, size_t>>& inner_map = outer_pair.second;
        for (auto iter = ++inner_map.begin(); iter->first < 1000; iter++) { //iterating through pressures (and their corresponding arrhenius data)
            // writelog_direct("TRUE2");
            ilow1_ = iter->second.first;
            ilow2_ = iter->second.second;       
            for (size_t i=0; i < 6; i++) { //iterating through our sample temperature array
                // writelog_direct("TRUE3");
                double k = 0;
                for (size_t p = ilow1_; p < ilow2_; p++) {
                    k += rates_[s].at(p).evalRate(log(T[i]), 1.0 / T[i]);
                }
                double eig0 = eig0_[s].evalRate(log(T[i]), 1.0/T[i]);
                // writelog("eig0={}\n",eig0);
                // writelog("k={}\n",k);
                if (!(k > 0)){ //flags error if k at a given T, P is not > 0
                    // writelog_direct("TRUE4a");
                    fmt_append(err_reactions1,"at P = {:.5g}, T = {:.1f}\n", std::exp(iter->first), T[i]);
                }
                else if (!(eig0>0)){ //flags error if eig0 at a given T is not > 0
                    // writelog_direct("TRUE4b");
                    fmt_append(err_reactions2,"at T = {:.1f}\n", T[i]);
                }
                auto it = eig0_extra_.find(s);
                if (it != eig0_extra_.end()) {
                    double eig0_extra = eig0_extra_[s].evalRate(log(T[i]), 1.0/T[i]);
                    if (!(eig0_extra > 0)){ //flags error if k at a given T, P is not > 0 
                        fmt_append(err_reactions3,"T = {:.1f}\n", std::exp(iter->first), T[i]);
                    }
                }
            }
        }
        if (err_reactions1.size()) {
            throw InputFileError("LmrRate::validate", m_input,
                    "\nInvalid rate coefficient, k, for reaction '{}'\n{}",equation, to_string(err_reactions1));
        }
        else if (err_reactions2.size()) {
            throw InputFileError("LmrRate::validate", m_input,
                    "\nInvalid rate coefficient, eig0 (k at low-pressure limit), for reaction '{}'\n{}",equation, to_string(err_reactions2));
        }
        else if (err_reactions3.size()) {
            throw InputFileError("LmrRate::validate", m_input,
                    "\nInvalid rate coefficient, eig0_extra (k at low-pressure limit), which represents one of the colliders without a rate constant table.");
        }
    }
}

double LmrRate::speciesPlogRate(const LmrData& shared_data){
    // auto iter = pressures_s_.upper_bound(logPeff_);
    auto iter = pressures_s_.upper_bound(logP_);
    // AssertThrowMsg(iter != pressures_s_.end(), "LmrRate::speciesPlogRate","Reduced-pressure out of range: {}", logPeff_);
    // AssertThrowMsg(iter != pressures_s_.begin(), "LmrRate::speciesPlogRate","Reduced-pressure out of range: {}", logPeff_); 
    AssertThrowMsg(iter != pressures_s_.end(), "LmrRate::speciesPlogRate","Log-Pressure out of range: {}", logP_);
    AssertThrowMsg(iter != pressures_s_.begin(), "LmrRate::speciesPlogRate","Log-Pressure out of range: {}", logP_); 
    logP2_ = iter->first;
    ihigh1_ = iter->second.first;
    ihigh2_ = iter->second.second;
    logP1_ = (--iter)->first;
    ilow1_ = iter->second.first;
    ilow2_ = iter->second.second;
    rDeltaP_ = 1.0 / (logP2_ - logP1_);
    double log_k1, log_k2;
    // writelog("log_k1 = {}\n", log_k1);
    if (ilow1_ == ilow2_) {
        log_k1 = rates_s_[ilow1_].evalLog(shared_data.logT, shared_data.recipT);
    } else {
        double k = 1e-300;
        for (size_t i = ilow1_; i < ilow2_; i++) {
            k += rates_s_[i].evalRate(shared_data.logT, shared_data.recipT);
        }
        log_k1 = std::log(k);
    }
    if (ihigh1_ == ihigh2_) {
        log_k2 = rates_s_[ihigh1_].evalLog(shared_data.logT, shared_data.recipT);
    } else {
        double k = 1e-300;
        for (size_t i = ihigh1_; i < ihigh2_; i++) {
            
            k += rates_s_[i].evalRate(shared_data.logT, shared_data.recipT);
        }
        log_k2 = std::log(k);
    }
    return exp(log_k1 + (log_k2-log_k1)*rDeltaP_*(logPeff_-logP1_));
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    k_LMR_=0;
    double eig0_mix=0;



    // double eig0_M=eig0_["M"].evalRate(shared_data.logT, shared_data.recipT);
    //writelog("eig0_M = {}\n", eig0_M);
    vector<double> moleFractions = shared_data.moleFractions;
    
    for (size_t i=0; i<shared_data.allSpecies_.size(); i++){ //testing each species listed at the top of yaml file
        double Xi = moleFractions[i];
        // writelog("Xi = {}\n", Xi);
        std::map<string, ArrheniusRate>::iterator it1 = eig0_.find(shared_data.allSpecies_[i]);
        std::map<string, ArrheniusRate>::iterator it2 = eig0_extra_.find(shared_data.allSpecies_[i]);
        if ((it1 != eig0_.end()) && (it2 == eig0_extra_.end())) {//key found, i.e. species has corresponding LMR data  
            // writelog("1. eig0 contributed to eig0mix = {}\n",Xi*eig0_[shared_data.allSpecies_[i]].evalRate(shared_data.logT, shared_data.recipT));
            // writelog("1. species = {}\n",shared_data.allSpecies_[i]);
            eig0_mix += Xi*eig0_[shared_data.allSpecies_[i]].evalRate(shared_data.logT, shared_data.recipT);
        } 
        else if ( (it1 == eig0_.end()) && (it2 != eig0_extra_.end()) ) {
            // writelog("2. eig0_extra contributed to eig0mix = {}\n",Xi*eig0_extra_[shared_data.allSpecies_[i]].evalRate(shared_data.logT, shared_data.recipT));
            // writelog("2. species = {}\n",shared_data.allSpecies_[i]);
            eig0_mix += Xi*eig0_extra_[shared_data.allSpecies_[i]].evalRate(shared_data.logT, shared_data.recipT);
        }
        else {//no LMR data for this species, so use M data as default
            // writelog("3. eig0_M contributed to eig0mix = {}\n",Xi*eig0_["M"].evalRate(shared_data.logT, shared_data.recipT));
            // writelog("3. species = {}\n",shared_data.allSpecies_[i]);
            eig0_mix += Xi*eig0_["M"].evalRate(shared_data.logT, shared_data.recipT);;
        }
    }
    // writelog("OVERALL eig0_mix = {}\n",eig0_mix);
    // writelog("eig0_M = {}\n",eig0_M);
    // writelog("eig0_mix = {}\n",eig0_mix);
    // writelog("eig0_mix = {}\n", eig0_mix);
    double log_eig0_mix = std::log(eig0_mix);
    //writelog("log_eig0_mix = {}\n", log_eig0_mix);
    
    double eig0; //eig0 val of a single species

    double eig0_M = eig0_["M"].evalRate(shared_data.logT, shared_data.recipT);

    pressures_s_=pressures_["M"];
    rates_s_=rates_["M"];
    logP_=shared_data.logP; 
    logPeff_=logP_+log_eig0_mix-log(eig0_M);
    double k_M = LmrRate::speciesPlogRate(shared_data); //single-component rate constant of species M    k_m(T,R_LMR)

    for (size_t i=0; i<shared_data.allSpecies_.size(); i++){ //testing each species listed at the top of yaml file
        double Xi = moleFractions[i];
        //writelog("2 Xi={}\n", Xi);
        // logP_=shared_data.logP;
        // logPeff_=logP_+log_eig0_mix-log(eig0_M);
        
        // std::map<string, ArrheniusRate>::iterator it1 = eig0_.find(shared_data.allSpecies_[i]);
        std::map<string, ArrheniusRate>::iterator it2 = eig0_extra_.find(shared_data.allSpecies_[i]);
        // if ( (it1 != eig0_.end()) && (it2 == eig0_extra_.end()) ) { //species is LMRR-specified with a PLOG table
        //     eig0 = eig0_[shared_data.allSpecies_[i]].evalRate(shared_data.logT, shared_data.recipT);
        //     // writelog("1. eig0 = {}\n",eig0);
        //     pressures_s_=pressures_[shared_data.allSpecies_[i]];
        //     rates_s_=rates_[shared_data.allSpecies_[i]];
        //     // logP_=shared_data.logP; 
        //     // logPeff_=logP_+log_eig0_mix-log(eig0);
        //     k_LMR_ += LmrRate::speciesPlogRate(shared_data)*eig0*Xi/eig0_mix;
        // } 
        if (it2 != eig0_extra_.end()) { //species is LMRR-specified without a PLOG table
            double eig0 = eig0_extra_[shared_data.allSpecies_[i]].evalRate(shared_data.logT, shared_data.recipT);
            // writelog("2. eig0 = {}\n",eig0);
            // for (const auto& entry : pressures_["M"]) {
            //     double modifiedKey = entry.first*30000;
            //     pressures_s_[modifiedKey] = entry.second;
            // rates_s_=rates_["M"];
            // logP_=shared_data.logP; 
            // logPeff_=logP_+log_eig0_mix-log(eig0);
            // k_LMR_ += LmrRate::speciesPlogRate(shared_data)*eig0*Xi/eig0_mix;
            double Xtilde = eig0*Xi/eig0_mix; 
            k_LMR_ += k_M*Xtilde;
            // writelog("2. k_M*Xtilde = {}\n",k_M*Xtilde);
        }
        else { //species is not LMRR-specified; use generic parameters for species "M"
            double eig0 = eig0_["M"].evalRate(shared_data.logT, shared_data.recipT);
            // writelog("3. eig0 = {}\n",eig0);
            // pressures_s_=pressures_["M"];
            // rates_s_=rates_["M"];
            // logP_=shared_data.logP; 
            // logPeff_=logP_+log_eig0_mix-log(eig0);
            // k_LMR_ += LmrRate::speciesPlogRate(shared_data)*eig0_M*Xi/eig0_mix;
            double Xtilde = eig0*Xi/eig0_mix;
            k_LMR_ += k_M*Xtilde;
            // writelog("3. k_M*Xtilde; = {}\n",k_M*Xtilde);
        }
        // writelog("logPeff_ = {}\n",logPeff_);
        // writelog("shared_data.allSpecies_[i] = {}\n",shared_data.allSpecies_[i]);
        // writelog("k_LMR_ (cumulative) = {}\n",k_LMR_);
    }
    // writelog("\n\n\n",logPeff_);
    return k_LMR_;
}

void LmrRate::getParameters(AnyMap& rateNode, const Units& rate_units) const{
    //Create copies of existing variables
    map<string, map<double, pair<size_t, size_t>>> pressures__ = pressures_;
    map<string, vector<ArrheniusRate>> rates__ = rates_;
    map<string,ArrheniusRate> eig0__ = eig0_;
    vector<string> colliderList;
    for (const auto& entry : eig0__) {
        colliderList.push_back(entry.first);
    }
    vector<AnyMap> topLevelList;
    for (size_t i=0; i<colliderList.size(); i++) {
        AnyMap speciesNode; //will be filled with all LMR data for a single collider (species)
        if (!valid()) { //WHERE TO PUT THIS CONDITIONAL STATEMENT?
            return;
        }
        string& s = colliderList[i];
        //1) Save name of species to "name"
        speciesNode["name"]=s;
        //2) Save single set of arrhenius params for eig0 to "low-P-rate-constant"
        AnyMap tempNode;
        eig0__[s].getRateParameters(tempNode);
        if (!tempNode.empty()){
            speciesNode["low-P-rate-constant"]=std::move(tempNode);
        }
        tempNode.clear();
        //3) Save list of rate constant params to "rate-constants"
        std::multimap<double, ArrheniusRate> rateMap;
        for (auto iter = ++pressures__[s].begin(); iter->first < 1000; ++iter) {
            for (size_t i = iter->second.first; i < iter->second.second; i++) {
                rateMap.insert({std::exp(iter->first), rates__[s][i]});
            }
        }
        vector<AnyMap> rateList;
        for (const auto& [pressure, rate] : rateMap) {
            tempNode["P"].setQuantity(pressure, "atm");
            rate.getRateParameters(tempNode);
            rateList.push_back(std::move(tempNode));
            tempNode.clear();
        }    
        speciesNode["rate-constants"] = std::move(rateList);
        topLevelList.push_back(std::move(speciesNode));
    }
    rateNode["collider-list"]=std::move(topLevelList);
}
}

