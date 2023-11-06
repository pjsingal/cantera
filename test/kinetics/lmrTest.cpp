// #include "gtest/gtest.h"
// #include "cantera/kinetics/Kinetics.h"
// #include "cantera/thermo/ThermoPhase.h"
// #include "cantera/base/Solution.h"
// #include "cantera/base/global.h"

// //SSH key fingerprint:   SHA256:pV/Ash29H8ZN7Al414mZV8cQOoYqPlZei3d21RIWnac pjsingal98@gmail.com

// // COMPILE/EXECUTE COMMANDS
// // scons test-kinetics verbose_tests=y -j4

// namespace Cantera
// {

// class lmrTest : public testing::Test
// {
// public:
//     lmrTest() {}

//     static void SetUpTestCase() {
//         soln_ = newSolution("../data/kineticsfromscratch_LMRtest.yaml"); //import kineticsfromscratch-lmrtest.yaml
//     }

//     static void TearDownTestCase() {
//         soln_.reset();
//     }

// protected:
//     static shared_ptr<Solution> soln_;
// };

// shared_ptr<Solution> lmrTest::soln_;

// TEST_F(lmrTest, reactionCounts) 
// {
//     EXPECT_EQ((size_t) 14, soln_->kinetics()->nReactions());
// }

// TEST_F(lmrTest, PlogLowPressure)
// {

//     string X = "H:1.0, O2:0.0";
//     soln_->thermo()->setState_TPX(900.0, 101325 * 8.0, X);

//     double ktest1=10;//placeholder
//     double ktest2=10;//placeholder


//     // Test that P-log reactions have the right low-pressure limit
//     //set_TP(500.0, 1e-7);
//     //vector<double> kf(7); //7 is number of rxns in that input file. But the solution obj can be used to automate this (the test above chaks that numb rxns is 7)
//     // soln_->kinetics()->getFwdRateConstants(&kf[0]);

//     // Pre-exponential factor decreases by 10^3 for second-order reaction
//     // when converting from cm + mol to m + kmol
//     // double kf0 = k(1.212400e+13, -0.5779, 10872.7);

//     // EXPECT_NEAR(kf0, kf[0], 1e-9 * kf0); //third arg is the relative tolerance (we have floating point numerals which are inexact)
//     EXPECT_DOUBLE_EQ(ktest1, ktest2);
//     //EXPECT_NEAR(kf3, kf[3], 1e-9 * kf3); //kf[3], e.g. selects rxn 3 from input file
// }
// } // namespace Cantera


#include "gtest/gtest.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/Solution.h"
#include "cantera/base/global.h"

namespace Cantera
{

class LmrTest : public testing::Test
{
public:
    LmrTest() {}

    static void SetUpTestCase() {
        soln_ = newSolution("../data/pdep-test.yaml");
    }

    static void TearDownTestCase() {
        soln_.reset();
    }

    void SetUp() override {
        string Xref = "H:1.0, R1A:1.0, R1B:1.0, R2:1.0, R3:1.0, R4:1.0, R5:1.0, R6:1.0";

        soln_->thermo()->setState_TPX(900.0, 101325 * 8.0, Xref);
    }

protected:
    static shared_ptr<Solution> soln_;

    void set_TP(double T, double P) {
        T_ = T;
        RT_ = GasConst_cal_mol_K * T;
        P_ = P;
        soln_->thermo()->setState_TP(T_, P_);
    }

    double k(double A, double n, double Ea) {
        return A * pow(T_, n) * exp(-Ea/RT_);
    }
    double T_, RT_, P_;
};

shared_ptr<Solution> LmrTest::soln_;

TEST_F(LmrTest, reactionCounts)
{
    EXPECT_EQ((size_t) 7, soln_->kinetics()->nReactions());
}

TEST_F(LmrTest, PlogLowPressure)
{
    // Test that P-log reactions have the right low-pressure limit
    set_TP(500.0, 1e-7);
    vector<double> kf(7);
    soln_->kinetics()->getFwdRateConstants(&kf[0]);

    // Pre-exponential factor decreases by 10^3 for second-order reaction
    // when converting from cm + mol to m + kmol
    double kf0 = k(1.212400e+13, -0.5779, 10872.7);
    double kf1 = k(1.230000e+05, 1.53, 4737.0);
    double kf2 = k(2.440000e+7, 1.04, 3980.0);
    double kf3 = k(2.889338e-17*(Avogadro/1e6), 1.98, 4521.0);

    EXPECT_NEAR(kf0, kf[0], 1e-9 * kf0);
    EXPECT_NEAR(kf1, kf[1], 1e-9 * kf1);
    EXPECT_NEAR(kf2, kf[2], 1e-9 * kf2);
    EXPECT_NEAR(kf3, kf[3], 1e-9 * kf3);
}

TEST_F(LmrTest, PlogHighPressure)
{
    // Test that P-log reactions have the right high-pressure limit
    set_TP(500.0, 1e10);
    vector<double> kf(7);
    soln_->kinetics()->getFwdRateConstants(&kf[0]);

    // Pre-exponential factor decreases by 10^3 for second-order reaction
    // when converting from cm + mol to m + kmol
    double kf0 = k(5.963200e+53, -11.529, 52599.6);
    double kf3 = k(2.889338e-17*(Avogadro/1e6), 1.98, 4521.0);

    EXPECT_NEAR(kf0, kf[0], 1e-9 * kf0);
    EXPECT_NEAR(kf3, kf[3], 1e-9 * kf3);
}

TEST_F(LmrTest, PlogDuplicatePressures)
{
    // Test that multiple rate expressions are combined when necessary
    set_TP(500.0, 1e10);
    vector<double> kf(7);

    soln_->kinetics()->getFwdRateConstants(&kf[0]);
    double kf1 = k(1.3700e+14, -0.79, 17603.0) + k(1.2800e+03, 1.71, 9774.0);
    double kf2 = k(-7.4100e+27, -5.54, 12108.0) + k(1.9000e+12, -0.29, 8306.0);

    EXPECT_NEAR(kf1, kf[1], 1e-9 * kf1);
    EXPECT_NEAR(kf2, kf[2], 1e-9 * kf2);
}

TEST_F(LmrTest, PlogCornerCases)
{
    // Test rate evaluation at the corner cases where the pressure
    // is exactly of the specified interpolation values
    set_TP(500.0, 101325);
    vector<double> kf(7);
    soln_->kinetics()->getFwdRateConstants(&kf[0]);

    double kf0 = k(4.910800e+28, -4.8507, 24772.8);
    double kf1 = k(1.2600e+17, -1.83, 15003.0) + k(1.2300e+01, 2.68, 6335.0);
    double kf2 = k(3.4600e+9, 0.442, 5463.0);

    EXPECT_NEAR(kf0, kf[0], 1e-9 * kf0);
    EXPECT_NEAR(kf1, kf[1], 1e-9 * kf1);
    EXPECT_NEAR(kf2, kf[2], 1e-9 * kf2);
}

TEST_F(LmrTest, PlogIntermediatePressure1)
{
    set_TP(1100.0, 20*101325);
    vector<double> ropf(7);
    soln_->kinetics()->getFwdRatesOfProgress(&ropf[0]);

    // Expected rates computed using Chemkin
    // ROP increases by 10**3 when converting from mol/cm3 to kmol/m3
    EXPECT_NEAR(3.100682e+05, ropf[0], 1e2);
    EXPECT_NEAR(2.006871e+05, ropf[1], 1e2);
    EXPECT_NEAR(4.468658e+06, ropf[2], 1e2);
    EXPECT_NEAR(1.774796e+06, ropf[3], 1e2);
}

TEST_F(LmrTest, PlogIntermediatePressure2)
{
    set_TP(1100.0, 0.5*101325);
    vector<double> ropf(7);
    soln_->kinetics()->getFwdRatesOfProgress(&ropf[0]);

    EXPECT_NEAR(5.244649e+02, ropf[0], 5e-2);
    EXPECT_NEAR(2.252537e+02, ropf[1], 2e-2);
    EXPECT_NEAR(2.985338e+03, ropf[2], 3e-1);
    EXPECT_NEAR(1.109248e+03, ropf[3], 1e-1);
}

TEST_F(LmrTest, PlogIntermediatePressure3)
{
    set_TP(800.0, 70*101325);
    vector<double> ropf(7);
    soln_->kinetics()->getFwdRatesOfProgress(&ropf[0]);

    EXPECT_NEAR(2.274501e+04, ropf[0], 1e+1);
    EXPECT_NEAR(2.307191e+05, ropf[1], 1e+2);
    EXPECT_NEAR(2.224601e+07, ropf[2], 1e+3);
    EXPECT_NEAR(1.007440e+07, ropf[3], 1e+3);
}

TEST_F(LmrTest, ChebyshevIntermediate1)
{
    // Test Chebyshev rates in the normal interpolation region
    vector<double> kf(7);

    set_TP(1100.0, 20 * 101325);
    soln_->kinetics()->getFwdRateConstants(&kf[0]);
    // Expected rates computed using RMG-py
    EXPECT_NEAR(3.130698657e+06, kf[4], 1e-1);
    EXPECT_NEAR(1.187949573e+00, kf[5], 1e-7);

    // Rate for a reaction specified as "molec" instead of "mol" should
    // be higher by a factor of the Avogadro constant (in mol, not kmol).
    // Accuracy is limited by the low precision used by ck2cti
    EXPECT_NEAR(kf[4], kf[6]/(Avogadro*1e-3), 5e2);
}

TEST_F(LmrTest, ChebyshevIntermediate2)
{
    // Test Chebyshev rates in the normal interpolation region
    vector<double> kf(7);

    set_TP(400.0, 0.1 * 101325);
    soln_->kinetics()->getFwdRateConstants(&kf[0]);
    // Expected rates computed using RMG-py
    EXPECT_NEAR(1.713599902e+05, kf[4], 1e-3);
    EXPECT_NEAR(9.581780687e-24, kf[5], 1e-31);
    EXPECT_NEAR(kf[4], kf[6]/(Avogadro*1e-3), 1e2);
}

TEST_F(LmrTest, ChebyshevIntermediateROP)
{
    set_TP(1100.0, 30 * 101325);
    vector<double> ropf(7);
    // Expected rates computed using Chemkin
    soln_->kinetics()->getFwdRatesOfProgress(&ropf[0]);
    EXPECT_NEAR(4.552930e+03, ropf[4], 1e-1);
    EXPECT_NEAR(4.877390e-02, ropf[5], 1e-5);
}

TEST_F(LmrTest, ChebyshevEdgeCases)
{
    vector<double> kf(7);

    // Minimum P
    set_TP(500.0, 1000.0);
    soln_->kinetics()->getFwdRateConstants(&kf[0]);
    EXPECT_NEAR(1.225785655e+06, kf[4], 1e-2);

    // Maximum P
    set_TP(500.0, 1.0e7);
    soln_->kinetics()->getFwdRateConstants(&kf[0]);
    EXPECT_NEAR(1.580981157e+03, kf[4], 1e-5);

    // Minimum T
    set_TP(300.0, 101325);
    soln_->kinetics()->getFwdRateConstants(&kf[0]);
    EXPECT_NEAR(5.405987017e+03, kf[4], 1e-5);

    // Maximum T
    set_TP(2000.0, 101325);
    soln_->kinetics()->getFwdRateConstants(&kf[0]);
    EXPECT_NEAR(3.354054351e+07, kf[4], 1e-1);
}

} // namespace Cantera