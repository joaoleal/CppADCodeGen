/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */
#include "CppADCGPatternTest.hpp"
#include "../models/cstr.hpp"


namespace CppAD {

    class CppADCGPatternCstrTest : public CppADCGPatternTest {
    public:
        typedef double Base;
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCGD;
    protected:
        static const size_t ns; // number of states in the CSTR model
        static const size_t nm; // number of controls in the CSTR model
        static const size_t npar; // number of parameters in the CSTR model

        static const size_t K; // order of collocation
        static const size_t m; // total number of equations in the collocation model 
        static const size_t na; // number of independent variables of the CSTR model
        static const size_t repeat; // number of time intervals
        std::vector<Base> xx; // default CSTR model values
        std::vector<Base> x; // values for the collocation model
    public:

        inline CppADCGPatternCstrTest(bool verbose = false, bool printValues = false) :
        CppADCGPatternTest(verbose, printValues),
        xx(na) {
            this->epsilonA_ = std::numeric_limits<Base>::epsilon() * 4e3;

            /**
             * CSTR model values
             */
            xx[0] = 0.3; // h
            xx[1] = 7.82e3; // Ca
            xx[2] = 304.65; // Tr
            xx[3] = 301.15; // Tj

            xx[4] = 2.3333e-04; // u1
            xx[5] = 6.6667e-05; // u2

            xx[6] = 6.2e14; // 
            xx[7] = 10080; //
            xx[8] = 2e3; //
            xx[9] = 10e3; //
            xx[10] = 1e-11; //
            xx[11] = 6.6667e-05; //
            xx[12] = 294.15; //
            xx[13] = 294.15; //
            xx[14] = 1000; //
            xx[15] = 4184; //Cp
            xx[16] = -33488; //deltaH
            xx[17] = 299.15; // Tj0
            xx[18] = 302.65; //   Tj2
            xx[19] = 7e5; // cwallj
            xx[20] = 1203; // csteam
            xx[21] = 3.22; //dsteam
            xx[22] = 950.0; //Ug
            xx[23] = 0.48649427192323; //vc6in
            xx[24] = 1000; //rhoj
            xx[25] = 4184; //Cpj
            xx[26] = 0.014; //Vj
            xx[27] = 1e-7; //cwallr

            /**
             * collocation model values
             */
            size_t nvarsk = ns;
            size_t nMstart = npar + nvarsk * K * repeat + nvarsk;

            x.resize(nMstart + repeat * nm, 1.0);
            xNorm_.resize(nMstart + repeat * nm, 1.0);
            // parameters
            for (size_t j = 0; j < npar; j++)
                xNorm_[j] = xx[ns + nm + j];

            size_t s = npar;

            // i = 0 K = 0
            // states
            for (size_t j = 0; j < ns; j++) {
                xNorm_[s++] = xx[j];
            }

            for (size_t i = 0; i < repeat; i++) {
                // controls
                for (size_t j = 0; j < nm; j++) {
                    xNorm_[nMstart + nm * i + j] = xx[ns + j];
                }

                // K = 1
                // states
                for (size_t j = 0; j < ns; j++) {
                    xNorm_[s++] = xx[j];
                }

                // K = 2
                // states
                for (size_t j = 0; j < ns; j++) {
                    xNorm_[s++] = xx[j];
                }

                // K = 3
                // states
                for (size_t j = 0; j < ns; j++) {
                    xNorm_[s++] = xx[j];
                }
            }

#if 0
            x = xNorm_;
            xNorm_.clear();
#else
            for (size_t j = 0; j < xNorm_.size(); j++) {
                xNorm_[j] = xNorm_[j] != 0.0 ? xNorm_[j] : 1.0;
            }

            eqNorm_.resize(repeat * m);
            size_t e = 0;
            for (size_t i = 0; i < repeat; i++) {
                for (size_t k = 0; k < K; k++) {
                    for (size_t j = 0; j < ns; j++) {
                        eqNorm_[e++] = xx[j];
                    }
                }
            }
#endif
        }

        static std::vector<ADCGD> modelCollocation(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<double>*>& atoms) {
            CGAbstractAtomicFun<Base>& atomicCstr = *atoms[0];

            size_t m2 = repeat * m;

            // dependent variable vector 
            std::vector<ADCGD> dep(m2);

            std::vector<ADCGD> dxikdt(ns);
            std::vector<ADCGD> xik(na);

            // parameters
            for (size_t j = 0; j < npar; j++)
                xik[ns + nm + j] = x[j];

            size_t s = npar;
            size_t nvarsk = ns;
            size_t nMstart = npar + nvarsk * K * repeat + nvarsk;
            size_t eq = 0;

            for (size_t i = 0; i < repeat; i++) {
                size_t s0 = s;

                // controls
                for (size_t j = 0; j < nm; j++) {
                    xik[ns + j] = x[nMstart + nm * i + j];
                }

                // K = 1
                for (size_t j = 0; j < ns; j++) {
                    xik[j] = x[s + j]; // states
                }
                s += nvarsk;
                // xik[ns + nm + npar] = x[s + ns];// time

                atomicCstr(xik, dxikdt); // ODE
                for (size_t j = 0; j < ns; j++) {
                    dep[eq + j] = dxikdt[j]
                            + 0.13797958971132715 * x[s0 + j]
                            + -0.10749149571305303 * x[s0 + nvarsk + j]
                            + -0.038928002823013501 * x[s0 + 2 * nvarsk + j]
                            + 0.008439908824739363 * x[s0 + 3 * nvarsk + j];
                }
                eq += ns;

                // K = 2
                for (size_t j = 0; j < ns; j++) {
                    xik[j] = x[s + j]; // states
                }
                s += nvarsk;
                // xik[ns + nm + npar] = x[s + ns];// time

                atomicCstr(xik, dxikdt); // ODE
                for (size_t j = 0; j < ns; j++) {
                    dep[eq + j] = dxikdt[j]
                            + -0.057979589711327127 * x[s0 + j]
                            + 0.11892800282301351 * x[s0 + nvarsk + j]
                            + -0.025841837620280327 * x[s0 + 2 * nvarsk + j]
                            + -0.035106575491406049 * x[s0 + 3 * nvarsk + j];
                }
                eq += ns;

                // K = 3
                for (size_t j = 0; j < ns; j++) {
                    xik[j] = x[s + j]; // states
                }
                s += nvarsk;
                // xik[ns + nm + npar] = x[s + ns];// time

                atomicCstr(xik, dxikdt); // ODE
                for (size_t j = 0; j < ns; j++) {
                    dep[eq + j] = dxikdt[j]
                            + 0.099999999999999978 * x[s0 + j]
                            + -0.18439908824739357 * x[s0 + nvarsk + j]
                            + 0.25106575491406025 * x[s0 + 2 * nvarsk + j]
                            + -0.16666666666666669 * x[s0 + 3 * nvarsk + j];
                }
                eq += ns;

            }

            return dep;
        }

        template<class T>
        static void atomicFunction(const std::vector<AD<T> >& x, std::vector<AD<T> >& y) {
            y = CstrFunc(x);
        }
    };

    const size_t CppADCGPatternCstrTest::ns = 4;
    const size_t CppADCGPatternCstrTest::nm = 2;
    const size_t CppADCGPatternCstrTest::npar = 22;

    const size_t CppADCGPatternCstrTest::K = 3;
    const size_t CppADCGPatternCstrTest::m = K * ns;
    const size_t CppADCGPatternCstrTest::na = ns + nm + npar;

    const size_t CppADCGPatternCstrTest::repeat = 6;
}

using namespace CppAD;

/**
 * @test test the usage of loops for the generation of a orthogonal collocation
 *       model used to integrate a CSTR ODE model
 *       with differential information for all variables (including parameters)
 */
TEST_F(CppADCGPatternCstrTest, AtomicAllVars) {
    using namespace CppAD;

    /**
     * create atomic function for the ODE
     */
    std::vector<AD<double> > ay(ns), ax(na);
    for (size_t j = 0; j < na; j++)
        ax[j] = xx[j];

    checkpoint<double> atomicfun("atomicFunc", atomicFunction<double>, ax, ay);
    std::vector<atomic_base<double>*> atomics(1);
    atomics[0] = &atomicfun;

    testPatternDetectionWithAtomics(modelCollocation, atomics, m, x, repeat);
    testLibCreationWithAtomics("modelCstrAtomicAllVars", modelCollocation, atomics, m, x, repeat);
}

/**
 * @test test the usage of loops for the generation of a orthogonal collocation
 *       model used to integrate a CSTR ODE model,
 *       only with differential information for states and controls
 */
TEST_F(CppADCGPatternCstrTest, Atomic) {

    using namespace CppAD;
    using namespace CppAD::extra;

    std::string modelName = "ctsr";

    // independent variables
    std::vector<ADCGD> ay(ns), ax(na);
    for (size_t j = 0; j < na; j++)
        ax[j] = xx[j];

    CppAD::Independent(ax);

    /**
     * create the CppAD tape as usual
     */
    atomicFunction(ax, ay);

    ADFun<CGD> fun;
    fun.Dependent(ay);

    /**
     * Create the dynamic library model
     */
    CLangCompileModelHelper<double> compHelp1(fun, modelName);
    compHelp1.setCreateForwardZero(true);
    compHelp1.setCreateForwardOne(true);
    compHelp1.setCreateReverseOne(true);
    compHelp1.setCreateReverseTwo(true);


    std::vector<std::set<size_t> > jacSparAll = jacobianSparsitySet<std::vector<std::set<size_t> > >(fun);
    std::vector<std::set<size_t> > jacSpar(jacSparAll.size());
    for (size_t i = 0; i < jacSparAll.size(); i++) {
        // only differential information for states and controls
        std::set<size_t>::const_iterator itEnd = jacSparAll[i].upper_bound(ns + nm - 1);
        if (itEnd != jacSparAll[i].begin())
            jacSpar[i].insert(jacSparAll[i].begin(), itEnd);
    }
    compHelp1.setCustomSparseJacobianElements(jacSpar);

    std::vector<std::set<size_t> > hessSparAll = hessianSparsitySet<std::vector<std::set<size_t> > >(fun);
    std::vector<std::set<size_t> > hessSpar(hessSparAll.size());
    for (size_t i = 0; i < ns + nm; i++) {
        std::set<size_t>::const_iterator it = hessSparAll[i].upper_bound(i); // only the lower left side
        if (it != hessSparAll[i].begin())
            hessSpar[i].insert(hessSparAll[i].begin(), it);
    }
    compHelp1.setCustomSparseHessianElements(hessSpar);

    /**
     * Create the dynamic library
     * (generate and compile source code)
     */
    GccCompiler<double> compiler1;
    std::vector<std::string> flags;
    flags.push_back("-O0");
    flags.push_back("-g");
    flags.push_back("-ggdb");
    flags.push_back("-D_FORTIFY_SOURCE=2");
    compiler1.setCompileFlags(flags);
    compiler1.setSourcesFolder("sources_cstr_atomiclib_" + modelName);

    CLangCompileDynamicHelper<double> compDynHelp(compHelp1);
    compDynHelp.setLibraryName("cstr");
    std::auto_ptr<DynamicLib<Base> > dynamicLib(compDynHelp.createDynamicLibrary(compiler1));
    std::auto_ptr<DynamicLibModel<Base> > modelLib(dynamicLib->model(modelName));

    std::vector<atomic_base<double>*> atomics(1);
    atomics[0] = &modelLib->asAtomic();

    testPatternDetectionWithAtomics(modelCollocation, atomics, m, x, repeat);
    testLibCreationWithAtomics("modelCstrAtomic", modelCollocation, atomics, m, x, repeat);

}