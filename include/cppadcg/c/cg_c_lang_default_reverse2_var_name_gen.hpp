#ifndef CPPAD_CG_C_LANG_DEFAULT_REVERSE2_VAR_NAME_GEN_INCLUDED
#define CPPAD_CG_C_LANG_DEFAULT_REVERSE2_VAR_NAME_GEN_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

namespace CppAD {

    /**
     * Creates variables names for the source code generated for second-order
     * reverse mode calculations.
     * The independent variables are considered to have been registered first,
     * followed by a first level of additional variables and then a second.
     * 
     * @author Joao Leal
     */
    template<class Base>
    class CLangDefaultReverse2VarNameGenerator : public VariableNameGenerator<Base> {
    protected:
        VariableNameGenerator<Base>* _nameGen;
        // the lowest variable ID used for the first independent variable level
        const size_t _minLevel1ID;
        // array name of the independent variables (1st level)
        const std::string _level1Name;
        // the lowest variable ID used for the second independent variable level
        const size_t _minLevel2ID;
        // array name of the independent variables (2nd level)
        const std::string _level2Name;
        // auxiliary string stream
        std::stringstream _ss;
    public:

        CLangDefaultReverse2VarNameGenerator(VariableNameGenerator<Base>* nameGen,
                                             size_t n,
                                             size_t n1) :
            _nameGen(nameGen),
            _minLevel1ID(n + 1),
            _level1Name("tx1"),
            _minLevel2ID(_minLevel1ID + n1),
            _level2Name("py2") {

            CPPADCG_ASSERT_KNOWN(_nameGen != NULL, "The name generator must not be null");

            initialize();
        }

        CLangDefaultReverse2VarNameGenerator(VariableNameGenerator<Base>* nameGen,
                                             size_t n,
                                             const std::string& level1Name,
                                             size_t n1,
                                             const std::string& level2Name) :
            _nameGen(nameGen),
            _minLevel1ID(n + 1),
            _level1Name(level1Name),
            _minLevel2ID(_minLevel1ID + n1),
            _level2Name(level2Name) {

            CPPADCG_ASSERT_KNOWN(_nameGen != NULL, "The name generator must not be null");
            CPPADCG_ASSERT_KNOWN(_level1Name.size() > 0, "The name for the first level must not be empty");
            CPPADCG_ASSERT_KNOWN(_level2Name.size() > 0, "The name for the second level must not be empty");

            initialize();
        }

        virtual const std::vector<FuncArgument>& getDependent() const {
            return _nameGen->getDependent();
        }

        virtual const std::vector<FuncArgument>& getTemporary() const {
            return _nameGen->getTemporary();
        }

        virtual size_t getMinTemporaryVariableID() const {
            return _nameGen->getMinTemporaryVariableID();
        }

        virtual size_t getMaxTemporaryVariableID() const {
            return _nameGen->getMaxTemporaryVariableID();
        }

        virtual std::string generateDependent(const CG<Base>& variable, size_t index) {
            return _nameGen->generateDependent(variable, index);
        }

        virtual std::string generateIndependent(const SourceCodeFragment<Base>& independent) {
            size_t id = independent.variableID();
            if (id < _minLevel1ID) {
                return _nameGen->generateIndependent(independent);
            } else {
                _ss.clear();
                _ss.str("");
                if (id < _minLevel2ID) {
                    _ss << _level1Name << "[" << (id - _minLevel1ID) << "]";
                } else {
                    _ss << _level2Name << "[" << (id - _minLevel2ID) << "]";
                }
                return _ss.str();
            }
        }

        virtual std::string generateTemporary(const SourceCodeFragment<Base>& variable) {
            return _nameGen->generateTemporary(variable);
        }

        virtual void setTemporaryVariableID(size_t minTempID, size_t maxTempID) {
            _nameGen->setTemporaryVariableID(minTempID, maxTempID);
        }

        inline virtual ~CLangDefaultReverse2VarNameGenerator() {
        }

    private:

        inline void initialize() {
            this->_independent = _nameGen->getIndependent(); // copy
            this->_independent.push_back(FuncArgument(_level1Name));
            this->_independent.push_back(FuncArgument(_level2Name));
        }

    };
}

#endif
