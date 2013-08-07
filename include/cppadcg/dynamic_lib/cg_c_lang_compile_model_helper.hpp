#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_INCLUDED
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
     * Useful class for generating source code for the creation of a dynamic
     * library.
     * 
     * @author Joao Leal
     */
    template<class Base>
    class CLangCompileModelHelper {
        typedef CppAD::CG<Base> CGBase;
        typedef CppAD::AD<CGBase> ADCG;
        typedef CppAD::vector<std::set<size_t> > SparsitySetType;
    public:
        static const std::string FUNCTION_FORWAD_ZERO;
        static const std::string FUNCTION_JACOBIAN;
        static const std::string FUNCTION_HESSIAN;
        static const std::string FUNCTION_FORWARD_ONE;
        static const std::string FUNCTION_REVERSE_ONE;
        static const std::string FUNCTION_REVERSE_TWO;
        static const std::string FUNCTION_SPARSE_JACOBIAN;
        static const std::string FUNCTION_SPARSE_HESSIAN;
        static const std::string FUNCTION_JACOBIAN_SPARSITY;
        static const std::string FUNCTION_HESSIAN_SPARSITY;
        static const std::string FUNCTION_HESSIAN_SPARSITY2;
        static const std::string FUNCTION_SPARSE_FORWARD_ONE;
        static const std::string FUNCTION_SPARSE_REVERSE_ONE;
        static const std::string FUNCTION_SPARSE_REVERSE_TWO;
        static const std::string FUNCTION_FORWARD_ONE_SPARSITY;
        static const std::string FUNCTION_REVERSE_ONE_SPARSITY;
        static const std::string FUNCTION_REVERSE_TWO_SPARSITY;
        static const std::string FUNCTION_INFO;
        static const std::string FUNCTION_ATOMIC_FUNC_NAMES;
    protected:
        static const std::string CONST;

        /**
         * Useful class for storing matrix indexes
         */
        class Position {
        public:
            bool defined;
            std::vector<size_t> row;
            std::vector<size_t> col;

            inline Position() :
                defined(false) {
            }

            inline Position(const std::vector<size_t>& r, const std::vector<size_t>& c) :
                defined(true),
                row(r),
                col(c) {
                CPPADCG_ASSERT_KNOWN(r.size() == c.size(), "The number of row indexes must be the same as the number of column indexes.");
            }

            template<class VectorSet>
            inline Position(const VectorSet& elements) :
                defined(true) {
                size_t nnz = 0;
                for (size_t i = 0; i < elements.size(); i++) {
                    nnz += elements[i].size();
                }
                row.resize(nnz);
                col.resize(nnz);

                nnz = 0;
                std::set<size_t>::const_iterator it;
                for (size_t i = 0; i < elements.size(); i++) {
                    for (it = elements[i].begin(); it != elements[i].end(); ++it) {
                        row[nnz] = i;
                        col[nnz] = *it;
                        nnz++;
                    }
                }
            }
        };

        /**
         * 
         */
        class LocalSparsityInfo {
        public:
            SparsitySetType sparsity;
            std::vector<size_t> rows;
            std::vector<size_t> cols;
        };

    protected:
        /**
         * the  model
         */
        ADFun<CGBase>& _funOrig;
        /**
         * Altered model with loops
         */
        ADFun<CGBase>* _funLoops;
        /**
         * work  model
         */
        ADFun<CGBase>* _fun;
        /**
         * the name of the model
         */
        std::string _name;
        /**
         * the name of the data type used in operations
         */
        const std::string _baseTypeName;
        /**
         * Typical values of the independent vector 
         */
        std::vector<Base> _x;
        bool _zero;
        bool _jacobian;
        bool _hessian;
        bool _sparseJacobian;
        bool _sparseHessian;
        bool _hessianByEquation;
        bool _forwardOne;
        bool _reverseOne;
        bool _reverseTwo;
        JacobianADMode _jacMode;
        /**
         * Custom Jacobian element indexes 
         */
        Position _custom_jac;
        LocalSparsityInfo _jacSparsity;
        /**
         * Custom Hessian element indexes 
         */
        Position _custom_hess;
        LocalSparsityInfo _hessSparsity;
        std::vector<LocalSparsityInfo> _hessSparsities;
        /**
         * The order of the atomic functions
         */
        std::vector<std::string> _atomicFunctions;
        /**
         * A string cache for code generation
         */
        std::ostringstream _cache;
        /**
         * maximum number of assignments per function (~ lines)
         */
        size_t _maxAssignPerFunc;
        /**
         * Whether or not to print some progress information to the standard 
         * output
         */
        bool _verbose;
        /**
         * auxiliary variable to measure the elapsed time
         */
        double _beginTime;
        /**
         * 
         */
        std::vector<std::set<size_t> > _relatedDepCandidates;
        /**
         * 
         */
        std::set<LoopAtomicFun<Base>* > _loopAtomics;
    public:

        /**
         * Creates a new C language compilation helper for a model
         * 
         * @param fun The ADFun with the taped model (should only be deleted
         *            after this object)
         * @param model The model name (must be a valid C function name)
         */
        CLangCompileModelHelper(ADFun<CppAD::CG<Base> >& fun, const std::string& model) :
            _funOrig(fun),
            _funLoops(NULL),
            _fun(&fun),
            _name(model),
            _baseTypeName(CLangCompileModelHelper<Base>::baseTypeName()),
            _zero(true),
            _jacobian(false),
            _hessian(false),
            _sparseJacobian(false),
            _sparseHessian(false),
            _hessianByEquation(false),
            _forwardOne(false),
            _reverseOne(false),
            _reverseTwo(false),
            _jacMode(AUTOMATIC),
            _maxAssignPerFunc(20000),
            _beginTime(0) {

            CPPADCG_ASSERT_KNOWN(!_name.empty(), "Model name cannot be empty");
            CPPADCG_ASSERT_KNOWN((_name[0] >= 'a' && _name[0] <= 'z') ||
                                 (_name[0] >= 'A' && _name[0] <= 'Z'),
                                 "Invalid model name character");
            for (size_t i = 1; i < _name.size(); i++) {
                char c = _name[i];
                CPPADCG_ASSERT_KNOWN((c >= 'a' && c <= 'z') ||
                                     (c >= 'A' && c <= 'Z') ||
                                     (c >= '0' && c <= '9') ||
                                     c == '_'
                                     , "Invalid model name character");
            }
        }

        /**
         * Provides the model name which should be a valid C function name.
         * 
         * @return the model name 
         */
        inline const std::string& getName() const {
            return _name;
        }

        /**
         * Defines typical values for the independent variable vector. These 
         * values can be usefull when there is a need to call atomic functions,
         * since they may allow to reduce some operations.
         * 
         * @param x The typical values. An empty vector removes the currently
         *          defined values.
         */
        template<class VectorBase>
        inline void setTypicalIndependentValues(const VectorBase& x) {
            CPPAD_ASSERT_KNOWN(x.size() == 0 || x.size() == _fun->Domain(),
                               "Invalid independent variable vector size");
            _x.resize(x.size());
            for (size_t i = 0; i < x.size(); i++) {
                _x[i] = x[i];
            }
        }

        inline void setRelatedDependents(const std::vector<std::set<size_t> >& relatedDepCandidates) {
            _relatedDepCandidates = relatedDepCandidates;
        }

        inline const std::vector<std::set<size_t> >& getRelatedDependents() const {
            return _relatedDepCandidates;
        }

        /**
         * Determines whether or not to generate source-code for a function
         * that evaluates a dense Hessian.
         * 
         * @return true if source-code for a dense Hessian should be created,
         *         false otherwise
         */
        inline bool isCreateHessian() const {
            return _hessian;
        }

        /**
         * Defines whether or not to generate source-code for a function
         * that evaluates a dense Hessian.
         * 
         * @param create true if source-code for a dense Hessian should be
         *               created, false otherwise
         */
        inline void setCreateHessian(bool create) {
            _hessian = create;
        }

        /**
         * Provides the Automatic Differentiation mode used to generate the
         * source code for the Jacobian
         * 
         * @return the Automatic Differentiation mode
         */
        inline JacobianADMode getJacobianADMode() const {
            return _jacMode;
        }

        /**
         * Defines the Automatic Differentiation mode used to generate the
         * source code for the Jacobian
         * 
         * @param mode the Automatic Differentiation mode
         */
        inline void setJacobianADMode(JacobianADMode mode) {
            _jacMode = mode;
        }

        /**
         * Determines whether or not to generate source-code for a function
         * that evaluates a dense Jacobian.
         * 
         * @return true if source-code for a dense Jacobian should be created,
         *         false otherwise
         */
        inline bool isCreateJacobian() const {
            return _jacobian;
        }

        /**
         * Defines whether or not to generate source-code for a function
         * that evaluates a dense Jacobian.
         * 
         * @param create true if source-code for a dense Jacobian should be
         *               created, false otherwise
         */
        inline void setCreateJacobian(bool create) {
            _jacobian = create;
        }

        /**
         * Determines whether or not to generate source-code for a function
         * that evaluates a sparse Hessian. If ReverseTwo is also enabled the
         * generated source-code will use the individual generated functions
         * from the second-order reverse mode. 
         * Enabling the generation of individuals functions for reverse-mode
         * can have a negative impact on the performe of the evaluation of the
         * sparse hessian since Hessian symmetry will not be exploited. To
         * improve performance one can request only the upper or lower elements
         * of the hessian using setCustomSparseHessianElements().
         * 
         * @see setCustomSparseHessianElements()
         * 
         * @return true if source-code for a sparse Hessian should be created,
         *         false otherwise
         */
        inline bool isCreateSparseHessian() const {
            return _sparseHessian;
        }

        /**
         * Defines whether or not to generate source-code for a function
         * that evaluates a sparse Hessian. If ReverseTwo is also enabled the
         * generated source-code will use the individual generated functions
         * from the second-order reverse mode. 
         * Enabling the generation of individuals functions for reverse-mode
         * can have a negative impact on the performe of the evaluation of the
         * sparse hessian since Hessian symmetry will not be exploited. To
         * improve performance one can request only the upper or lower elements
         * of the hessian using setCustomSparseHessianElements().
         * 
         * @see setCustomSparseHessianElements()
         * 
         * @param create true if source-code for a sparse Hessian should be
         *               created, false otherwise
         */
        inline void setCreateSparseHessian(bool create) {
            _sparseHessian = create;
        }

        /**
         * Determines whether or not to generate source-code for a function that 
         * provides the Hessian sparsity pattern for each equation/dependent,
         * when the sparse hessian creation is enabled.
         * Even if this flag is set to false the function can still be generated
         * if the second-order reverse mode is enabled.
         * 
         * @return true if source-code for a Hessians sparsities patterns should
         *         be created, false otherwise
         */
        inline bool isCreateHessianSparsityByEquation() const {
            return _hessianByEquation;
        }

        /**
         * Defines whether or not to generate source-code for a function that 
         * provides the Hessian sparsity pattern for each equation/dependent,
         * when the sparse hessian creation is enabled.
         * Even if this flag is set to false the function can still be generated
         * if the second-order reverse mode is enabled.
         * 
         * @param create true if source-code for a Hessians sparsities should be
         *               created, false otherwise
         */
        inline void setCreateHessianSparsityByEquation(bool create) {
            _hessianByEquation = create;
        }

        /**
         * Determines whether or not to generate source-code for a function
         * that evaluates a sparse Jacobian. If ReverseOne or ForwardOne 
         * functions are enabled, then the sparse Jacobian evaluation might
         * use those functions.
         * Enabling the generation of individuals functions for reverse-mode
         * can have a small negative impact on the performe of the evaluation of
         * the parse Jacobian.

         * @see setCustomSparseJacobianElements()
         * 
         * @return true if source-code for a sparse Jacobian should be created,
         *         false otherwise
         */
        inline bool isCreateSparseJacobian() const {
            return _sparseJacobian;
        }

        /**
         * Defines whether or not to generate source-code for a function
         * that evaluates a sparse Jacobian. If ReverseOne or ForwardOne 
         * functions are enabled, then the sparse Jacobian evaluation might
         * use those functions.
         * Enabling the generation of individuals functions for reverse-mode
         * can have a small negative impact on the performe of the evaluation of
         * the parse Jacobian.

         * @see setCustomSparseJacobianElements()
         * 
         * @param create true if source-code for a sparse Jacobian should be
         *               created, false otherwise
         */
        inline void setCreateSparseJacobian(bool create) {
            _sparseJacobian = create;
        }

        /**
         * Determines whether or not to generate source-code for a function
         * that evaluates the original model.
         * 
         * @return true if source-code for the original model should be created,
         *         false otherwise
         */
        inline bool isCreateForwardZero() const {
            return _zero;
        }

        /**
         * Defines whether or not to generate source-code for a function
         * that evaluates the original model.
         * 
         * @return create true if source-code for the original model should be
         *                created, false otherwise
         */
        inline void setCreateForwardZero(bool create) {
            _zero = create;
        }

        /**
         * Determines whether or not to generate source-code for the
         * first-order forward mode that is used for the evaluation of the
         * Jacobian when the model is used through a user defined atomic
         * AD function.
         * Enabling the generation of individuals functions for forward-mode
         * might have a small negative impact on the performe of the evaluation
         * of the sparse Jacobian (if forward mode is selected).
         * 
         * @see isCreateSparseJacobian()
         * 
         * @return true if the generation of the source for first-order forward
         *         mode is enabled, false otherwise.
         */
        inline bool isCreateSparseForwardOne() const {
            return _forwardOne;
        }

        /**
         * Defines whether or not to generate source-code for the
         * first-order forward mode that is used for the evaluation of the
         * Jacobian when the model is used through a user defined atomic
         * AD function.
         * Enabling the generation of individuals functions for forward-mode
         * might have a small negative impact on the performe of the evaluation
         * of the sparse Jacobian (if forward-mode is selected).
         * 
         * @see setCreateSparseJacobian()
         * 
         * @param create true if the generation of the source for first-order 
         *               forward mode is enabled, false otherwise.
         */
        inline void setCreateForwardOne(bool create) {
            _forwardOne = create;
        }

        /**
         * Determines whether or not to generate source-code for the
         * first-order reverse mode that is used for the evaluation of the
         * Jacobian when the model is used through a user defined atomic
         * AD function.
         * Enabling the generation of individuals functions for reverse-mode
         * might have a small negative impact on the performe of the evaluation
         * of the sparse Jacobian (if reverse-mode is selected).
         * 
         * @see isCreateSparseJacobian()
         * 
         * @return true if the generation of the source for first-order reverse
         *         mode is enabled, false otherwise.
         */
        inline bool isCreateReverseOne() const {
            return _reverseOne;
        }

        /**
         * Determines whether or not to generate source-code for the
         * first-order reverse mode that is used for the evaluation of the
         * Jacobian when the model is used through a user defined atomic
         * AD function.
         * Enabling the generation of individuals functions for reverse-mode
         * might have a small negative impact on the performe of the evaluation
         * of the sparse Jacobian (if reverse-mode is selected).
         * 
         * @see setCreateSparseJacobian()
         * 
         * @return true if the generation of the source for first-order reverse
         *         mode is enabled, false otherwise.
         */
        inline void setCreateReverseOne(bool create) {
            _reverseOne = create;
        }

        /**
         * Determines whether or not to generate source-code for the
         * second-order reverse mode that is used for the evaluation of the
         * hessian when the model is used through a user defined atomic
         * AD function.
         * Enabling the generation of individuals functions for reverse-mode
         * can have a negative impact on the performe of the evaluation of the
         * sparse hessian.
         * 
         * Warning: only the values for px[j * (k+1)] will be defined, since
         *          px[j * (k+1) + 1] is not used during the hessian evaluation.
         * 
         * @return true if the generation of the source for second-order reverse
         *         mode is enabled, false otherwise.
         */
        inline bool isCreateReverseTwo() const {
            return _reverseOne;
        }

        /**
         * Defines whether or not to enable the generation of the source-code 
         * for the second-order reverse mode that is used for the evaluation
         * of the hessian when the model is used through a user defined atomic
         * AD function.
         * Enabling the generation of individuals functions for reverse-mode
         * can have a negative impact on the performe of the evaluation of the
         * hessian. To improve performance one can request only the upper or 
         * lower elements of the hessian using setCustomSparseHessianElements()
         * and later only request those elements through the outer model (ADFun).
         * 
         * Warning: only the values for px[j * (k+1)] will be defined, since
         *          px[j * (k+1) + 1] is not used during the hessian evaluation.
         * 
         * @param create true to enable the generation of the source for
         *               second-order reverse mode, false otherwise.
         */
        inline void setCreateReverseTwo(bool create) {
            _reverseTwo = create;
        }

        inline void setCustomSparseJacobianElements(const std::vector<size_t>& row,
                                                    const std::vector<size_t>& col) {
            _custom_jac = Position(row, col);
        }

        template<class VectorSet>
        inline void setCustomSparseJacobianElements(const VectorSet& elements) {
            _custom_jac = Position(elements);
        }

        inline void setCustomSparseHessianElements(const std::vector<size_t>& row,
                                                   const std::vector<size_t>& col) {
            _custom_hess = Position(row, col);
        }

        template<class VectorSet>
        inline void setCustomSparseHessianElements(const VectorSet& elements) {
            _custom_hess = Position(elements);
        }

        inline size_t getMaxAssignmentsPerFunc() const {
            return _maxAssignPerFunc;
        }

        inline void setMaxAssignmentsPerFunc(size_t maxAssignPerFunc) {
            _maxAssignPerFunc = maxAssignPerFunc;
        }

        inline virtual ~CLangCompileModelHelper() {
            delete _funLoops;

            typename std::set<LoopAtomicFun<Base>* >::const_iterator it;
            for (it = _loopAtomics.begin(); it != _loopAtomics.end(); ++it) {
                delete *it;
            }
        }

        static inline std::string baseTypeName();

    protected:

        inline bool isVerbose() const {
            return _verbose;
        }

        inline void setVerbose(bool verbose) {
            _verbose = verbose;
        }

        virtual VariableNameGenerator<Base>* createVariableNameGenerator(const std::string& depName,
                                                                         const std::string& indepName,
                                                                         const std::string& tmpName,
                                                                         const std::string& tmpArrayName);

        virtual void compileSources(CLangCompiler<Base>& compiler, bool posIndepCode);

        virtual void generateLoops();

        virtual void generateInfoSource(std::map<std::string, std::string>& sources);

        virtual void generateAtomicFuncNames(std::map<std::string, std::string>& sources);

        virtual void generateZeroSource(std::map<std::string, std::string>& sources);


        /***********************************************************************
         * Jacobian
         **********************************************************************/

        virtual void generateJacobianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseJacobianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseJacobianSource(std::map<std::string, std::string>& sources,
                                                  bool forward);

        /**
         * Loops
         */
        virtual void prepareSparseJacobianWithLoops(CodeHandler<Base>& handler,
                                                    std::vector<CGBase>& jac);

        virtual OperationNode<Base>* findSparseJacLoopResult(CodeHandler<Base>& handler,
                                                             size_t e,
                                                             size_t nnz,
                                                             OperationNode<Base>* jacNode,
                                                             std::map<LoopAtomicFun<Base>*, std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > >& jacIndexPatterns,
                                                             std::vector<JacOrigElementLoopInfo<Base> >& garbageCollection);

        virtual OperationNode<Base>* handleSparseJacLoopResult(CodeHandler<Base>& handler,
                                                               size_t e,
                                                               size_t nnz,
                                                               OperationNode<Base>* jacNode,
                                                               std::map<LoopAtomicFun<Base>*, std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > >& jacIndexPatterns,
                                                               Argument<Base> arg,
                                                               std::vector<JacOrigElementLoopInfo<Base> >& garbageCollection);

        virtual void generateSparseJacobianForRevSource(std::map<std::string, std::string>& sources,
                                                        bool forward);

        /***********************************************************************
         * Hessian
         **********************************************************************/

        virtual void generateHessianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseHessianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseHessianSourceDirectly(std::map<std::string, std::string>& sources);

        virtual void generateSparseHessianSourceFromRev2(std::map<std::string, std::string>& sources);

        virtual void determineSecondOrderElements4Eval(std::vector<size_t>& userRows,
                                                       std::vector<size_t>& userCols);

        /***********************************************************************
         * Sparsities for forward/reverse
         **********************************************************************/

        virtual void generateSparsity1DSource(const std::string& function,
                                              const std::vector<size_t>& sparsity);

        virtual void generateSparsity2DSource(const std::string& function,
                                              const LocalSparsityInfo& sparsity);

        virtual void generateSparsity2DSource2(const std::string& function,
                                               const std::vector<LocalSparsityInfo>& sparsities);

        virtual void generateSparsity1DSource2(const std::string& function,
                                               const std::map<size_t, std::vector<size_t> >& rows);

        /***********************************************************************
         * Forward 1 mode
         **********************************************************************/

        virtual void generateSparseForwardOneSources(std::map<std::string, std::string>& sources);

        virtual void generateForwardOneSources(std::map<std::string, std::string>& sources);

        /***********************************************************************
         * Reverse 1 mode
         **********************************************************************/

        virtual void generateSparseReverseOneSources(std::map<std::string, std::string>& sources);

        virtual void generateReverseOneSources(std::map<std::string, std::string>& sources);

        /***********************************************************************
         * Reverse 2 mode
         **********************************************************************/

        virtual void generateSparseReverseTwoSources(std::map<std::string, std::string>& sources);

        virtual void generateReverseTwoSources(std::map<std::string, std::string>& sources);

        virtual void generateGlobalDirectionalFunctionSource(const std::string& function,
                                                             const std::string& function2_suffix,
                                                             const std::string& function_sparsity,
                                                             const std::map<size_t, std::vector<size_t> >& elements,
                                                             std::map<std::string, std::string>& sources);

        virtual void generateFunctionDeclarationSource(std::ostringstream& cache,
                                                       const std::string& model_function,
                                                       const std::string& suffix,
                                                       const std::map<size_t, std::vector<size_t> >& elements,
                                                       const std::string& argsDcl);

        /***********************************************************************
         * Sparsities
         **********************************************************************/

        virtual void determineJacobianSparsity();

        virtual void generateJacobianSparsitySource(std::map<std::string, std::string>& sources);

        virtual void determineHessianSparsity();

        virtual void generateHessianSparsitySource(std::map<std::string, std::string>& sources);


        static inline std::map<size_t, std::vector<std::set<size_t> > > determineOrderByCol(const std::map<size_t, std::vector<size_t> >& elements,
                                                                                            const LocalSparsityInfo& sparsity);

        static inline std::map<size_t, std::vector<std::set<size_t> > > determineOrderByCol(const std::map<size_t, std::vector<size_t> >& elements,
                                                                                            const std::vector<size_t>& userRows,
                                                                                            const std::vector<size_t>& userCols);

        static inline std::map<size_t, std::vector<std::set<size_t> > > determineOrderByRow(const std::map<size_t, std::vector<size_t> >& elements,
                                                                                            const LocalSparsityInfo& sparsity);

        static inline std::map<size_t, std::vector<std::set<size_t> > > determineOrderByRow(const std::map<size_t, std::vector<size_t> >& elements,
                                                                                            const std::vector<size_t>& userRows,
                                                                                            const std::vector<size_t>& userCols);

    private:
        void inline startingGraphCreation(const std::string& jobName);

        void inline finishedGraphCreation();

        CLangCompileModelHelper(const CLangCompileModelHelper&); // not implemented

        CLangCompileModelHelper& operator=(const CLangCompileModelHelper&); // not implemented

        friend class
        CLangCompileDynamicHelper<Base>;
    };

}

#endif
