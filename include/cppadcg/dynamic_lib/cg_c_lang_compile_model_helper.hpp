#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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
        typedef std::pair<size_t, size_t> TapeVarType; // tape independent -> reference orig independent (temporaries only)
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
         * Saves sparsity information in more than one format
         */
        class LocalSparsityInfo {
        public:
            /**
             * Calculated sparsity from the model
             * (may differ from the requested sparsity)
             */
            SparsitySetType sparsity;
            // rows (in a custom order)
            std::vector<size_t> rows;
            // columns (in a custom order)
            std::vector<size_t> cols;
        };

        /**
         * Used for coloring
         */
        class Color {
        public:
            /// all row with this color
            std::set<size_t> rows;
            /// maps column indexes to the corresponding row
            std::map<size_t, size_t> column2Row;
            /// maps row indexes to the corresponding columns
            std::map<size_t, std::set<size_t> > row2Columns;
            /// used columns
            std::set<size_t> forbiddenRows;
        };

    protected:
        /**
         * the orignal  model
         */
        ADFun<CGBase>& _fun;
        /**
         * Altered model without the loop equations and with extra dependents
         * for the non-indexed temporary variables used by loops
         */
        LoopFreeModel<Base>* _funNoLoops;
        /**
         * loop models
         */
        std::set<LoopModel<Base>*> _loopTapes;
        /**
         * the name of the model
         */
        std::string _name;
        /**
         * the name of the data type used in operations
         */
        const std::string _baseTypeName;
        /**
         * the maximum precision used to print values
         */
        size_t _parameterPrecision;
        /**
         * Typical values of the independent vector 
         */
        std::vector<Base> _x;
        /// generate source code for the zero order model evaluation
        bool _zero;
        bool _zeroEvaluated;
        /// generate source code for a dense Jacobian
        bool _jacobian;
        /// generate source code for a dense Hessian
        bool _hessian;
        /// generate source code for a sparse Jacobian
        bool _sparseJacobian;
        /// generate source code for a sparse Hessian
        bool _sparseHessian;
        /**
         * generate source-code the Hessian sparsity pattern for each 
         * equation/dependent
         */
        bool _hessianByEquation;
        /// generate source code for forward first order mode
        bool _forwardOne;
        /// generate source code for reverse first order mode
        bool _reverseOne;
        /// generate source code for reverse second order mode
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
        /**
         * hessian sparsity from the model for each equation
         */
        std::vector<LocalSparsityInfo> _hessSparsities;
        /**
         * The order of the atomic functions
         */
        std::vector<std::string> _atomicFunctions;
        /**
         * Maps each atomic function ID to the independent variable indexes 
         * which affect a call of that atomic function
         */
        std::map<size_t, std::set<size_t> >* _atomicsIndeps;
        /**
         * A string cache for code generation
         */
        std::ostringstream _cache;
        /**
         * maximum number of assignments per function (~ lines)
         */
        size_t _maxAssignPerFunc;
        /**
         * 
         */
        std::vector<std::set<size_t> > _relatedDepCandidates;
        /**
         * Maps the column groups of each loop model to the set of columns
         * (loop->group->{columns->{compressed forward 1 position} })
         */
        std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > > _loopFor1Groups;
        /**
         * Jacobian columns with a contribution from non loop equations
         *  [var]{compressed forward 1 position}
         */
        std::map<size_t, std::set<size_t> > _nonLoopFor1Elements;
        /**
         * Maps the column groups of each loop model to the set of columns
         * (loop->group->{row->{compressed reverse 1 position} })
         */
        std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > > _loopRev1Groups;
        /**
         * Jacobian rows with a contribution from non loop equations
         *  [eq]{compressed reverse 1 position}
         */
        std::map<size_t, std::set<size_t> > _nonLoopRev1Elements;
        /**
         * Maps the row groups of each loop model to the set of rows
         * (loop->group->{rows->{compressed reverse 2 position} })
         */
        std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > > _loopRev2Groups;
        /**
         * Hessian rows with a contribution from non loop equations
         *  [var]{compressed reverse 2 position}
         */
        std::map<size_t, std::set<size_t> > _nonLoopRev2Elements;
        /**
         * 
         */
        JobTimer* _jobTimer;
    public:

        /**
         * Creates a new C language compilation helper for a model
         * 
         * @param fun The ADFun with the taped model (should only be deleted
         *            after this object)
         * @param model The model name (must be a valid C function name)
         */
        CLangCompileModelHelper(ADFun<CppAD::CG<Base> >& fun, const std::string& model) :
            _fun(fun),
            _funNoLoops(NULL),
            _name(model),
            _baseTypeName(CLangCompileModelHelper<Base>::baseTypeName()),
            _parameterPrecision(std::numeric_limits<Base>::digits10),
            _zero(true),
            _zeroEvaluated(false),
            _jacobian(false),
            _hessian(false),
            _sparseJacobian(false),
            _sparseHessian(false),
            _hessianByEquation(false),
            _forwardOne(false),
            _reverseOne(false),
            _reverseTwo(false),
            _jacMode(AUTOMATIC),
            _atomicsIndeps(NULL),
            _maxAssignPerFunc(20000),
            _jobTimer(NULL) {

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
            CPPAD_ASSERT_KNOWN(x.size() == 0 || x.size() == _fun.Domain(),
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
         * Provides the maximum precision used to print constant values in the
         * generated source code
         * 
         * @return the maximum number of digits
         */
        virtual size_t getParameterPrecision() const {
            return _parameterPrecision;
        }

        /**
         * Defines the maximum precision used to print constant values in the
         * generated source code
         * 
         * @param p the maximum number of digits
         */
        virtual void setParameterPrecision(size_t p) {
            _parameterPrecision = p;
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
            delete _funNoLoops;
            delete _atomicsIndeps;

            typename std::set<LoopModel<Base>* >::const_iterator it;
            for (it = _loopTapes.begin(); it != _loopTapes.end(); ++it) {
                delete *it;
            }
        }

    public:
        static inline std::string baseTypeName();

        template<class T>
        static void generateFunctionDeclarationSource(std::ostringstream& cache,
                                                      const std::string& model_function,
                                                      const std::string& suffix,
                                                      const std::map<size_t, T>& elements,
                                                      const std::string& argsDcl);

    protected:

        virtual VariableNameGenerator<Base>* createVariableNameGenerator(const std::string& depName,
                                                                         const std::string& indepName,
                                                                         const std::string& tmpName,
                                                                         const std::string& tmpArrayName);

        virtual void compileSources(CLangCompiler<Base>& compiler,
                                    bool posIndepCode,
                                    JobTimer* timer = NULL) throw (CGException);

        virtual void generateLoops() throw (CGException);

        virtual void generateInfoSource(std::map<std::string, std::string>& sources);

        virtual void generateAtomicFuncNames(std::map<std::string, std::string>& sources);

        virtual bool isAtomicsUsed();

        virtual const std::map<size_t, std::set<size_t> >& getAtomicsIndeps();

        /***********************************************************************
         * zero order (the original model)
         **********************************************************************/

        virtual void generateZeroSource(std::map<std::string, std::string>& sources);

        /**
         * Generates the operation graph for the zero order model with loops
         */
        virtual vector<CGBase> prepareForward0WithLoops(CodeHandler<Base>& handler,
                                                        const vector<CGBase>& x);

        /***********************************************************************
         * Jacobian
         **********************************************************************/

        virtual void generateJacobianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseJacobianSource(std::map<std::string, std::string>& sources) throw (CGException);

        virtual void generateSparseJacobianSource(std::map<std::string, std::string>& sources,
                                                  bool forward) throw (CGException);

        virtual void generateSparseJacobianForRevSource(std::map<std::string, std::string>& sources,
                                                        bool forward);

        /**
         * Generates a sparse Jacobian using loops.
         * 
         * The original model is split into two models: 
         *   - one for the repeated equations
         * \f[ y_i = f(x_{l(j)}, x_v, z_k) \f]
         *   - and another for the equations which do not belong in a loop and the 
         *   non-indexed temporary variables (\f$z\f$) used by \f$f\f$.
         * \f[ z_k = g_k(x_v) \f]
         * 
         * The jacobian elements for the equations in loops are evaluated as:
         * \f[ \frac{\mathrm{d} y_i}{\mathrm{d} x_v} = 
         *        \sum_k \left( \frac{\partial f_i}{\partial z_k} \frac{\partial z_k}{\partial x_v} \right) +
         *        \sum_j \left( \frac{\partial f_i}{\partial x_{l(j)}} \frac{\partial x_{l(j)}}{\partial x_v} \right) +
         *        \frac{\partial f_i}{\partial x_v} 
         * \f]
         * 
         * @param handler The operation graph handler
         * @param indVars The independent variables
         * @return the operation graph for the compressed jacobin with loops
         */
        virtual vector<CGBase> prepareSparseJacobianWithLoops(CodeHandler<Base>& handler,
                                                              const vector<CGBase>& x,
                                                              bool forward) throw (CGException);

        inline void prepareSparseJacobianRowWithLoops(CodeHandler<Base>& handler,
                                                      LoopModel<Base>& lModel,
                                                      size_t tapeI,
                                                      const loops::JacobianWithLoopsRowInfo& rowInfo,
                                                      const std::vector<std::map<size_t, CGBase> >& dyiDxtape,
                                                      const std::vector<std::map<size_t, CGBase> >& dzDx,
                                                      const CGBase& py,
                                                      IndexOperationNode<Base>& iterationIndexOp,
                                                      vector<loops::IfElseInfo<Base> >& ifElses,
                                                      size_t& jacLE,
                                                      vector<std::pair<CG<Base>, IndexPattern*> >& indexedLoopResults,
                                                      std::set<size_t>& allLocations);

        inline void analyseSparseJacobianWithLoops(const std::vector<size_t>& rows,
                                                   const std::vector<size_t>& cols,
                                                   const std::vector<size_t>& location,
                                                   vector<std::set<size_t> >& noLoopEvalSparsity,
                                                   vector<std::map<size_t, std::set<size_t> > >& noLoopEvalLocations,
                                                   std::map<LoopModel<Base>*, vector<std::set<size_t> > >& loopsEvalSparsities,
                                                   std::map<LoopModel<Base>*, std::vector<loops::JacobianWithLoopsRowInfo> >& loopEqInfo) throw (CGException);


        virtual void generateSparseJacobianWithLoopsSourceFromForRev(std::map<std::string, std::string>& sources,
                                                                     const std::map<size_t, std::vector<std::set<size_t> > >& userJacElLocation,
                                                                     const std::map<size_t, bool>& keyOrdered,
                                                                     size_t maxCompressedSize,
                                                                     const std::string& localFunctionTypeName,
                                                                     const std::string& suffix,
                                                                     const std::string& keyName,
                                                                     const std::map<size_t, std::set<size_t> >& nonLoopElements,
                                                                     const std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > >& loopGroups,
                                                                     void (*generateLocalFunctionName)(std::ostringstream& cache, const std::string& modelName, const LoopModel<Base>& loop, size_t g)) throw (CGException);

        inline virtual void generateFunctionNameLoopFor1(std::ostringstream& cache,
                                                         const LoopModel<Base>& loop,
                                                         size_t g);

        inline static void generateFunctionNameLoopFor1(std::ostringstream& cache,
                                                        const std::string& modelName,
                                                        const LoopModel<Base>& loop,
                                                        size_t g);

        inline virtual void generateFunctionNameLoopRev1(std::ostringstream& cache,
                                                         const LoopModel<Base>& loop,
                                                         size_t i);

        inline static void generateFunctionNameLoopRev1(std::ostringstream& cache,
                                                        const std::string& modelName,
                                                        const LoopModel<Base>& loop,
                                                        size_t i);

        /***********************************************************************
         * Hessian
         **********************************************************************/

        virtual void generateHessianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseHessianSource(std::map<std::string, std::string>& sources) throw (CGException);

        virtual void generateSparseHessianSourceDirectly(std::map<std::string, std::string>& sources) throw (CGException);

        virtual void generateSparseHessianSourceFromRev2(std::map<std::string, std::string>& sources) throw (CGException);

        virtual void determineSecondOrderElements4Eval(std::vector<size_t>& userRows,
                                                       std::vector<size_t>& userCols);

        /**
         * Loops
         */
        /**
         * Generates a sparse Hessian using loops.
         * 
         * The original model is split into two models: 
         *   - one for the repeated equations
         * \f[ y_i = f(x_{l(j)}, x_v, z_k) \f]
         *   - and another for the equations which do not belong in a loop and 
         *     the non-indexed temporary variables (\f$z\f$) used by \f$f\f$.
         * \f[ z_k = g_k(x_v) \f]
         * 
         * The Hessian elements for the equations in loops are evaluated as:
         * \f[ \frac{\mathrm{d}^2 y_i}{\partial x_w \partial x_v} = 
         *        \sum_k \left( \frac{\partial^2 f_i}{\partial x_w \partial z_k} \frac{\partial z_k}{\partial x_v} + 
         *                      \frac{\partial f_i}{\partial z_k} \frac{\partial^2 z_k}{\partial x_w \partial x_v}
         *               \right) +
         *        \sum_j \left( \frac{\partial^2 f_i}{\partial x_w \partial x_{l(j)}} \frac{\partial x_{l(j)}}{\partial x_v} \right) +
         *        \frac{\partial^2 f_i}{\partial x_w \partial x_v} 
         * \f]
         * 
         * @param handler The operation graph handler
         * @param indVars The independent variables
         * @return the operation graph for the compressed jacobin with loops
         */
        virtual vector<CGBase> prepareSparseHessianWithLoops(CodeHandler<Base>& handler,
                                                             vector<CGBase>& indVars,
                                                             vector<CGBase>& w,
                                                             const std::vector<size_t>& lowerHessRows,
                                                             const std::vector<size_t>& lowerHessCols,
                                                             const std::vector<size_t>& lowerHessOrder,
                                                             const std::map<size_t, size_t>& duplicates);

        inline void analyseSparseHessianWithLoops(const std::vector<size_t>& lowerHessRows,
                                                  const std::vector<size_t>& lowerHessCols,
                                                  const std::vector<size_t>& lowerHessOrder,
                                                  vector<std::set<size_t> >& noLoopEvalJacSparsity,
                                                  vector<std::set<size_t> >& noLoopEvalHessSparsity,
                                                  vector<std::map<size_t, std::set<size_t> > >& noLoopEvalHessLocations,
                                                  std::map<LoopModel<Base>*, loops::HessianWithLoopsInfo<Base> >& loopHessInfo,
                                                  bool useSymmetry);

        inline virtual void generateSparseHessianWithLoopsSourceFromRev2(std::map<std::string, std::string>& sources,
                                                                         const std::map<size_t, std::vector<std::set<size_t> > >& userHessElLocation,
                                                                         const std::map<size_t, bool>& ordered,
                                                                         size_t maxCompressedSize) throw (CGException);

        inline virtual void generateFunctionNameLoopRev2(std::ostringstream& cache,
                                                         const LoopModel<Base>& loop,
                                                         size_t g);

        static inline void generateFunctionNameLoopRev2(std::ostringstream& cache,
                                                        const std::string& modelName,
                                                        const LoopModel<Base>& loop,
                                                        size_t g);

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

        virtual void generateSparseForwardOneSourcesWithAtomics(std::map<std::string, std::string>& sources,
                                                                const std::map<size_t, std::vector<size_t> >& elements);

        virtual void generateSparseForwardOneSourcesNoAtomics(std::map<std::string, std::string>& sources,
                                                              const std::map<size_t, std::vector<size_t> >& elements);

        virtual void generateForwardOneSources(std::map<std::string, std::string>& sources);

        virtual void prepareSparseForwardOneWithLoops(std::map<std::string, std::string>& sources,
                                                      const std::map<size_t, std::vector<size_t> >& elements) throw (CGException);

        virtual void createForwardOneWithLoopsNL(CodeHandler<Base>& handler,
                                                 size_t j,
                                                 vector<CG<Base> >& jacCol,
                                                 std::map<std::string, std::string>& sources);


        inline static std::map<size_t, std::map<size_t, CG<Base> > > generateLoopFor1Jac(ADFun<CGBase>& fun,
                                                                                         const vector<std::set<size_t> >& sparsity,
                                                                                         const vector<std::set<size_t> >& evalSparsity,
                                                                                         const vector<CGBase>& xl,
                                                                                         bool constainsAtomics);

        /***********************************************************************
         * Reverse 1 mode
         **********************************************************************/

        virtual void generateSparseReverseOneSources(std::map<std::string, std::string>& sources);

        virtual void generateSparseReverseOneSourcesWithAtomics(std::map<std::string, std::string>& sources,
                                                                const std::map<size_t, std::vector<size_t> >& elements);

        virtual void generateSparseReverseOneSourcesNoAtomics(std::map<std::string, std::string>& sources,
                                                              const std::map<size_t, std::vector<size_t> >& elements);

        virtual void generateReverseOneSources(std::map<std::string, std::string>& sources);

        virtual void prepareSparseReverseOneWithLoops(std::map<std::string, std::string>& sources,
                                                      const std::map<size_t, std::vector<size_t> >& elements) throw (CGException);

        virtual void createReverseOneWithLoopsNL(CodeHandler<Base>& handler,
                                                 size_t i,
                                                 vector<CG<Base> >& jacRow,
                                                 std::map<std::string, std::string>& sources);

        inline static std::vector<std::map<size_t, CGBase> > generateLoopRev1Jac(ADFun<CGBase>& fun,
                                                                                 const vector<std::set<size_t> >& sparsity,
                                                                                 const vector<std::set<size_t> >& evalSparsity,
                                                                                 const vector<CGBase>& xl,
                                                                                 bool constainsAtomics);

        /***********************************************************************
         * Reverse 2 mode
         **********************************************************************/

        virtual void generateSparseReverseTwoSources(std::map<std::string, std::string>& sources) throw (CGException);

        virtual void generateSparseReverseTwoSourcesWithAtomics(std::map<std::string, std::string>& sources,
                                                                const std::map<size_t, std::vector<size_t> >& elements) throw (CGException);

        virtual void generateSparseReverseTwoSourcesNoAtomics(std::map<std::string, std::string>& sources,
                                                              const std::map<size_t, std::vector<size_t> >& elements,
                                                              const std::vector<size_t>& evalRows,
                                                              const std::vector<size_t>& evalCols) throw (CGException);

        virtual void generateReverseTwoSources(std::map<std::string, std::string>& sources);

        virtual void generateGlobalDirectionalFunctionSource(const std::string& function,
                                                             const std::string& function2_suffix,
                                                             const std::string& function_sparsity,
                                                             const std::map<size_t, std::vector<size_t> >& elements,
                                                             std::map<std::string, std::string>& sources);

        /**
         * Loops
         */
        virtual void prepareSparseReverseTwoWithLoops(std::map<std::string, std::string>& sources,
                                                      const std::map<size_t, std::vector<size_t> >& elements) throw (CGException);

        /***********************************************************************
         * Sparsities
         **********************************************************************/

        virtual void determineJacobianSparsity();

        virtual void generateJacobianSparsitySource(std::map<std::string, std::string>& sources);

        virtual void determineHessianSparsity();

        /**
         * Determines groups of rows from a sparsity pattern which do not share
         * the same columns
         * 
         * @param columns the column indexes of interest (all others are ignored);
         *                an empty set means all columns
         * @param sparsity The sparsity pattern to color
         * @return the colors
         */
        inline vector<CLangCompileModelHelper<Base>::Color> colorByRow(const std::set<size_t>& columns,
                                                                       const SparsitySetType& sparsity);

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
        /**
         * 
         */
        inline void startingJob(const std::string& jobName);

        inline void finishedJob();

    private:
        CLangCompileModelHelper(const CLangCompileModelHelper&); // not implemented

        CLangCompileModelHelper& operator=(const CLangCompileModelHelper&); // not implemented

        friend class
        CLangCompileDynamicHelper<Base>;
    };

}

#endif
