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
        /**
         * hessian sparsity from the model for each equation
         */
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
        std::map<LoopModel<Base>*, std::map<TapeVarType, vector<GroupLoopRev2ColInfo<Base>* > > > _loopRev2Groups;
        typedef std::pair<size_t, size_t> Compressed2JColType;
        // elements[var]{compressed location, original column index}
        std::map<size_t, std::vector<Compressed2JColType> > _nonLoopRev2Elements;
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

            typename std::set<LoopModel<Base>* >::const_iterator it;
            for (it = _loopTapes.begin(); it != _loopTapes.end(); ++it) {
                delete *it;
            }

            typename std::map<LoopModel<Base>*, std::map<TapeVarType, vector<GroupLoopRev2ColInfo<Base>* > > >::const_iterator itljg;
            for (itljg = _loopRev2Groups.begin(); itljg != _loopRev2Groups.end(); ++itljg) {
                typename std::map<TapeVarType, vector<GroupLoopRev2ColInfo<Base>* > >::const_iterator itjg;
                for (itjg = itljg->second.begin(); itjg != itljg->second.end(); ++itjg) {
                    const vector<GroupLoopRev2ColInfo<Base>*>& groups = itjg->second;
                    for (size_t g = 0; g < groups.size(); g++) {
                        delete groups[g];
                    }
                }
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

        /***********************************************************************
         * zero order (the orginal model)
         **********************************************************************/

        virtual void generateZeroSource(std::map<std::string, std::string>& sources);

        /**
         * Generates the operation graph for the zero order model with loops
         */
        virtual vector<CGBase> prepareForward0WithLoops(CodeHandler<Base>& handler,
                                                        const vector<CGBase>& x);

        inline vector<CGBase> createIndexedIndependents(CodeHandler<Base>& handler,
                                                        LoopModel<Base>& loop,
                                                        IndexOperationNode<Base>& iterationIndexOp);

        inline vector<CGBase> createLoopIndependentVector(CodeHandler<Base>& handler,
                                                          LoopModel<Base>& loop,
                                                          const vector<CGBase>& indexedIndeps,
                                                          const vector<CGBase>& nonIndexedIndeps,
                                                          const vector<CGBase>& nonIndexedTmps);
        /***********************************************************************
         * Jacobian
         **********************************************************************/

        virtual void generateJacobianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseJacobianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseJacobianSource(std::map<std::string, std::string>& sources,
                                                  bool forward);

        virtual void generateSparseJacobianForRevSource(std::map<std::string, std::string>& sources,
                                                        bool forward);

        /**
         * Estimates the work load of forward vs reverse mode
         * 
         * @return true if the foward mode should be used, false for the reverse mode
         */
        static bool estimateBestJacobianADMode(const std::vector<size_t>& jacRows,
                                               const std::vector<size_t>& jacCols);
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

        virtual void analyseSparseHessianWithLoops(const std::vector<size_t>& lowerHessRows,
                                                   const std::vector<size_t>& lowerHessCols,
                                                   const std::vector<size_t>& lowerHessOrder,
                                                   vector<std::set<size_t> >& noLoopEvalJacSparsity,
                                                   vector<std::set<size_t> >& noLoopEvalHessSparsity,
                                                   vector<std::map<size_t, std::set<size_t> > >& noLoopEvalHessLocations,
                                                   std::map<LoopModel<Base>*, loops::HessianWithLoopsInfo<Base> >& loopHessInfo);

        void generateGlobalReverseTwoWithLoopsFunctionSource(const std::map<size_t, std::vector<size_t> >& elements,
                                                             std::map<std::string, std::string>& sources);

        inline virtual void generateSparseHessianWithLoopsSourceFromRev2(std::map<std::string, std::string>& sources,
                                                                         const std::map<size_t, std::vector<std::set<size_t> > >& userHessElLocation,
                                                                         const std::map<size_t, bool>& ordered,
                                                                         size_t maxCompressedSize);

        inline virtual void generateFunctionDeclarationSourceLoopRev2(std::ostringstream& cache,
                                                                      CLanguage<Base>& langC);

        inline virtual void generateFunctionNameLoopRev2(std::ostringstream& cache,
                                                         const LoopModel<Base>& loop,
                                                         const TapeVarType& jTape1,
                                                         size_t g);

        inline static vector<CG<Base> > createLoopDependentVector(CodeHandler<Base>& handler,
                                                                  LoopModel<Base>& loop,
                                                                  IndexOperationNode<Base>& iterationIndexOp);

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

        template<class T>
        void generateFunctionDeclarationSource(std::ostringstream& cache,
                                               const std::string& model_function,
                                               const std::string& suffix,
                                               const std::map<size_t, T>& elements,
                                               const std::string& argsDcl);

        /**
         * Loops
         */
        virtual void prepareSparseReverseTwoWithLoops(std::map<std::string, std::string>& sources,
                                                      const std::map<size_t, std::vector<size_t> >& elements);

        virtual void prepareSparseReverseTwoSourcesForLoop(std::map<std::string, std::string>& sources,
                                                           CodeHandler<Base>& handler,
                                                           LoopModel<Base>& loop,
                                                           std::map<size_t, std::vector<LoopRev2ValInfo<Base> > >& hess,
                                                           const CGBase& tx1);
        /*
                std::string generateSparseReverseTwoWithLoopsVarGroupSource(const std::string& functionName,
                                                                            const std::string& jobName,
                                                                            LoopModel<Base>& loop,
                                                                            CodeHandler<Base>& handler,
                                                                            const std::pair<size_t, size_t>& jTape1,
                                                                            const GroupLoopRev2ColInfo<Base>& group,
                                                                            const Index& indexJrow,
                                                                            const Index& indexLocalIt,
                                                                            const Index& indexLocalItCount,
                                                                            const IndexPattern& itPattern,
                                                                            const IndexPattern* itCountPattern,
                                                                            const std::map<TapeVarType, Plane2DIndexPattern*>& loopDepIndexes,
                                                                            std::map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> >& evaluations,
                                                                            const CGBase& tx1);
         */
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

        /***********************************************************************
         * Loops
         **********************************************************************/
        /*
                static inline void prepareLoops(CodeHandler<Base>& handler,
                                                std::vector<CGBase>& jac,
                                                std::map<LoopModel<Base>*, std::map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> > >& evaluations,
                                                std::map<LoopModel<Base>*, vector<IndexedDependentLoopInfo<Base>* > >& dependentIndexes,
                                                size_t assignOrAdd = 0);
         */
        static inline LoopEndOperationNode<Base>* createLoopEnd(CodeHandler<Base>& handler,
                                                                LoopStartOperationNode<Base>& loopStart,
                                                                const vector<std::pair<CG<Base>, IndexPattern*> >& indexedLoopResults,
                                                                const std::set<IndexOperationNode<Base>*>& indexesOps,
                                                                const LoopNodeInfo<Base>& loopInfo,
                                                                size_t assignOrAdd);

        static inline void moveNonIndexedOutsideLoop(LoopStartOperationNode<Base>& loopStart,
                                                     LoopEndOperationNode<Base>& loopEnd,
                                                     const Index& loopIndex);

        static inline bool findNonIndexedNodes(OperationNode<Base>& node,
                                               std::set<OperationNode<Base>*>& nonIndexed,
                                               const Index& loopIndex);

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
