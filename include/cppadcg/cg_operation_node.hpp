#ifndef CPPAD_CG_EXPRESSION_NODE_INCLUDED
#define CPPAD_CG_EXPRESSION_NODE_INCLUDED
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
     * An operation
     * 
     * @author Joao Leal
     */
    template<class Base>
    class OperationNode {
    public:
        static const std::set<CGOpCode> CUSTOM_NODE_CLASS;
    private:
        // the operations used to create this variable (temporary variables only)
        CGOpCode operation_;
        // additional operation information
        std::vector<size_t> info_;
        // the code blocks this block depends upon (empty for independent 
        // variables and possibly for the 1st assignment of a dependent variable)
        std::vector<Argument<Base> > arguments_;
        // variable ID that was altered/assigned in this source code
        // (zero means that no variable is assigned)
        size_t var_id_;
        //
        size_t evaluation_order_;
        // the total number of times the result of this operation is used
        size_t total_use_count_;
        // the number of times the result of this operation has been used
        size_t use_count_;
        // the last source code order in the call graph that uses the result of
        // this operation as an argument
        size_t last_usage_order_;
        //
        size_t color_;
        // generated variable name
        std::string* name_;
    public:

        inline OperationNode(CGOpCode op) :
            operation_(op),
            var_id_(0),
            evaluation_order_(0),
            total_use_count_(0),
            use_count_(0),
            last_usage_order_(0),
            color_(0),
            name_(nullptr) {
        }

        inline OperationNode(CGOpCode op,
                             const Argument<Base>& arg) :
            operation_(op),
            arguments_{arg},
            var_id_(0),
            evaluation_order_(0),
            total_use_count_(0),
            use_count_(0),
            last_usage_order_(0),
            color_(0),
            name_(NULL) {
        }

        inline OperationNode(CGOpCode op,
                             std::vector<Argument<Base> >&& args) :
            operation_(op),
            arguments_(args),
            var_id_(0),
            evaluation_order_(0),
            total_use_count_(0),
            use_count_(0),
            last_usage_order_(0),
            color_(0),
            name_(nullptr) {
        }

        inline OperationNode(CGOpCode op,
                             std::vector<size_t>&& info,
                             std::vector<Argument<Base> >&& args) :
            operation_(op),
            info_(info),
            arguments_(args),
            var_id_(0),
            evaluation_order_(0),
            total_use_count_(0),
            use_count_(0),
            last_usage_order_(0),
            color_(0),
            name_(nullptr) {
        }

        inline OperationNode(CGOpCode op,
                             const std::vector<size_t>& info,
                             const std::vector<Argument<Base> >& args) :
            operation_(op),
            info_(info),
            arguments_(args),
            var_id_(0),
            evaluation_order_(0),
            total_use_count_(0),
            use_count_(0),
            last_usage_order_(0),
            color_(0),
            name_(nullptr) {
        }

        inline void makeAlias(const Argument<Base>& other) {
            CPPADCG_ASSERT_UNKNOWN(CUSTOM_NODE_CLASS.find(operation_) == CUSTOM_NODE_CLASS.end()); // TODO: consider relaxing this check

            operation_ = CGAliasOp;
            arguments_.resize(1);
            arguments_[0] = other;
            var_id_ = 0;
            delete name_;
            name_ = nullptr;
        }

        inline CGOpCode getOperationType() const {
            return operation_;
        }

        inline void setOperation(CGOpCode op, const std::vector<Argument<Base> >& arguments = std::vector<Argument<Base> >()) {
            CPPADCG_ASSERT_UNKNOWN(op == operation_ || CUSTOM_NODE_CLASS.find(op) == CUSTOM_NODE_CLASS.end()); // cannot transform into a node with a custom class

            operation_ = op;
            arguments_ = arguments;
        }

        /**
         * Provides the arguments used in the operation represented by this
         * code fragment.
         * @return the arguments for the operation in this code fragment
         */
        inline const std::vector<Argument<Base> >& getArguments() const {
            return arguments_;
        }

        inline std::vector<Argument<Base> >& getArguments() {
            return arguments_;
        }

        /**
         * Provides additional information used in the operation.
         * @return the additional operation information/options
         */
        inline const std::vector<size_t>& getInfo() const {
            return info_;
        }

        inline std::vector<size_t>& getInfo() {
            return info_;
        }

        /**
         * Provides the variable ID that was altered/assigned in this source 
         * code (zero means that no variable is assigned).
         * @return the variable ID
         */
        inline size_t getVariableID() const {
            return var_id_;
        }

        /**
         * Specifies a variable ID to the result of this source code
         * (zero means that no variable is created).
         */
        inline void setVariableID(size_t var_id) {
            var_id_ = var_id;
        }

        inline size_t getEvaluationOrder() const {
            return evaluation_order_;
        }

        inline void setEvaluationOrder(size_t evaluation_order) {
            evaluation_order_ = evaluation_order;
        }

        /**
         * Provides the total number of times the result of this operation is
         * being used as an argument for another operation.
         * @return the total usage count
         */
        inline size_t getTotalUsageCount() const {
            return total_use_count_;
        }

        inline void setTotalUsageCount(size_t cout) {
            total_use_count_ = cout;
        }

        inline void increaseTotalUsageCount() {
            total_use_count_++;
        }

        inline void resetTotalUsageCount() {
            total_use_count_ = 0;
        }

        /**
         * Provides the number of times the result of this operation has been 
         * used as an argument for another operation.
         * 
         * @return the current usage count
         */
        inline size_t getUsageCount() const {
            return use_count_;
        }

        inline void increaseUsageCount() {
            use_count_++;
        }

        inline void resetUsageCount() {
            use_count_ = 0;
        }

        inline size_t getLastUsageEvaluationOrder() const {
            return last_usage_order_;
        }

        inline void setLastUsageEvaluationOrder(size_t last) {
            last_usage_order_ = last;
            if (operation_ == CGArrayElementOp) {
                OperationNode<Base>* array = arguments_[0].getOperation();
                CPPADCG_ASSERT_UNKNOWN(array->getOperationType() == CGArrayCreationOp);
                if (array->getLastUsageEvaluationOrder() < last) {
                    array->setLastUsageEvaluationOrder(last);
                }
            } else if (operation_ == CGTmpOp) {
                OperationNode<Base>* declr = arguments_[0].getOperation();
                CPPADCG_ASSERT_UNKNOWN(declr->getOperationType() == CGTmpDclOp);
                if (declr->getLastUsageEvaluationOrder() < last) {
                    declr->setLastUsageEvaluationOrder(last);
                }
            }
        }

        inline void resetHandlerCounters() {
            total_use_count_ = 0;
            use_count_ = 0;
            var_id_ = 0;
            evaluation_order_ = 0;
            last_usage_order_ = 0;
        }

        inline const std::string* getName() const {
            return name_;
        }

        inline void setName(const std::string& name) {
            if (name_ != nullptr)
                *name_ = name;
            else
                name_ = new std::string(name);
        }

        inline void clearName() {
            delete name_;
            name_ = nullptr;
        }

        inline size_t getColor() const {
            return color_;
        }

        inline void setColor(size_t color) {
            color_ = color;
        }

        inline virtual ~OperationNode() {
            delete name_;
        }

    private:

        static inline std::set<CGOpCode> makeCustomNodeClassesSet();

        OperationNode(const OperationNode& orig) :
            operation_(orig.operation_),
            info_(orig.info_),
            arguments_(orig.arguments_),
            var_id_(0),
            evaluation_order_(0),
            total_use_count_(0),
            use_count_(0),
            last_usage_order_(orig.last_usage_order_),
            color_(orig.color_),
            name_(orig.name_ != nullptr ? new std::string(*orig.name_) : nullptr) {
        }

        friend class CodeHandler<Base>;

    };

    template<class Base>
    inline std::set<CGOpCode> OperationNode<Base>::makeCustomNodeClassesSet() {
        std::set<CGOpCode> s;
        s.insert(CGIndexAssignOp);
        s.insert(CGIndexOp);
        s.insert(CGLoopStartOp);
        s.insert(CGLoopEndOp);
        s.insert(CGPriOp);
        return s;
    }

    template<class Base>
    const std::set<CGOpCode> OperationNode<Base>::CUSTOM_NODE_CLASS = makeCustomNodeClassesSet();

    template<class Base>
    inline std::ostream& operator <<(
            std::ostream& os, //< stream to write to
            const CppAD::OperationNode<Base>& c) {
        CGOpCode op = c.getOperationType();
        switch (op) {
            case CGArrayCreationOp:
                os << "new $1[" << c.getArguments().size() << "]";
                break;
            case CGSparseArrayCreationOp:
                os << "new $1[" << c.getInfo()[0] << "]";
                break;
            case CGArrayElementOp:
                os << "$1[" << c.getInfo()[0] << "]";
                break;
            case CGAtomicForwardOp:
                os << "atomicFunction.forward(" << c.getInfo()[0] << ", " << c.getInfo()[1] << ", vx, vy, $1, $2)";
                break;
            case CGAtomicReverseOp:
                os << "atomicFunction.reverse(" << c.getInfo()[0] << ", $1, $2, $3, $4)";
                break;
            case CGSignOp:
                os << "if($1 > 0) { 1 } else if($1 == 0) { 0 } else { -1 }";
                break;

            default:
                os << op;
        }

        return os;
    }

}

#endif
