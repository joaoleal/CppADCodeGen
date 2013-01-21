#ifndef CPPAD_CG_BIPARTITE_INCLUDED
#define CPPAD_CG_BIPARTITE_INCLUDED
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

#include <assert.h>
#include <iostream>
#include <set>

namespace CppAD {

    /**
     * Bipartite graph node
     */
    template<class Base>
    class BiPGraphNode {
    protected:
        size_t index_; // location of node
        bool colored_; // node visited
    public:

        inline BiPGraphNode(size_t index) :
            index_(index),
            colored_(false) {
        }

        inline void color() {
            colored_ = true;

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "      Coloured " << nodeType() << " " << index_ << "\n";
#endif
        }

        inline void uncolor() {
            colored_ = false;
        }

        inline bool isColored() const {
            return colored_;
        }

        inline size_t index() const {
            return index_;
        }

        virtual std::string nodeType() = 0;

        inline virtual ~BiPGraphNode() {
        }
    };

    template<class Base>
    class Vnode; // forward declaration

    /**
     * Equation nodes
     */
    template<class Base>
    class Enode : public BiPGraphNode<Base> {
    protected:
        static const std::string TYPE;
        // original variables present in this equation which where not deleted
        std::set<Vnode<Base>*> vnodes_orig_;
        // variables present in this equation which where not deleted
        std::set<Vnode<Base>*> vnodes_;
        /**
         * the differentiated equation used to produce this one
         *  (B in Pantelides algorithm)
         */
        Enode<Base>* differentiation_;
        /**
         * 
         */
        Enode<Base>* differentiationOf_;
        /**
         * true if this equation is a purely algebraic equation
         */
        bool algebraic_;
    public:

        inline Enode(size_t index) :
            BiPGraphNode<Base>(index),
            differentiation_(NULL),
            differentiationOf_(NULL),
            algebraic_(true) {
        }

        inline Enode(size_t index, Enode<Base>* differentiationOf) :
            BiPGraphNode<Base>(index),
            differentiation_(NULL),
            differentiationOf_(differentiationOf),
            algebraic_(true) {

            differentiationOf_->setDerivative(this);
        }

        inline bool isAlgebraic() const {
            return algebraic_;
        }

        inline const std::set<Vnode<Base>*>& variables() const {
            return vnodes_;
        }

        inline const std::set<Vnode<Base>*>& originalVariables() const {
            return vnodes_orig_;
        }

        inline void addVariable(Vnode<Base>* j) {
            vnodes_orig_.insert(j);
            if (j->derivativeOf() != NULL) {
                algebraic_ = false;
            }
            if (!j->isDeleted()) {
                vnodes_.insert(j);
                j->addEquation(this);
            }
        }

        /**
         * \return the equation that was derived by differentiating this 
         * equation
         */
        inline Enode<Base>* derivative() const {
            return differentiation_;
        }

        inline Enode<Base>* derivativeOf() const {
            return differentiationOf_;
        }

        inline Enode<Base>* originalEquation() {
            if (differentiationOf_ == NULL) {
                return this;
            } else {
                return differentiationOf_->originalEquation();
            }
        }

        inline void deleteNode(Vnode<Base>* j) {
            vnodes_.erase(j);
        }

        virtual std::string nodeType() {
            return TYPE;
        }

        inline virtual ~Enode() {
        }

    protected:

        inline void setDerivative(Enode<Base>* difEq) {
            differentiation_ = difEq;
        }
    };

    template<class Base>
    inline std::ostream& operator <<(std::ostream& os, const Enode<Base>& i) {
        if (i.derivativeOf() != NULL) {
            os << "Diff(" << *i.derivativeOf() << ")";
        } else {
            os << "Equation " << i.index();
        }

        return os;
    }

    template<class Base>
    const std::string Enode<Base>::TYPE = "Equation";

    /**
     * Variable nodes
     */
    template<class Base>
    class Vnode : public BiPGraphNode<Base> {
    protected:
        static const std::string TYPE;
        /**
         * 
         */
        bool deleted_;
        /**
         * Whether or not this variable is time dependent
         */
        bool parameter_;
        /**
         * equations that use this variable
         */
        std::set<Enode<Base>*> enodes_;
        /**
         * 
         */
        Enode<Base>* assign_;
        /**
         *  the time derivative variable of this variable
         *  (A in Pantelides algorithm)
         */
        Vnode<Base>* derivative_;
        /**
         *  the variable which was differentiated to create this one
         */
        Vnode<Base>* derivativeOf_;
        /**
         * The index in the tape
         */
        size_t tapeIndex_;
        /**
         * name
         */
        std::string name_;

    public:

        inline Vnode(size_t index, int tapeIndex, const std::string& name) :
            BiPGraphNode<Base>(index),
            deleted_(false),
            parameter_(false),
            assign_(NULL),
            derivative_(NULL),
            derivativeOf_(NULL),
            tapeIndex_(tapeIndex),
            name_(name) {

        }

        inline Vnode(size_t index, size_t tapeIndex, Vnode<Base>* derivativeOf, const std::string& name = "") :
            BiPGraphNode<Base>(index),
            deleted_(false),
            parameter_(false),
            assign_(NULL),
            derivative_(NULL),
            derivativeOf_(derivativeOf),
            tapeIndex_(tapeIndex),
            name_(name.empty() ? "d" + derivativeOf->name() + "dt" : name) {
            assert(derivativeOf_ != NULL);

            derivativeOf_->setDerivative(this);
        }

        const std::string& name() const {
            return name_;
        }

        size_t tapeIndex() const {
            return tapeIndex_;
        }

        inline const std::set<Enode<Base>*>& equations() const {
            return enodes_;
        }

        /**
         * \return the time derivative variable
         */
        inline Vnode<Base>* derivative() const {
            return derivative_;
        }

        /**
         * \return the variable which was differentiated to create this one
         */
        inline Vnode<Base>* derivativeOf() const {
            return derivativeOf_;
        }

        inline Vnode<Base>* originalVariable() {
            if (derivativeOf_ == NULL) {
                return this;
            } else {
                return derivativeOf_->originalVariable();
            }
        }

        inline Vnode<Base>* originalVariable(size_t origVarSize) {
            if (derivativeOf_ == NULL || this->index_ < origVarSize) {
                return this;
            } else {
                return derivativeOf_->originalVariable();
            }
        }

        inline bool isDeleted() const {
            return deleted_;
        }

        inline void makeParameter() {
            parameter_ = true;
            deleteNode();
        }

        inline bool isParameter() const {
            return parameter_;
        }

        inline void deleteNode() {
#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "Deleting " << *this << "\n";
#endif

            deleted_ = true;
            for (typename std::set<Enode<Base>*>::iterator i = enodes_.begin(); i != enodes_.end(); ++i) {
                (*i)->deleteNode(this);
            }
            enodes_.clear();
        }

        inline Enode<Base>* assigmentEquation() const {
            return assign_;
        }

        inline void setAssigmentEquation(Enode<Base>& i) {
#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "      Assigning " << *this << " to " << i << "\n";
#endif
            assign_ = &i;
        }

        virtual std::string nodeType() {
            return TYPE;
        }

        unsigned int order() const {
            if (derivativeOf_ == NULL) {
                return 0u;
            } else {
                return derivativeOf_->order() + 1u;
            }
        }

        inline virtual ~Vnode() {
        }

    protected:

        inline void addEquation(Enode<Base>* i) {
            if (!deleted_) {
                enodes_.insert(i);
            }
        }

        inline void setDerivative(Vnode<Base>* div) {
            derivative_ = div;
        }

        friend class Enode<Base>;
    };

    template<class Base>
    inline std::ostream& operator <<(std::ostream& os, const Vnode<Base>& j) {
        if (j.derivativeOf() != NULL) {
            os << "Diff(" << *j.derivativeOf() << ")";
        } else {
            os << "Variable " << j.name();
        }
        return os;
    }

    template<class Base>
    const std::string Vnode<Base>::TYPE = "Variable";

}

#endif

