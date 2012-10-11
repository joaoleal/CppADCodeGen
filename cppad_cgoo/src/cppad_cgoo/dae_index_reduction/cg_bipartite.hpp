#ifndef CPPAD_CG_BIPARTITE_INCLUDED
#define	CPPAD_CG_BIPARTITE_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

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
    public:

        inline Enode(size_t index, Enode<Base>* differentiationOf = NULL, Enode<Base>* diff = NULL) :
            BiPGraphNode<Base>(index),
            differentiation_(diff),
            differentiationOf_(differentiationOf) {

            if (differentiationOf_ != NULL) {
                differentiationOf_->setDerivative(this);
            }
        }

        inline const std::set<Vnode<Base>*>& variables() const {
            return vnodes_;
        }

        inline const std::set<Vnode<Base>*>& originalVariables() const {
            return vnodes_orig_;
        }

        inline void addVariable(Vnode<Base>* j) {
            vnodes_orig_.insert(j);
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

        inline void deleteNode(Vnode<Base>* j) {
            vnodes_.erase(j);
        }

        virtual std::string nodeType() {
            return TYPE;
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
        bool deleted_;
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
    public:

        inline Vnode(size_t index, Vnode<Base>* derivativeOf = NULL, Vnode<Base>* diff = NULL) :
            BiPGraphNode<Base>(index),
            deleted_(false),
            assign_(NULL),
            derivative_(diff),
            derivativeOf_(derivativeOf) {

            if (derivativeOf_ != NULL) {
                derivativeOf_->setDerivative(this);
            }
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

        inline bool isDeleted() const {
            return deleted_;
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
            os << "Variable " << j.index();
        }
        return os;
    }

    template<class Base>
    const std::string Vnode<Base>::TYPE = "Variable";

}

#endif

