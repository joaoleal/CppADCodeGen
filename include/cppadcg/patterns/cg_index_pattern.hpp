#ifndef CPPAD_CG_INDEX_PATTERN_INCLUDED
#define CPPAD_CG_INDEX_PATTERN_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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
     * Index pattern types
     */
    enum IndexPatternType {
        linear,
        linearNConst,
        random
    };

    /**
     * Generic index pattern
     */
    class IndexPattern {
    public:
        virtual IndexPatternType getType() const = 0;

        inline virtual ~IndexPattern() {
        }

        static inline IndexPattern* detect(const vector<size_t>& indexes);
    };

    /**
     * Linear pattern
     */
    class LinearIndexPattern : public IndexPattern {
    protected:
        long a_;
        size_t b_;
    public:

        inline LinearIndexPattern(long a, size_t b) :
            a_(a), b_(b) {
        }

        inline long getLinearSlope() const {
            return a_;
        }

        inline size_t getLinearConstantTerm() const {
            return b_;
        }

        inline virtual IndexPatternType getType() const {
            return linear;
        }

        inline virtual ~LinearIndexPattern() {
        }
    };

    /**
     * Linear pattern followed by a constant index
     */
    class LinearNconstIndexPattern : public LinearIndexPattern {
    protected:
        size_t max_; // end of the linear section
        size_t i_; // value after max
    public:

        inline LinearNconstIndexPattern(long a, size_t b, size_t max, size_t i) :
            LinearIndexPattern(a, b),
            max_(max), i_(i) {
        }

        inline size_t getLinearSectionEnd() const {
            return max_;
        }

        inline size_t getConstSectionValue() const {
            return i_;
        }

        inline virtual IndexPatternType getType() const {
            return linearNConst;
        }

        inline virtual ~LinearNconstIndexPattern() {
        }
    };

    /**
     * Random pattern
     */
    class RandomIndexPattern : public IndexPattern {
    public:

        inline virtual IndexPatternType getType() const {
            return random;
        }

        inline virtual ~RandomIndexPattern() {
        }
    };

    IndexPattern* IndexPattern::detect(const vector<size_t>& indexes) {
        assert(indexes.size() > 1);

        size_t b = indexes[0];
        long a = long(indexes[1]) - indexes[0];
        size_t lastLinear = indexes.size();
        for (size_t pos = 2; pos < indexes.size(); pos++) {
            if (indexes[pos] != a * pos + b) {
                lastLinear = pos;
                break;
            }
        }

        if (lastLinear == indexes.size()) {
            return new LinearIndexPattern(a, b);
        }

        // maybe it is constant after a linear part
        bool constIndex = true;
        for (size_t pos = lastLinear; pos < indexes.size(); pos++) {
            if (indexes[pos] != indexes[lastLinear]) {
                constIndex = false;
                break;
            }
        }

        if (constIndex) {
            return new LinearNconstIndexPattern(a, b, lastLinear, indexes[lastLinear]);
        }

        return new RandomIndexPattern();
    }

}

#endif