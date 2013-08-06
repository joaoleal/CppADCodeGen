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
        linear2Sections,
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

        template<class VectorSizeT>
        static inline IndexPattern* detect(const VectorSizeT& indexes);
    };

    /**
     * Linear pattern
     */
    class LinearIndexPattern : public IndexPattern {
    protected:
        long a_;
        long b_;
    public:

        inline LinearIndexPattern(long a, long b) :
            a_(a), b_(b) {
        }

        inline long getLinearSlope() const {
            return a_;
        }

        inline long getLinearConstantTerm() const {
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
    class Linear2SectionsIndexPattern : public IndexPattern {
    protected:
        LinearIndexPattern linear1_;
        LinearIndexPattern linear2_;
        size_t itSplit_; // the start of the second linear section
    public:

        inline Linear2SectionsIndexPattern(long a1, size_t b1,
                                           long a2, size_t b2,
                                           size_t itSplit) :
            linear1_(a1, b1),
            linear2_(a2, b2),
            itSplit_(itSplit) {
        }

        inline size_t getItrationSplit() const {
            return itSplit_;
        }

        const LinearIndexPattern& getLinearSection1() const {
            return linear1_;
        }

        const LinearIndexPattern& getLinearSection2() const {
            return linear2_;
        }

        inline virtual IndexPatternType getType() const {
            return linear2Sections;
        }

        inline virtual ~Linear2SectionsIndexPattern() {
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

    template<class VectorSizeT>
    IndexPattern* IndexPattern::detect(const VectorSizeT& indexes) {
        assert(indexes.size() > 1);

        long a = long(indexes[1]) - indexes[0];
        long b = indexes[0];
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

        // maybe there are 2 linear sections
        long a2 = long(indexes[lastLinear + 1]) - indexes[lastLinear];
        long b2 = indexes[lastLinear] - a2 * lastLinear;
        size_t lastLinear2 = indexes.size();
        for (size_t pos = lastLinear + 2; pos < indexes.size(); pos++) {
            if (indexes[pos] != a2 * pos + b2) {
                lastLinear = pos;
                break;
            }
        }
        if (lastLinear2 == indexes.size()) {
            return new Linear2SectionsIndexPattern(a, b, a2, b2, lastLinear);
        }

        return new RandomIndexPattern();
    }

}

#endif