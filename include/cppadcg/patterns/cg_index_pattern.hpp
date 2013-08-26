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
     * Generic index pattern
     */
    class IndexPattern {
    public:
        virtual IndexPatternType getType() const = 0;

        inline virtual ~IndexPattern() {
        }

        template<class VectorSizeT>
        static inline IndexPattern* detect(const VectorSizeT& indexes);

        IndexPattern* detect(const std::map<size_t, size_t>& indexes);
    };

    /**
     * Linear pattern y = a * x + b
     */
    class LinearIndexPattern : public IndexPattern {
    protected:
        // slope
        long a_;
        // constant term
        long b_;
    public:

        inline LinearIndexPattern() {
        }

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
            return LINEAR;
        }

        inline virtual ~LinearIndexPattern() {
        }
    };

    /**
     * Linear pattern followed by a constant index
     */
    class LinearSectionsIndexPattern : public IndexPattern {
    protected:
        /**
         * maps the start of the linear section (first x) to the linear pattern
         */
        std::map<size_t, LinearIndexPattern> sections_;
    public:

        inline LinearSectionsIndexPattern(const std::map<size_t, LinearIndexPattern>& sections) :
            sections_(sections) {
        }

        const std::map<size_t, LinearIndexPattern>& getLinearSections() const {
            return sections_;
        }

        inline virtual IndexPatternType getType() const {
            return LINEARSECTIONS;
        }

        inline virtual ~LinearSectionsIndexPattern() {
        }
    };

    /**
     * Random pattern
     */
    class RandomIndexPattern : public IndexPattern {
    public:

        inline virtual IndexPatternType getType() const {
            return RANDOM;
        }

        inline virtual ~RandomIndexPattern() {
        }
    };

    template<class VectorSizeT>
    IndexPattern* IndexPattern::detect(const VectorSizeT& indexes) {
        assert(indexes.size() > 1);

        std::map<size_t, LinearIndexPattern> linearSections;
        size_t xStart = 0;
        while (xStart != indexes.size()) {
            long a, b;
            size_t lastLinear;
            if (xStart + 1 == indexes.size()) {
                a = 0;
                b = indexes[xStart];
                lastLinear = xStart + 1;
            } else {
                a = long(indexes[xStart + 1]) - indexes[xStart];
                b = long(indexes[xStart]) - a * xStart;
                lastLinear = indexes.size();
                for (size_t x = xStart + 2; x < indexes.size(); x++) {
                    if (indexes[x] != a * x + b) {
                        lastLinear = x;
                        break;
                    }
                }
            }

            linearSections[xStart] = LinearIndexPattern(a, b);
            xStart = lastLinear;
        }

        if (linearSections.size() == 1) {
            return new LinearIndexPattern(linearSections.begin()->second);
        } else if (linearSections.size() <= 3 ||
                (linearSections.size() < indexes.size() / 4 && linearSections.size() < 10)) {
            return new LinearSectionsIndexPattern(linearSections);
        } else {
            throw CGException("Random index patterns not implemented yet!");
            return new RandomIndexPattern();
        }

    }

    IndexPattern* IndexPattern::detect(const std::map<size_t, size_t>& indexes) {
        assert(!indexes.empty());

        std::map<size_t, LinearIndexPattern> linearSections;

        std::map<size_t, size_t>::const_iterator pStart = indexes.begin();
        while (pStart != indexes.end()) {
            std::map<size_t, size_t>::const_iterator pNextSection = indexes.end();
            std::map<size_t, size_t>::const_iterator p1 = pStart;
            ++p1;
            long a, b;
            if (p1 == indexes.end()) {
                a = 0;
                b = pStart->second;
                pNextSection = p1;
            } else {
                a = long(p1->first) - pStart->first;
                b = long(pStart->second) - a * pStart->first;

                for (std::map<size_t, size_t>::const_iterator itp = p1; itp != indexes.end(); ++itp) {
                    size_t x = itp->first;
                    size_t y = itp->second;
                    if (y != a * x + b) {
                        pNextSection = itp;
                        break;
                    }
                }
            }

            linearSections[pStart->first] = LinearIndexPattern(a, b);
            pStart = pNextSection;
        }

        if (linearSections.size() == 1) {
            return new LinearIndexPattern(linearSections.begin()->second);
        } else if (linearSections.size() <= 3 ||
                (linearSections.size() < indexes.size() / 4 && linearSections.size() < 10)) {
            return new LinearSectionsIndexPattern(linearSections);
        } else {
            throw CGException("Random index patterns not implemented yet!");
            return new RandomIndexPattern();
        }
    }

}

#endif