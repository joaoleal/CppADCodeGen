#ifndef CPPAD_CG_C_LANGUAGE_INDEX_PATTERNS_INCLUDED
#define CPPAD_CG_C_LANGUAGE_INDEX_PATTERNS_INCLUDED
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

    template<class Base>
    inline void CLanguage<Base>::createIndexDeclaration() {
        if (_indexes->empty())
            return;
        
        _ss << _spaces << "unsigned long";
        std::set<const Index*>::const_iterator iti;
        for (iti = _indexes->begin(); iti != _indexes->end(); ++iti) {
            if (iti != _indexes->begin())
                _ss << ",";
            _ss << " " << (*iti)->getName();
        }
        _ss << ";\n";
    }

    template<class Base>
    inline std::string CLanguage<Base>::createIndexPattern(const IndexPattern& ip) {
        std::stringstream ss;
        switch (ip.getType()) {
            case LINEAR: // y = x * a + b
            {
                const LinearIndexPattern& lip = static_cast<const LinearIndexPattern&> (ip);
                return createLinearIndexPattern(lip);
            }
            case LINEAR2: // y = x / a + b
            {
                const Linear2IndexPattern& lip = static_cast<const Linear2IndexPattern&> (ip);
                return createLinear2IndexPattern(lip);
            }
            case SECTIONED:
            {
                const SectionedIndexPattern* lip = static_cast<const SectionedIndexPattern*> (&ip);
                const std::map<size_t, IndexPattern*>& sections = lip->getLinearSections();
                size_t sSize = sections.size();
                assert(sSize > 1);

                std::map<size_t, IndexPattern*>::const_iterator its = sections.begin();
                for (size_t s = 0; s < sSize - 1; s++) {
                    const IndexPattern* lp = its->second;
                    ++its;
                    size_t xStart = its->first;

                    assert(!lp->getIndexes().empty());
                    ss << "(" << lp->getIndexes()[0]->getName() << "<" << xStart << ")? "
                            << createIndexPattern(*lp) << ": ";
                }
                ss << createIndexPattern(*its->second);

                return ss.str();
            }

            case PLANE2D: // y = f(x) + f(z)
            {
                std::string index;
                const Plane2DIndexPattern& pip = static_cast<const Plane2DIndexPattern&> (ip);
                if (pip.getPattern1() != NULL) {
                    index += createIndexPattern(*pip.getPattern1());
                }

                if (pip.getPattern2() != NULL) {
                    if (pip.getPattern1() != NULL)
                        index += " + ";
                    index += createIndexPattern(*pip.getPattern2()); // indexName!!!!!!!!!
                    throw CGException("Not implemented yet");
                }

                return index;
            }
            case RANDOM:
                throw CGException("Random indexes not implemented yet");
                //return ss.str();
            default:
                assert(false); // should never reach this
                return "";
        }
    }

    template<class Base>
    inline std::string CLanguage<Base>::createLinearIndexPattern(const LinearIndexPattern& lip) {
        assert(lip.getIndexes().size() == 1);

        std::stringstream ss;
        if (lip.getLinearSlope() != 0) {
            if (lip.getLinearSlope() != 1) {
                ss << lip.getLinearSlope() << " * ";
            }
            ss << lip.getIndexes()[0]->getName();
        }

        if (lip.getLinearConstantTerm() != 0) {
            if (lip.getLinearSlope() != 0)
                ss << " + ";
            ss << lip.getLinearConstantTerm();
        }
        return ss.str();
    }

    template<class Base>
    inline std::string CLanguage<Base>::createLinear2IndexPattern(const Linear2IndexPattern& lip) {
        assert(lip.getIndexes().size() == 1);

        std::stringstream ss;

        ss << lip.getIndexes()[0]->getName();
        if (lip.getInverseLinearSlope() != 1) {
            ss << " / " << lip.getInverseLinearSlope();
        }

        if (lip.getLinearConstantTerm() != 0) {
            ss << " + " << lip.getLinearConstantTerm();
        }
        return ss.str();
    }
}

#endif