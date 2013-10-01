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

        std::set<const IndexDclrOperationNode<Base>*> funcArgs(_funcArgIndexes.begin(), _funcArgIndexes.end());

        bool first = true;

        _ss << _spaces << "unsigned long";
        typename std::set<const IndexDclrOperationNode<Base>*>::const_iterator iti;
        for (iti = _indexes->begin(); iti != _indexes->end(); ++iti) {

            if (funcArgs.find(*iti) == funcArgs.end()) {
                if (first) first = false;
                else _ss << ",";

                _ss << " " << (*(*iti)->getName());
            }

        }
        _ss << ";\n";
    }

    template<class Base>
    inline std::string CLanguage<Base>::indexPattern2String(const IndexPattern& ip,
                                                           const IndexDclrOperationNode<Base>& index) {
        std::vector<const IndexDclrOperationNode<Base>*>indexes(1);
        indexes[0] = &index;
        return indexPattern2String(ip, indexes);
    }

    template<class Base>
    inline std::string CLanguage<Base>::indexPattern2String(const IndexPattern& ip,
                                                            const std::vector<const IndexDclrOperationNode<Base>*>& indexes) {
        std::stringstream ss;
        switch (ip.getType()) {
            case LINEAR: // y = x * a + b
            {
                CPPADCG_ASSERT_KNOWN(indexes.size() == 1, "Invalid number of indexes");
                const LinearIndexPattern& lip = static_cast<const LinearIndexPattern&> (ip);
                return linearIndexPattern2String(lip, *indexes[0]);
            }
            case SECTIONED:
            {
                CPPADCG_ASSERT_KNOWN(indexes.size() == 1, "Invalid number of indexes");
                const SectionedIndexPattern* lip = static_cast<const SectionedIndexPattern*> (&ip);
                const std::map<size_t, IndexPattern*>& sections = lip->getLinearSections();
                size_t sSize = sections.size();
                assert(sSize > 1);

                std::map<size_t, IndexPattern*>::const_iterator its = sections.begin();
                for (size_t s = 0; s < sSize - 1; s++) {
                    const IndexPattern* lp = its->second;
                    ++its;
                    size_t xStart = its->first;

                    ss << "(" << (*indexes[0]->getName()) << "<" << xStart << ")? "
                            << indexPattern2String(*lp, *indexes[0]) << ": ";
                }
                ss << indexPattern2String(*its->second, *indexes[0]);

                return ss.str();
            }

            case PLANE2D: // y = f(x) + f(z)
            {
                CPPADCG_ASSERT_KNOWN(indexes.size() >= 1, "Invalid number of indexes");
                std::string indexExpr;
                const Plane2DIndexPattern& pip = static_cast<const Plane2DIndexPattern&> (ip);
                if (pip.getPattern1() != NULL) {
                    indexExpr += indexPattern2String(*pip.getPattern1(), *indexes[0]);
                }

                if (pip.getPattern2() != NULL) {
                    if (pip.getPattern1() != NULL)
                        indexExpr += " + ";
                    indexExpr += indexPattern2String(*pip.getPattern2(), *indexes.back());
                }

                return indexExpr;
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
    inline std::string CLanguage<Base>::linearIndexPattern2String(const LinearIndexPattern& lip,
                                                                 const IndexDclrOperationNode<Base>& index) {
        long dy = lip.getLinearSlopeDy();
        long dx = lip.getLinearSlopeDx();
        long b = lip.getLinearConstantTerm();
        long xOffset = lip.getXOffset();

        std::stringstream ss;
        if (dy != 0) {
            if (xOffset != 0) {
                ss << "(";
            }
            ss << (*index.getName());
            if (xOffset != 0) {
                ss << " - " << xOffset << ")";
            }

            if (dx != 1) {
                ss << " / " << dx;
            }
            if (dy != 1) {
                ss << " * " << dy;
            }
        } else if (b == 0) {
            ss << "0"; // when dy == 0 and b == 0
        }

        if (b != 0) {
            if (dy != 0)
                ss << " + ";
            ss << b;
        }
        return ss.str();
    }

}

#endif