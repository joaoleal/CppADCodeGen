#ifndef CPPAD_CG_C_LANGUAGE_INDEX_PATTERNS_INCLUDED
#define CPPAD_CG_C_LANGUAGE_INDEX_PATTERNS_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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

    template<class Base>
    inline void CLanguage<Base>::generateNames4RandomIndexPatterns(const std::set<RandomIndexPattern*>& randomPatterns) {
        std::ostringstream os;

        std::set<std::string> usedNames;

        // save existing names so that they are not to overridden 
        // (independent variable names might have already used them)
        std::set<RandomIndexPattern*>::const_iterator it;
        for (it = randomPatterns.begin(); it != randomPatterns.end(); ++it) {
            RandomIndexPattern* ip = *it;
            if (!ip->getName().empty()) {
                usedNames.insert(ip->getName());
            }
        }

        // create new names for the index pattern arrays without a name
        size_t c = 0;
        for (it = randomPatterns.begin(); it != randomPatterns.end(); ++it) {
            RandomIndexPattern* ip = *it;
            if (ip->getName().empty()) {
                // new name required
                std::string name;
                do {
                    os << _C_STATIC_INDEX_ARRAY << c;
                    name = os.str();
                    os.str("");
                    c++;
                } while (usedNames.find(name) != usedNames.end());

                ip->setName(name);
            }
        }

    }

    template<class Base>
    inline void CLanguage<Base>::printRandomIndexPatternDeclaration(std::ostringstream& os,
                                                                    const std::string& identation,
                                                                    const std::set<RandomIndexPattern*>& randomPatterns) {
        typename std::set<RandomIndexPattern*>::const_iterator itr;
        for (itr = randomPatterns.begin(); itr != randomPatterns.end(); ++itr) {
            RandomIndexPattern* ip = *itr;
            if (ip->getType() == RANDOM1D) {
                /**
                 * 1D
                 */
                Random1DIndexPattern* ip1 = static_cast<Random1DIndexPattern*> (ip);
                const std::map<size_t, size_t>& x2y = ip1->getValues();

                std::vector<size_t> y(x2y.rbegin()->first + 1);
                std::map<size_t, size_t>::const_iterator it;
                for (it = x2y.begin(); it != x2y.end(); ++it)
                    y[it->first] = it->second;

                os << identation;
                printStaticIndexArray(os, ip->getName(), y);
            } else {
                assert(ip->getType() == RANDOM2D);
                /**
                 * 2D
                 */
                Random2DIndexPattern* ip2 = static_cast<Random2DIndexPattern*> (ip);
                os << identation;
                printStaticIndexMatrix(os, ip->getName(), ip2->getValues());
            }
        }
    }

    template<class Base>
    inline void CLanguage<Base>::createIndexDeclaration() {
        if (_indexes->empty())
            return;

        std::set<const IndexDclrOperationNode<Base>*> funcArgs(_funcArgIndexes.begin(), _funcArgIndexes.end());

        bool first = true;

        printRandomIndexPatternDeclaration(_ss, _spaces, *_indexRandomPatterns);

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
    void CLanguage<Base>::printStaticIndexArray(std::ostringstream& os,
                                                const std::string& name,
                                                const std::vector<size_t>& values) {
        os << "static unsigned long const " << name << "[" << values.size() << "] = {";
        if (!values.empty()) {
            os << values[0];
            for (size_t i = 1; i < values.size(); i++) {
                os << "," << values[i];
            }
        }
        os << "};\n";
    }

    template<class Base>
    void CLanguage<Base>::printStaticIndexMatrix(std::ostringstream& os,
                                                 const std::string& name,
                                                 const std::map<size_t, std::map<size_t, size_t> >& values) {
        size_t m = 0;
        size_t n = 0;

        std::map<size_t, std::map<size_t, size_t> >::const_iterator it;
        std::map<size_t, size_t>::const_iterator ity2z;

        if (!values.empty()) {
            m = values.rbegin()->first + 1;

            for (it = values.begin(); it != values.end(); ++it) {
                if (!it->second.empty())
                    n = std::max(n, it->second.rbegin()->first + 1);
            }
        }

        os << "static unsigned long const " << name << "[" << m << "][" << n << "] = {";
        size_t x = 0;
        for (it = values.begin(); it != values.end(); ++it) {
            if (it->first != x) {
                while (it->first != x) {
                    os << "{},";
                    x++;
                }
            }

            os << "{";

            size_t y = 0;
            for (ity2z = it->second.begin(); ity2z != it->second.end(); ++ity2z) {
                if (ity2z->first != y) {
                    while (ity2z->first != y) {
                        os << "0,";
                        y++;
                    }
                }

                os << ity2z->second;
                if (ity2z->first != it->second.rbegin()->first) os << ",";

                y++;
            }

            os << "}";
            if (it->first != values.rbegin()->first) os << ",";

            x++;
        }
        os << "};\n";
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
            case RANDOM1D:
            {
                CPPADCG_ASSERT_KNOWN(indexes.size() == 1, "Invalid number of indexes");
                const Random1DIndexPattern& rip = static_cast<const Random1DIndexPattern&> (ip);
                CPPADCG_ASSERT_KNOWN(!rip.getName().empty(), "Invalid name for array");
                return rip.getName() + "[" + (*indexes[0]->getName()) + "]";
            }
            case RANDOM2D:
            {
                CPPADCG_ASSERT_KNOWN(indexes.size() == 2, "Invalid number of indexes");
                const Random2DIndexPattern& rip = static_cast<const Random2DIndexPattern&> (ip);
                CPPADCG_ASSERT_KNOWN(!rip.getName().empty(), "Invalid name for array");
                return rip.getName() + "[" + (*indexes[0]->getName()) + "][" + (*indexes[1]->getName()) + "]";
            }
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