#ifndef CPPAD_CG_LANGUAGE_C_ARRAYS_INCLUDED
#define CPPAD_CG_LANGUAGE_C_ARRAYS_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2014 Ciengis
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
namespace cg {

template<class Base>
void LanguageC<Base>::printArrayCreationOp(OperationNode<Base>& array) {
    CPPADCG_ASSERT_KNOWN(array.getArguments().size() > 0, "Invalid number of arguments for array creation operation");
    const size_t id = array.getVariableID();
    const std::vector<Argument<Base> >& args = array.getArguments();
    const size_t argSize = args.size();

    size_t startPos = id - 1;

    bool firstElement = true;
    for (size_t i = 0; i < argSize; i++) {
        bool newValue = !isSameArgument(args[i], _tmpArrayValues[startPos + i]);

        if (newValue) {
            if (firstElement) {
                _code << _indentation << auxArrayName_ << " = " << _nameGen->generateTemporaryArray(array) << "; // size: " << args.size() << "\n";
                firstElement = false;
            }

            // try to use a loop for element assignment
            size_t newI = printArrayCreationUsingLoop(startPos, array, i, _tmpArrayValues);

            if (newI == i) {
                // individual element assignment
                _code << _indentation << auxArrayName_ << "[" << i << "] = ";
                print(args[i]);
                _code << ";\n";

                _tmpArrayValues[startPos + i] = &args[i];

            } else {
                i = newI - 1;
            }
        }
    }
}

template<class Base>
void LanguageC<Base>::printSparseArrayCreationOp(OperationNode<Base>& array) {
    const std::vector<size_t>& info = array.getInfo();
    CPPADCG_ASSERT_KNOWN(info.size() > 0, "Invalid number of information elements for sparse array creation operation");

    const std::vector<Argument<Base> >& args = array.getArguments();
    const size_t argSize = args.size();

    CPPADCG_ASSERT_KNOWN(info.size() == argSize + 1, "Invalid number of arguments for sparse array creation operation");

    if (argSize == 0)
        return; // empty array

    const size_t id = array.getVariableID();
    size_t startPos = id - 1;

    bool firstElement = true;
    for (size_t i = 0; i < argSize; i++) {
        bool newValue = !isSameArgument(args[i], _tmpSparseArrayValues[startPos + i]);

        if (newValue) {
            if (firstElement) {
                _code << _indentation << auxArrayName_ << " = " << _nameGen->generateTemporarySparseArray(array)
                        << "; // nnz: " << args.size() << "  size:" << info[0] << "\n";
                firstElement = false;
            }

            // try to use a loop for element assignment
            size_t newI = printArrayCreationUsingLoop(startPos, array, i, _tmpSparseArrayValues);

            if (newI == i) {
                // individual element assignment
                _code << _indentation << auxArrayName_ << "[" << i << "] = ";
                print(args[i]);
                _code << "; ";
                // print indexes (location of values)
                _code << _C_SPARSE_INDEX_ARRAY << "[";
                if (startPos != 0) _code << startPos << "+";
                _code << i << "] = " << info[i + 1] << ";\n";

                _tmpSparseArrayValues[startPos + i] = &args[i];

            } else {
                // print indexes (location of values)
                for (size_t j = i; j < newI; j++) {
                    _code << _indentation << _C_SPARSE_INDEX_ARRAY << "[";
                    if (startPos != 0) _code << startPos << "+";
                    _code << j << "] = " << info[j + 1] << ";\n";
                }

                i = newI - 1;
            }


        } else {
            // print indexes (location of values)
            _code << _indentation
                    << _C_SPARSE_INDEX_ARRAY << "[";
            if (startPos != 0) _code << startPos << "+";
            _code << i << "] = " << info[i + 1] << ";\n";
        }

    }
}

template<class Base>
inline size_t LanguageC<Base>::printArrayCreationUsingLoop(size_t startPos,
                                                           OperationNode<Base>& array,
                                                           size_t starti,
                                                           std::vector<const Argument<Base>*>& tmpArrayValues) {
    const std::vector<Argument<Base> >& args = array.getArguments();
    const size_t argSize = args.size();
    size_t i = starti + 1;

    std::ostringstream arrayAssign;

    const Argument<Base>& ref = args[starti];
    if (ref.getOperation() != nullptr) {
        // 
        const OperationNode<Base>& refOp = *ref.getOperation();
        CGOpCode op = refOp.getOperationType();
        if (op == CGOpCode::Inv) {
            /**
             * from independents array
             */
            for (; i < argSize; i++) {
                if (isSameArgument(args[i], tmpArrayValues[startPos + i]))
                    break; // no assignment needed

                if (args[i].getOperation() == nullptr ||
                        args[i].getOperation()->getOperationType() != CGOpCode::Inv ||
                        !_nameGen->isConsecutiveInIndepArray(*args[i - 1].getOperation(), *args[i].getOperation())) {
                    break;
                }
            }

            if (i - starti < 3)
                return starti;

            // use loop
            const std::string& indep = _nameGen->getIndependentArrayName(refOp);
            size_t start = _nameGen->getIndependentArrayIndex(refOp);
            long offset = long(start) - starti;
            if (offset == 0)
                arrayAssign << indep << "[i]";
            else
                arrayAssign << indep << "[" << offset << " + i]";

        } else if (op == CGOpCode::LoopIndexedIndep) {
            /**
             * from independents array in a loop
             */
            size_t pos = refOp.getInfo()[1];
            IndexPattern* refIp = (*_info)->loopIndependentIndexPatterns[pos];

            LinearIndexPattern* refLIp = nullptr;
            SectionedIndexPattern* refSecp = nullptr;

            if (refIp->getType() == IndexPatternType::Linear) {
                refLIp = static_cast<LinearIndexPattern*> (refIp);
            } else if (refIp->getType() == IndexPatternType::Sectioned) {
                refSecp = static_cast<SectionedIndexPattern*> (refIp);
            } else {
                return starti; // cannot determine consecutive elements
            }

            for (; i < argSize; i++) {
                if (isSameArgument(args[i], tmpArrayValues[startPos + i]))
                    break; // no assignment needed

                if (args[i].getOperation() == nullptr ||
                        args[i].getOperation()->getOperationType() != CGOpCode::LoopIndexedIndep) {
                    break; // not an independent index pattern
                }

                if (!_nameGen->isInSameIndependentArray(refOp, *args[i].getOperation()))
                    break;

                pos = args[i].getOperation()->getInfo()[1];
                const IndexPattern* ip = (*_info)->loopIndependentIndexPatterns[pos];

                if (!isOffsetBy(ip, refIp, long(i) - long(starti))) {
                    break; // different pattern type
                }
            }

            if (i - starti < 3)
                return starti;

            std::unique_ptr<Plane2DIndexPattern> p2dip;
            if (refLIp != nullptr) {
                p2dip.reset(encapsulateIndexPattern(*refLIp, starti));
            } else {
                assert(refSecp != nullptr);
                p2dip.reset(encapsulateIndexPattern(*refSecp, starti));
            }

            IndexDclrOperationNode<Base> indexI("i");
            IndexOperationNode<Base> iterationIndexOp(indexI);

            OperationNode<Base> op2(CGOpCode::LoopIndexedIndep, refOp.getInfo(), refOp.getArguments()); //clone
            op2.getInfo()[1] = std::numeric_limits<size_t>::max(); // just to be safe (this would be the index pattern id in the handler)
            op2.getArguments().push_back(iterationIndexOp);

            arrayAssign << _nameGen->generateIndexedIndependent(op2, *p2dip);
        } else if (refOp.getVariableID() >= this->_minTemporaryVarID && op != CGOpCode::LoopIndexedDep && op != CGOpCode::LoopIndexedTmp && op != CGOpCode::Tmp) {
            /**
             * from temporary variable array
             */
            for (; i < argSize; i++) {
                if (isSameArgument(args[i], tmpArrayValues[startPos + i]))
                    break; // no assignment needed
                else if (args[i].getOperation() == nullptr)
                    break;

                const OperationNode<Base>& opNode2 = *args[i].getOperation();
                if (opNode2.getVariableID() < this->_minTemporaryVarID)
                    break;

                CGOpCode op2 = opNode2.getOperationType();
                if (op2 == CGOpCode::LoopIndexedIndep || op2 == CGOpCode::LoopIndexedDep || op2 == CGOpCode::LoopIndexedTmp || op2 == CGOpCode::Tmp)
                    break;

                if (!_nameGen->isConsecutiveInTemporaryVarArray(*args[i - 1].getOperation(), *args[i].getOperation()))
                    break;
            }

            if (i - starti < 3)
                return starti;

            // use loop
            const std::string& tmpName = _nameGen->getTemporaryVarArrayName(refOp);
            size_t start = _nameGen->getTemporaryVarArrayIndex(refOp);
            long offset = long(start) - starti;
            if (offset == 0)
                arrayAssign << tmpName << "[i]";
            else
                arrayAssign << tmpName << "[" << offset << " + i]";

        } else {
            // no loop used
            return starti;
        }
    } else {
        /**
         * constant value?
         */
        const Base& value = *args[0].getParameter();
        for (; i < argSize; i++) {
            if (args[i].getParameter() == nullptr || *args[i].getParameter() != value) {
                break; // not the same constant value
            }

            const Argument<Base>* oldArg = tmpArrayValues[startPos + i];
            if (oldArg != nullptr && oldArg->getParameter() != nullptr && *oldArg->getParameter() == value) {
                break; // values are the same (no need to redefine)
            }
        }

        if (i - starti < 3)
            return starti;

        arrayAssign << value;
    }

    /**
     * print the loop
     */
    _code << _indentation << "for(i = " << starti << "; i < " << i << "; i++) "
            << auxArrayName_ << "[i] = " << arrayAssign.str() << ";\n";

    /**
     * update values in the global temporary array
     */
    for (size_t ii = starti; ii < i; ii++) {
        tmpArrayValues[startPos + ii] = &args[ii];
    }

    return i;
}

template<class Base>
inline std::string LanguageC<Base>::getTempArrayName(const OperationNode<Base>& op) {
    if (op.getOperationType() == CGOpCode::ArrayCreation)
        return _nameGen->generateTemporaryArray(op);
    else
        return _nameGen->generateTemporarySparseArray(op);
}

template<class Base>
void LanguageC<Base>::printArrayElementOp(OperationNode<Base>& op) {
    CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 2, "Invalid number of arguments for array element operation");
    CPPADCG_ASSERT_KNOWN(op.getArguments()[0].getOperation() != nullptr, "Invalid argument for array element operation");
    CPPADCG_ASSERT_KNOWN(op.getInfo().size() == 1, "Invalid number of information indexes for array element operation");

    OperationNode<Base>& arrayOp = *op.getArguments()[0].getOperation();
    std::string arrayName;
    if (arrayOp.getOperationType() == CGOpCode::ArrayCreation)
        arrayName = _nameGen->generateTemporaryArray(arrayOp);
    else
        arrayName = _nameGen->generateTemporarySparseArray(arrayOp);

    _code << "(" << arrayName << ")[" << op.getInfo()[0] << "]";
}

template<class Base>
inline void LanguageC<Base>::printArrayStructInit(const std::string& dataArrayName,
                                                  size_t pos,
                                                  const std::vector<OperationNode<Base>*>& arrays,
                                                  size_t k) {
    _ss.str("");
    _ss << dataArrayName << "[" << pos << "]";
    printArrayStructInit(_ss.str(), *arrays[k]);
}

template<class Base>
inline void LanguageC<Base>::printArrayStructInit(const std::string& dataArrayName,
                                                  OperationNode<Base>& array) {
    const std::string& aName = createVariableName(array);

    if (array.getOperationType() == CGOpCode::ArrayCreation) {
        size_t size = array.getArguments().size();
        if (size > 0)
            _code << dataArrayName << ".data = " << aName << "; ";
        else
            _code << dataArrayName << ".data = NULL; ";
        _code << dataArrayName << ".size = " << size << "; "
                << dataArrayName << ".sparse = " << false << ";";
    } else {
        CPPADCG_ASSERT_KNOWN(array.getOperationType() == CGOpCode::SparseArrayCreation, "Invalid node type");
        size_t nnz = array.getArguments().size();
        if (nnz > 0)
            _code << dataArrayName << ".data = " << aName << "; ";
        else
            _code << dataArrayName << ".data = NULL; ";
        _code << dataArrayName << ".size = " << array.getInfo()[0] << "; "
                << dataArrayName << ".sparse = " << true << "; "
                << dataArrayName << ".nnz = " << nnz << "; ";
        if (nnz > 0) {
            size_t id = array.getVariableID();
            _code << dataArrayName << ".idx = &(" << _C_SPARSE_INDEX_ARRAY << "[" << (id - 1) << "]);";
        }
    }
    _code << "\n";
}

template<class Base>
inline void LanguageC<Base>::markArrayChanged(OperationNode<Base>& ty) {
    size_t id = ty.getVariableID();
    size_t tySize = ty.getArguments().size();

    if (ty.getOperationType() == CGOpCode::ArrayCreation) {
        for (size_t i = 0; i < tySize; i++) {
            _tmpArrayValues[id - 1 + i] = nullptr;
        }
    } else {
        for (size_t i = 0; i < tySize; i++) {
            _tmpSparseArrayValues[id - 1 + i] = nullptr;
        }
    }
}

} // END cg namespace
} // END CppAD namespace

#endif