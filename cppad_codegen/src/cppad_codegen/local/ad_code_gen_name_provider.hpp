/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#ifndef CPPAD_CODEGEN_AD_CODE_GEN_NAME_PROVIDER_INCLUDED
#define	CPPAD_CODEGEN_AD_CODE_GEN_NAME_PROVIDER_INCLUDED

#include <vector>
#include <map>
#include <string>

CPPAD_BEGIN_NAMESPACE

typedef struct VarID {
    size_t order;
    addr_t taddr;
} VarID;

template<class Base>
class CodeGenNameProvider {
public:
    virtual std::string generateVarName(size_t d, addr_t var) = 0;

    virtual std::string generatePartialName(size_t d, addr_t var) = 0;

    std::string generateVarName(size_t d, const AD<Base>& var) {
        return generateVarName(d, getTapeAddr(var));
    }

    virtual void clearUsedVariables() = 0;

    virtual void clearUsedPartials() = 0;

    virtual const std::string& baseTypeName() const = 0;

    virtual const std::string& zero() const = 0;

    virtual const std::string& one() const = 0;

    virtual const std::string& endl() const = 0;

    virtual const std::string& compareChangeCounter() const = 0;

    virtual const std::string& tempIntegerVarName() const = 0;

    virtual const std::string& tempBaseVarName() const = 0;

    virtual std::vector<VarID> getUsedVariables() const = 0;

    virtual std::vector<VarID> getUsedPartials() const = 0;

    virtual const std::string& PrintStreamName() const = 0;

    virtual const std::string& StartComment() const = 0;

    virtual const std::string& EndComment() const = 0;

    virtual std::string PrintBase(const Base& value) const = 0;

    virtual std::string CastToBase(const int value) const = 0;

    template <typename T>
    std::string toString(T v) const {
        std::stringstream str;
        // make sure all digits of floating point values are printed
        str << std::setprecision(99) << v;
        std::string result;
        str >> result;
        return result;
    }

protected:

    addr_t getTapeAddr(const AD<Base>& var) const {
        return var.taddr_;
    }
};

//template <>
//std::string ADCodeGenNameProvider::toString<double>(double v) const {
//    std::stringstream str;
//    str << std::setprecision(std::numeric_limits<double>::digits10 + 2)
//            << v;
//    std::string result;
//    str >> result;
//    return result;
//}
//
//template <>
//std::string ADCodeGenNameProvider::toString<float>(float v) const {
//    std::stringstream str;
//    str << std::setprecision(std::numeric_limits<float>::digits10 + 2)
//            << v;
//    std::string result;
//    str >> result;
//    return result;
//}

template<class Base>
class DefaultCCodeGenNameProvider : public CodeGenNameProvider<Base> {
protected:
    /// base variable type name for code generation
    std::string _typeName;
    std::string _endl;
    std::string _integerVarName; // temporary integer variable name
    std::string _baseVarName; // temporary base variable name
    std::string _zero; // value zero
    std::string _one; // value 1
    std::string _compCounter;
    std::string _printStreamName;
    std::string _startComment;
    std::string _endComment;
    std::map<addr_t, std::vector<bool> > _usedVariables;
    std::map<addr_t, std::vector<bool> > _usedPartials;
public:
    /// default constructor

    DefaultCCodeGenNameProvider() {
        _typeName = "double";
        _endl = ";\n";
        _integerVarName = "k";
        _baseVarName = "v_aux";
        _zero = this->PrintBase(Base(0));
        _one = this->PrintBase(Base(1));
        _compCounter = "comparison";
        _printStreamName = "stdout";
        _startComment = "/**";
        _endComment = "*/";
    }

    /// destructor

    virtual ~DefaultCCodeGenNameProvider() {
    }

    virtual std::string generateVarName(size_t d, addr_t var) {
        std::map<addr_t, std::vector<bool> >::iterator it = _usedVariables.find(var);
        if (it == _usedVariables.end()) {
            _usedVariables[var].resize(d + 1);
        } else if (it->second.size() < d + 1) {
            it->second.resize(d + 1);
        }
        _usedVariables[var][d] = true;

        std::string name = "v_" + this->toString(d) + "_" + this->toString(var);
        return name;
    }

    virtual std::string generatePartialName(size_t d, addr_t var) {
        std::map<addr_t, std::vector<bool> >::iterator it = _usedPartials.find(var);
        if (it == _usedPartials.end()) {
            _usedPartials[var].resize(d + 1);
        } else if (it->second.size() < d + 1) {
            it->second.resize(d + 1);
        }
        _usedPartials[var][d] = true;

        std::string name = "p_" + this->toString(d) + "_" + this->toString(var);
        return name;
    }

    virtual std::vector<VarID> getUsedVariables() const {
        std::vector<VarID> used;
        std::map<addr_t, std::vector<bool> >::const_iterator it;
        for (it = _usedVariables.begin(); it != _usedVariables.end(); ++it) {
            const std::vector<bool>& v = it->second;
            for (size_t i = 0; i < v.size(); i++) {
                if (v[i]) {
                    VarID id;
                    id.order = i;
                    id.taddr = it->first;
                    used.push_back(id);
                }
            }
        }

        return used;
    }

    virtual std::vector<VarID> getUsedPartials() const {
        std::vector<VarID> used;
        std::map<addr_t, std::vector<bool> >::const_iterator it;
        for (it = _usedPartials.begin(); it != _usedPartials.end(); ++it) {
            const std::vector<bool>& v = it->second;
            for (size_t i = 0; i < v.size(); i++) {
                if (v[i]) {
                    VarID id;
                    id.order = i;
                    id.taddr = it->first;
                    used.push_back(id);
                }
            }
        }

        return used;
    }

    virtual void clearUsedVariables() {
        _usedVariables.clear();
    }

    virtual void clearUsedPartials() {
        _usedPartials.clear();
    }

    virtual const std::string& baseTypeName() const {
        return _typeName;
    }

    virtual const std::string& zero() const {
        return _zero;
    }

    virtual const std::string& one() const {
        return _one;
    }

    virtual const std::string& endl() const {
        return _endl;
    }

    virtual const std::string& tempIntegerVarName() const {
        return _integerVarName;
    }

    virtual const std::string& tempBaseVarName() const {
        return _baseVarName;
    }

    virtual const std::string& compareChangeCounter() const {
        return _compCounter;
    }

    virtual const std::string& PrintStreamName() const {
        return _printStreamName;
    }

    virtual const std::string& StartComment() const {
        return _startComment;
    }

    virtual const std::string& EndComment() const {
        return _endComment;
    }

    virtual std::string PrintBase(const Base& value) const {
        std::stringstream str;
        // make sure all digits of floating point values are printed

        str << std::scientific << std::setprecision(std::numeric_limits< Base >::digits10 + 2) << value;
        std::string result;
        str >> result;
        return result;
    }

    virtual std::string CastToBase(const int value) const {
        return CastStringToBase(this->toString(value));
    }

protected:

    virtual std::string CastStringToBase(const std::string& value) const {
        return "((" + this->baseTypeName() + ") " + value + ")";
    }
};

CPPAD_END_NAMESPACE

#endif	/* CPPAD_CODEGEN_AD_CODE_GEN_NAME_PROVIDER_INCLUDED */

