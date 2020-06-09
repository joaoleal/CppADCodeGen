#ifndef CPPAD_CG_SYSTEM_INCLUDED
#define CPPAD_CG_SYSTEM_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *    Copyright (C) 2020 Joao Leal
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

/**
 * System dependent functions
 */
namespace CppAD::cg::system {

template<class T = int >
class SystemInfo {
public:
    static const std::string DYNAMIC_LIB_EXTENSION;
    static const std::string STATIC_LIB_EXTENSION;
};

/**
 * Calls an external executable (system dependent).
 * In the case of an error during execution an exception will be thrown.
 * 
 * @param executable the executable path
 * @param args the command line arguments to the executable
 * @param stdOutErrMessage standard output and standard error message
 *                         from the executable
 * @param stdInMessage information to pass as standard input to the executable
 * @throws CGException on failure to call the executable
 */
inline void callExecutable(const std::filesystem::path& executable,
                           const std::vector<std::string>& args,
                           std::string* stdOutErrMessage = nullptr,
                           const std::string* stdInMessage = nullptr);

} // END namespace

#endif