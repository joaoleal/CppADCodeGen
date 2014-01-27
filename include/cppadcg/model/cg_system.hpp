#ifndef CPPAD_CG_SYSTEM_INCLUDED
#define CPPAD_CG_SYSTEM_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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

    /**
     * System dependent functions
     */
    namespace system {

        template<class T = int >
        class SystemInfo {
        public:
            static const std::string DYNAMIC_LIB_EXTENSION;
            static const std::string STATIC_LIB_EXTENSION;
        };

        /**
         * creates a new folder (system dependent)
         * 
         * @param folder the path to the folder
         */
        inline void createFolder(const std::string& folder);

        /**
         * Creates a new path (system dependent)
         * 
         * @param baseFolder the path to the base folder
         * @param file the file or folder name inside the base folder
         * @return the new path
         */
        inline std::string createPath(const std::string& baseFolder, const std::string& file);

        /**
         * Escapes a file or folder path (system dependent)
         * 
         * @param path the file/folder path
         * @return the escaped file/folder path
         */
        inline std::string escapePath(const std::string& path);

        inline std::string filenameFromPath(const std::string& path);

        /**
         * Calls an external executable (system dependent)
         * 
         * @param executable the executable path
         * @param args the command line arguments to the executable
         * @param pipe whether or not to create a pipe to the executable
         * @param message the information to pass in the pipe
         */
        inline void callExecutable(const std::string& executable,
                                   const std::vector<std::string>& args,
                                   bool pipe = false,
                                   const std::string& message = "") throw (CGException);

        inline double currentTime();
    }
}

#endif
