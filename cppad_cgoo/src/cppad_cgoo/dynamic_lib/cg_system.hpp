#ifndef CPPAD_CG_SYSTEM_INCLUDED
#define	CPPAD_CG_SYSTEM_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    /**
     * System dependent functions
     */
    namespace system {
        
        /**
         * creates a new folder (system dependent)
         * 
         * \param folder the path to the folder
         */
        inline void createFolder(const std::string& folder);

        /**
         * Creates a new path (system dependent)
         * 
         * \param baseFolder the path to the base folder
         * \param file the file or folder name inside the base folder
         * \return the new path
         */
        inline std::string createPath(const std::string& baseFolder, const std::string& file);

        /**
         * Calls an external executable (system dependent)
         * 
         * \param executable the executable path
         * \param args the command line arguments to the executable
         * \param pipe whether or not to create a pipe to the executable
         * \param message the information to pass in the pipe
         */
        inline void callExecutable(const std::string& executable,
                            const std::vector<std::string>& args,
                            bool pipe = false,
                            const std::string& message = "");
        /**
         * Escapes a file or folder path (system dependent)
         * 
         * \param path the file/folder path
         * \return the escaped file/folder path
         */
        inline std::string escapePath(const std::string& path);

        inline double currentTime();
    }
}

#endif
