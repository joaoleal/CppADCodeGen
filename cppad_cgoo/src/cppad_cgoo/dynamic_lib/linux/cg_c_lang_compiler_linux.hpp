#ifndef CPPAD_CG_C_LANG_COMPILER_LINUX_INCLUDED
#define	CPPAD_CG_C_LANG_COMPILER_LINUX_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#ifdef __linux__
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>

namespace CppAD {

    template<class Base>
    void CLangCompiler<Base>::createFolder(const std::string& folder) {
        int ret = mkdir(folder.c_str(), 0755);
        if (ret == -1) {
            if (errno != EEXIST) {
                const char* error = strerror(errno);
                CPPADCG_ASSERT_KNOWN(false, error);
            }
        }
    }

    template<class Base>
    std::string CLangCompiler<Base>::createPath(const std::string& baseFolder, const std::string& file) {
        return baseFolder + "/" + file;
    }

    template<class Base>
    void CLangCompiler<Base>::callExecutable(const std::string& executable,
                                             const std::vector<std::string>& args,
                                             bool createPipe,
                                             const std::string& pipeMessage) {

        int fd[2];

        if (createPipe) {
            //Create pipe for piping source to the compiler
            CPPADCG_ASSERT_KNOWN(pipe(fd) == 0, "Failed to create pipe");
        }

        //Fork the compiler, pipe source to it, wait for the compiler to exit
        pid_t pid = fork();
        CPPADCG_ASSERT_KNOWN(pid >= 0, "Failed to fork program");

        if (pid == 0) {
            //  Child process

            if (createPipe) {
                // close write end of pipe
                close(fd[1]);
                // Send pipe input to stdin
                close(STDIN_FILENO);
                dup2(fd[0], STDIN_FILENO);
            }

            std::vector<char*> args2(args.size() + 1);
            for (size_t i = 0; i < args.size(); i++) {
                const size_t s = args[i].size() + 1;
                args2[i] = new char[s];
                for (size_t c = 0; c < s - 1; c++) {
                    args2[i][c] = args[i].at(c);
                }
                args2[i][s - 1] = '\0';
            }
            args2[args.size()] = (char *) NULL; // END             

            execv(executable.c_str(), &args2[0]);

            for (size_t i = 0; i < args.size(); i++) {
                delete [] args2[i];
            }

            exit(0);
        }

        // Parent process
        if (createPipe) {
            // close read end of pipe
            close(fd[0]);
            //Pipe source to the executable
            write(fd[1], pipeMessage.c_str(), pipeMessage.size());
            close(fd[1]);
        }

        //Wait for the executable to exit
        int status;
        CPPADCG_ASSERT_KNOWN(wait(&status) >= 0, ("Failed while waiting for '" + executable + "'").c_str());
    }

    template<class Base>
    std::string CLangCompiler<Base>::escapePath(const std::string& path) {
        return std::string("\"") + path + "\"";
    }
}
#endif

#endif
