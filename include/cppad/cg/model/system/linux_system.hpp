#ifndef CPPAD_CG_LINUX_SYSTEM_INCLUDED
#define CPPAD_CG_LINUX_SYSTEM_INCLUDED
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

#if CPPAD_CG_SYSTEM_LINUX
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>

namespace CppAD {
namespace cg {

/**
 * Linux system dependent functions
 */
namespace system {

namespace {

/**
 * Utility class used to close automatically file descriptors
 */
class FDHandler {
public:
    int fd;
    bool closed;
public:

    inline FDHandler() : closed(true) {
    }

    inline explicit FDHandler(int fd) : fd(fd), closed(false) {
    }

    inline void close() {
        if (!closed) {
            ::close(fd);
            closed = true;
        }
    }

    inline ~FDHandler() {
        close();
    }
};

/**
 * Utility class for pipes
 */
class PipeHandler {
public:
    FDHandler read;
    FDHandler write;
public:

    inline void create() throw (CGException) {
        int fd[2]; /** file descriptors used to communicate between processes*/
        if (pipe(fd) < 0) {
            throw CGException("Failed to create pipe");
        }
        read.fd = fd[0];
        read.closed = false;
        write.fd = fd[1];
        write.closed = false;
    }
};

}

template<class T>
const std::string SystemInfo<T>::DYNAMIC_LIB_EXTENSION = ".so";

template<class T>
const std::string SystemInfo<T>::STATIC_LIB_EXTENSION = ".a";

inline std::string getWorkingDirectory() throw (CGException) {
    char buffer[1024];

    char* ret = getcwd(buffer, 1024);
    if (ret == nullptr) {
        const char* error = strerror(errno);
        throw CGException("Failed to get current working directory: ", error);
    }

    return buffer;
}

inline void createFolder(const std::string& folder) throw (CGException) {
    int ret = mkdir(folder.c_str(), 0755);
    if (ret == -1) {
        if (errno != EEXIST) {
            const char* error = strerror(errno);
            throw CGException("Failed to create directory '", folder + "': ", error);
        }
    }
}

inline std::string createPath(const std::string& baseFolder,
                              const std::string& file) {
    return baseFolder + "/" + file;
}

inline std::string escapePath(const std::string& path) {
    return std::string("\"") + path + "\"";
}

inline std::string filenameFromPath(const std::string& path) {
    size_t pos = path.rfind('/');
    if (pos != std::string::npos) {
        if (pos == path.size() - 1) {
            return "";
        } else {
            return path.substr(pos + 1);
        }
    } else {
        return path;
    }
}

inline void callExecutable(const std::string& executable,
                           const std::vector<std::string>& args,
                           bool pipeSource,
                           const std::string& pipeMessage) throw (CGException) {
    std::string execName = filenameFromPath(executable);

    PipeHandler pipeMsg; // file descriptors used to communicate between processes
    pipeMsg.create();

    PipeHandler pipeSrc;
    if (pipeSource) {
        //Create pipe for piping source to the compiler
        pipeSrc.create();
    }

    //Fork the compiler, pipe source to it, wait for the compiler to exit
    pid_t pid = fork();
    if (pid < 0) {
        throw CGException("Failed to fork process");
    }

    if (pid == 0) {
        /***********************************************************************
         * Child process
         **********************************************************************/
        pipeMsg.read.close();

        if (pipeSource) {
            pipeSrc.write.close(); // close write end of pipe
            // Send pipe input to stdin
            close(STDIN_FILENO);
            dup2(pipeSrc.read.fd, STDIN_FILENO);
        }

        auto toCharArray = [](const std::string & args) {
            const size_t s = args.size() + 1;
            char* args2 = new char[s];
            for (size_t c = 0; c < s - 1; c++) {
                args2[c] = args.at(c);
            }
            args2[s - 1] = '\0';
            return args2;
        };

        std::vector<char*> args2(args.size() + 2);
        args2[0] = toCharArray(execName);
        for (size_t i = 0; i < args.size(); i++) {
            args2[i + 1] = toCharArray(args[i]);
        }
        args2.back() = (char *) nullptr; // END             

        int eCode = execv(executable.c_str(), &args2[0]);

        for (size_t i = 0; i < args.size(); i++) {
            delete [] args2[i];
        }

        if (eCode < 0) {
            char buf[512];
            std::string error = executable + ": " + strerror_r(errno, buf, 511); // thread safe
            ssize_t size = error.size() + 1;
            if (write(pipeMsg.write.fd, error.c_str(), size) != size) {
                std::cerr << "Failed to send message to parent process" << std::endl;
            }
            std::cerr << "*** ERROR: exec failed" << std::endl;
            exit(EXIT_FAILURE);
        }

        exit(EXIT_SUCCESS);
    }

    /***************************************************************************
     * Parent process
     **************************************************************************/
    pipeMsg.write.close();

    auto readCErrorMsg = []() {
        int error = errno;
        errno = 0;
        char buf[512];
        return std::string(strerror_r(error, buf, 512));
    };

    std::string writeError;
    if (pipeSource) {
        // close read end of pipe
        pipeSrc.read.close();
        //Pipe source to the executable
        ssize_t writeFlag = write(pipeSrc.write.fd, pipeMessage.c_str(), pipeMessage.size());
        if (writeFlag == -1)
            writeError = readCErrorMsg() + " ";
        pipeSrc.write.close();
    }

    //Wait for the executable to exit
    int status;
    // Read message from the child
    std::string message;
    char buffer[128];
    do {
        while (read(pipeMsg.read.fd, buffer, sizeof (buffer)) > 0) {
            message += buffer;
            if (message.size() > 1e4) break;
        }

        if (waitpid(pid, &status, 0) < 0) {
            throw CGException("Waitpid failed for pid ", pid, " [", readCErrorMsg(), "]");
        }
    } while (!WIFEXITED(status) && !WIFSIGNALED(status));

    pipeMsg.read.close();

    if (!writeError.empty()) {
        std::stringstream s;
        s << "Failed to write to pipe";
        if (!message.empty()) s << ": " << message;
        else s << ": " << writeError;
        throw CGException(s.str());
    }

    if (WIFEXITED(status)) {
        if (WEXITSTATUS(status) != EXIT_SUCCESS) {
            std::stringstream s;
            s << "Executable '" << executable << "' (pid " << pid << ") exited with code " << WEXITSTATUS(status);
            if (!message.empty()) s << ": " << message;
            throw CGException(s.str());
        }
    } else if (WIFSIGNALED(status)) {
        std::stringstream s;
        s << "Executable '" << executable << "' (pid " << pid << ") terminated by signal " << WTERMSIG(status);
        if (!message.empty()) s << ": " << message;
        throw CGException(s.str());
    }
}

} // END system namespace

} // END cg namespace
} // END CppAD namespace

#endif
#endif