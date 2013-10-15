#ifndef CPPAD_CG_JOB_TIME_INCLUDED
#define CPPAD_CG_JOB_TIME_INCLUDED
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

    /**
     * Utility class used to print elapsed times of jobs
     */
    class JobTime {
    protected:
        /**
         * Whether or not to print progress information to the standard 
         * output
         */
        bool _verbose;
        /**
         * saves the current job names
         */
        std::vector<std::string> _jobNames;
        /**
         * Whether or not there are jobs inside other jobs
         */
        std::vector<bool> _nestedJobs;
        /**
         * auxiliary variable to measure the elapsed time for each job
         */
        std::vector<double> _beginTimes;
        /**
         * 
         */
        size_t _maxLineWidth;
    private:
        /**
         * 
         */
        std::ostringstream _os;
    public:

        JobTime() :
            _maxLineWidth(80) {
        }

        inline bool isVerbose() const {
            return _verbose;
        }

        inline void setVerbose(bool verbose) {
            _verbose = verbose;
        }

        inline void startingJob(const std::string& jobName) {
            if (!_verbose) {
                return;
            }
            
            size_t ident = 0;
            if (!_jobNames.empty()) {
                if (!_nestedJobs.back())
                    std::cout << "\n";
                ident = 3 * _jobNames.size();
            }

            _os.str("");
            if (ident > 0)
                _os << std::string(ident, ' ');
            _os << "generating " << jobName << " ...";

            std::cout << std::setw(_maxLineWidth) << std::setfill('.') << std::left << _os.str();
            std::cout.flush();

            std::fill(_nestedJobs.begin(), _nestedJobs.end(), true);

            _nestedJobs.push_back(false);
            _jobNames.push_back(jobName);
            _beginTimes.push_back(system::currentTime());

        }

        inline void finishedJob() {
            if (!_verbose) {
                return;
            }
            
            assert(_nestedJobs.size() > 0);

            double beginTime = _beginTimes.back();
            double endTime = system::currentTime();
            double elapsed = endTime - beginTime;
            if (_nestedJobs.back()) {
                _os.str("");
                if (!_jobNames.empty())
                    _os << std::string(3 * (_jobNames.size() - 1), ' ');
                _os << "generated " << _jobNames.back() << " ...";

                std::cout << std::setw(_maxLineWidth) << std::setfill('.') << std::left << _os.str();
            }

            std::cout << " done [" << std::fixed << std::setprecision(3) << elapsed << "]" << std::endl;

            _nestedJobs.pop_back();
            _jobNames.pop_back();
            _beginTimes.pop_back();
        }
    };
}

#endif