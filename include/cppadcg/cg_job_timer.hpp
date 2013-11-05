#ifndef CPPAD_CG_JOB_TIMER_INCLUDED
#define CPPAD_CG_JOB_TIMER_INCLUDED
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

    class Job {
    private:
        /**
         * Job name
         */
        std::string _name;
        /**
         * Job starting time (seconds)
         */
        double _beginTime;
        /**
         * Whether or not there are/were other jobs inside
         */
        bool _nestedJobs;
    public:

        inline Job(const std::string& name) :
            _name(name),
            _beginTime(system::currentTime()),
            _nestedJobs(false) {
        }

        inline const std::string& name() const {
            return _name;
        }

        inline double beginTime() const {
            return _beginTime;
        }

        friend class JobTimer;
    };

    /**
     * Utility class used to print elapsed times of jobs
     */
    class JobTimer {
    protected:
        /**
         * Whether or not to print progress information to the standard 
         * output
         */
        bool _verbose;
    private:
        /**
         * saves the current job names
         */
        std::vector<Job> _jobs;
        /**
         * 
         */
        size_t _maxLineWidth;
        /**
         * 
         */
        std::ostringstream _os;
    public:

        JobTimer() :
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

            _jobs.push_back(Job(jobName));
            Job& job = _jobs.back();

            size_t ident = 0;
            if (_jobs.size() > 1) {
                Job& parent = _jobs[_jobs.size() - 2]; // must be after adding job
                if (!parent._nestedJobs) {
                    parent._nestedJobs = true;
                    std::cout << "\n";
                }
                ident = 3 * (_jobs.size() - 1);
            }

            _os.str("");
            if (ident > 0)
                _os << std::string(ident, ' ');
            _os << "generating " << job.name() << " ...";

            char f = std::cout.fill();
            std::cout << std::setw(_maxLineWidth) << std::setfill('.') << std::left << _os.str();
            std::cout.flush();
            std::cout.fill(f); // restore fill character
        }

        inline void finishedJob() {
            if (!_verbose) {
                return;
            }

            assert(_jobs.size() > 0);
            Job& job = _jobs.back();

            double elapsed = system::currentTime() - job.beginTime();
            if (job._nestedJobs) {
                _os.str("");
                if (!_jobs.empty())
                    _os << std::string(3 * (_jobs.size() - 1), ' ');
                _os << "generated " << job.name() << " ...";

                char f = std::cout.fill();
                std::cout << std::setw(_maxLineWidth) << std::setfill('.') << std::left << _os.str();
                std::cout.fill(f); // restore fill character
            }

            std::cout << " done [" << std::fixed << std::setprecision(3) << elapsed << "]" << std::endl;

            _jobs.pop_back();
        }

    };
}

#endif