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

    class JobListener {
    public:
        virtual void jobStarted(const std::vector<Job>& job) = 0;

        virtual void jobEndended(const std::vector<Job>& job, double elapsed) = 0;
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
         * number of spaces per indentation level
         */
        size_t _indent;
        /**
         * 
         */
        std::ostringstream _os;
        /**
         * 
         */
        std::string _action;
        /**
         * 
         */
        std::string _actionEnd;
        /**
         * 
         */
        std::set<JobListener*> _listeners;
    public:

        JobTimer() :
            _maxLineWidth(80),
            _indent(2),
            _action("generating"),
            _actionEnd("generated") {
        }

        inline bool isVerbose() const {
            return _verbose;
        }

        inline void setVerbose(bool verbose) {
            _verbose = verbose;
        }

        inline const std::string& getActionName() const {
            return _action;
        }

        inline void setActionName(const std::string& action) {
            _action = action;
        }

        inline const std::string& getActionEndName()const {
            return _actionEnd;
        }

        inline void setActionEndName(const std::string& actionEnd) {
            _actionEnd = actionEnd;
        }

        inline size_t getMaxLineWidth() const {
            return _maxLineWidth;
        }

        inline void setMaxLineWidth(size_t width) {
            _maxLineWidth = width;
        }

        /**
         * Provides the number of currently running jobs
         * 
         * @return the number of running jobs
         */
        inline size_t getJobCount() const {
            return _jobs.size();
        }

        inline void addListener(JobListener& l) {
            _listeners.insert(&l);
        }

        inline bool removeListener(JobListener& l) {
            return _listeners.erase(&l) > 0;
        }

        inline void startingJob(const std::string& jobName,
                                const std::string& prefix = "") {
            if (!_verbose) {
                return;
            }

            _jobs.push_back(Job(jobName));
            Job& job = _jobs.back();

            size_t indent = 0;
            if (_jobs.size() > 1) {
                Job& parent = _jobs[_jobs.size() - 2]; // must be after adding job
                if (!parent._nestedJobs) {
                    parent._nestedJobs = true;
                    std::cout << "\n";
                }
                indent = _indent * (_jobs.size() - 1);
            }

            _os.str("");
            if (indent > 0) _os << std::string(indent, ' ');
            if (!prefix.empty()) _os << prefix << " ";
            _os << _action << " " << job.name() << " ...";

            char f = std::cout.fill();
            std::cout << std::setw(_maxLineWidth) << std::setfill('.') << std::left << _os.str();
            std::cout.flush();
            std::cout.fill(f); // restore fill character

            // notify listeners
            std::set<JobListener*>::const_iterator itl;
            for (itl = _listeners.begin(); itl != _listeners.end(); ++itl) {
                (*itl)->jobStarted(_jobs);
            }
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
                    _os << std::string(_indent * (_jobs.size() - 1), ' ');
                _os << _actionEnd << " " << job.name() << " ...";

                char f = std::cout.fill();
                std::cout << std::setw(_maxLineWidth) << std::setfill('.') << std::left << _os.str();
                std::cout.fill(f); // restore fill character
            }

            std::cout << " done [" << std::fixed << std::setprecision(3) << elapsed << "]" << std::endl;

            // notify listeners
            std::set<JobListener*>::const_iterator itl;
            for (itl = _listeners.begin(); itl != _listeners.end(); ++itl) {
                (*itl)->jobEndended(_jobs, elapsed);
            }

            _jobs.pop_back();
        }

    };
}

#endif