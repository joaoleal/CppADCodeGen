/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2016 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Johan Hanssen Seferidis
 */

/**
 * This file was adapted by Joao Leal from
 *  https://github.com/Pithikos/C-Thread-Pool/blob/master/thpool.c
 */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <errno.h>
#include <time.h>
#if defined(__linux__)
#include <sys/prctl.h>
#include <time.h>
#include <sys/time.h>
#define __USE_GNU /* required before including  resource.h */
#include <sys/resource.h>
#endif

enum group_strategy {SINGLE_JOB, MULTI_JOB};

typedef struct thpool_* threadpool;
typedef void (* thpool_function_type)(void*);

static volatile threadpool cppadcg_pool;
static volatile int cppadcg_pool_n_threads = 2;
static volatile int cppadcg_pool_disabled = 0; // false
static volatile int cppadcg_pool_verbose = 0; // false
static volatile float cppadcg_pool_multijob_maxgroupwork = 0.75;

static volatile enum group_strategy group_gen_strategy = SINGLE_JOB;

/* ==================== INTERNAL HIGH LEVEL API  ====================== */

static threadpool thpool_init(int num_threads);

static int thpool_add_job(threadpool,
                          thpool_function_type function,
                          void* arg,
                          float* elapsed);

static int thpool_add_jobs(threadpool,
                           thpool_function_type functions[],
                           void* args[],
                           float elapsed[],
                           int order[],
                           int nJobs);

static void thpool_wait(threadpool);

static void thpool_destroy(threadpool);

/* ========================== STRUCTURES ============================ */
/* Binary semaphore */
typedef struct bsem {
    pthread_mutex_t mutex;
    pthread_cond_t   cond;
    int v;
} bsem;


/* Job */
typedef struct job {
    struct job*  prev;                   /* pointer to previous job   */
    thpool_function_type function;       /* function pointer          */
    void*  arg;                          /* function's argument       */
    float* elapsed;                     /* elapsed time              */
} job;

/* Job group */
typedef struct job_group {
    struct job* jobs;                    /* jobs                      */
    int size;                            /* number of jobs            */
} job_group;

/* Job queue */
typedef struct jobqueue {
    pthread_mutex_t rwmutex;             /* used for queue r/w access */
    job  *front;                         /* pointer to front of queue */
    job  *rear;                          /* pointer to rear  of queue */
    bsem *has_jobs;                      /* flag as binary semaphore  */
    int   len;                           /* number of jobs in queue   */
    float total_time;                   /* total expected time to complete the work */
    float highest_expected_return;      /* the time when the last running thread is expected to request new work */
} jobqueue;


/* Thread */
typedef struct thread {
    int       id;                        /* friendly id               */
    pthread_t pthread;                   /* pointer to actual thread  */
    struct thpool_* thpool_p;            /* access to thpool          */
} thread;


/* Threadpool */
typedef struct thpool_ {
    thread**   threads;                  /* pointer to threads        */
    int num_threads;                     /* total number of threads   */
    volatile int num_threads_alive;      /* threads currently alive   */
    volatile int num_threads_working;    /* threads currently working */
    pthread_mutex_t  thcount_lock;       /* used for thread count etc */
    pthread_cond_t  threads_all_idle;    /* signal to thpool_wait     */
    jobqueue*  jobqueue_p;               /* pointer to the job queue  */
    volatile int threads_keepalive;
    volatile int threads_on_hold;
} thpool_;

/* ========================== PUBLIC API ============================ */

void cppadcg_thpool_set_threads(int n) {
    cppadcg_pool_n_threads = n;
}

int cppadcg_thpool_get_threads() {
    return cppadcg_pool_n_threads;
}

void cppadcg_thpool_set_scheduler_strategy(enum group_strategy s) {
    if(cppadcg_pool != NULL) {
        pthread_mutex_lock(&cppadcg_pool->jobqueue_p->rwmutex);
        group_gen_strategy = s;
        pthread_mutex_unlock(&cppadcg_pool->jobqueue_p->rwmutex);
    } else {
        // pool not yet created
        group_gen_strategy = s;
    }
}

enum group_strategy cppadcg_thpool_get_scheduler_strategy() {
    if(cppadcg_pool != NULL) {
        enum group_strategy e;
        pthread_mutex_lock(&cppadcg_pool->jobqueue_p->rwmutex);
        e = group_gen_strategy;
        pthread_mutex_unlock(&cppadcg_pool->jobqueue_p->rwmutex);
        return e;
    } else {
        // pool not yet created
        return group_gen_strategy;
    }
}

void cppadcg_thpool_set_disabled(int disabled) {
    cppadcg_pool_disabled = disabled;
}

int cppadcg_thpool_is_disabled() {
    return cppadcg_pool_disabled;
}

void cppadcg_thpool_set_multijob_maxgroupwork(float v) {
    if(cppadcg_pool != NULL) {
        pthread_mutex_lock(&cppadcg_pool->jobqueue_p->rwmutex);
        cppadcg_pool_multijob_maxgroupwork = v;
        pthread_mutex_unlock(&cppadcg_pool->jobqueue_p->rwmutex);
    } else {
        // pool not yet created
        cppadcg_pool_multijob_maxgroupwork = v;
    }
}

float cppadcg_thpool_get_multijob_maxgroupwork() {
    if(cppadcg_pool != NULL) {
        float r;
        pthread_mutex_lock(&cppadcg_pool->jobqueue_p->rwmutex);
        r = cppadcg_pool_multijob_maxgroupwork;
        pthread_mutex_unlock(&cppadcg_pool->jobqueue_p->rwmutex);
        return r;
    } else {
        // pool not yet created
        return cppadcg_pool_multijob_maxgroupwork;
    }
}

void cppadcg_thpool_set_verbose(int v) {
    cppadcg_pool_verbose = v;
}

int cppadcg_thpool_is_verbose() {
    return cppadcg_pool_verbose;
}

void cppadcg_thpool_prepare() {
    if(cppadcg_pool == NULL) {
        cppadcg_pool = thpool_init(cppadcg_pool_n_threads);
    }
}

void cppadcg_thpool_add_job(thpool_function_type function,
                            void* arg,
                            float* elapsed) {
    if (!cppadcg_pool_disabled) {
        cppadcg_thpool_prepare();
        if (cppadcg_pool != NULL) {
            thpool_add_job(cppadcg_pool, function, arg, elapsed);
            return;
        }
    }

    // thread pool not used
    (*function)(arg);
}

void cppadcg_thpool_add_jobs(thpool_function_type functions[],
                             void* args[],
                             float elapsed[],
                             int order[],
                             int nJobs) {
    int i;
    if (!cppadcg_pool_disabled) {
        cppadcg_thpool_prepare();
        if (cppadcg_pool != NULL) {
            thpool_add_jobs(cppadcg_pool, functions, args, elapsed, order, nJobs);
            return;
        }
    }

    // thread pool not used
    for (i = 0; i < nJobs; ++i) {
        (*functions[i])(args[i]);
    }
}

void cppadcg_thpool_wait() {
    if(cppadcg_pool != NULL) {
        thpool_wait(cppadcg_pool);
    }
}

typedef struct pair_double_int {
    float val;
    int index;
} pair_double_int;

static int comparePair(const void* a, const void* b) {
    if (((pair_double_int*) a)->val < ((pair_double_int*) b)->val)
        return -1;
    if (((pair_double_int*) a)->val == ((pair_double_int*) b)->val)
        return 0;
    return 1;
}

void cppadcg_thpool_update_order(float elapsed[],
                                 int order[],
                                 int nJobs) {
    if(nJobs == 0 || elapsed == NULL || order == NULL)
        return;

    struct pair_double_int elapsed2[nJobs];
    int i;
    int nonZero = 0; // false

    for(i = 0; i < nJobs; ++i) {
        if(elapsed[i] != 0) {
            nonZero = 1;
            break;
        }
    }

    if (nonZero) {
        for(i = 0; i < nJobs; ++i) {
            elapsed2[i].val = elapsed[i];
            elapsed2[i].index = i;
        }

        qsort(elapsed2, nJobs, sizeof(struct pair_double_int), comparePair);

        for (i = 0; i < nJobs; ++i) {
            order[elapsed2[i].index] = nJobs - i - 1; // descending order
        }

        if (cppadcg_pool_verbose) {
            fprintf(stdout, "new order:\n");
            for (i = 0; i < nJobs; ++i) {
                fprintf(stdout, " original: %i   new: %i   time: %e s\n", i, order[i], elapsed[i]);
            }
        }
    } else {
        if (cppadcg_pool_verbose) {
            fprintf(stdout, "order not updated: all times are zero\n");
        }
    }

}

void cppadcg_thpool_shutdown() {
    if(cppadcg_pool != NULL) {
        thpool_destroy(cppadcg_pool);
        cppadcg_pool = NULL;
    }
}

/* ========================== PROTOTYPES ============================ */

static int  thread_init(thpool_* thpool_p,
                        struct thread** thread_p,
                        int id);
static void* thread_do(struct thread* thread_p);
static void  thread_destroy(struct thread* thread_p);

static int   jobqueue_init(thpool_* thpool_p);
static void  jobqueue_clear(thpool_* thpool_p);
static void  jobqueue_push(thpool_* thpool_p,
                           struct job* newjob_p);
static void jobqueue_multipush(thpool_* thpool_p,
                               struct job* newjob[],
                               int nJobs);
static struct job_group* jobqueue_pull(thpool_* thpool_p);
static void  jobqueue_destroy(thpool_* thpool_p);

static void  bsem_init(struct bsem *bsem_p, int value);
static void  bsem_reset(struct bsem *bsem_p);
static void  bsem_post(struct bsem *bsem_p);
static void  bsem_post_all(struct bsem *bsem_p);
static void  bsem_wait(struct bsem *bsem_p);

/* ========================== THREADPOOL ============================ */

/**
 * @brief  Initialize threadpool
 *
 * Initializes a threadpool. This function will not return untill all
 * threads have initialized successfully.
 *
 * @example
 *
 *    ..
 *    threadpool thpool;                     //First we declare a threadpool
 *    thpool = thpool_init(4);               //then we initialize it to 4 threads
 *    ..
 *
 * @param  num_threads   number of threads to be created in the threadpool
 * @return threadpool    created threadpool on success,
 *                       NULL on error
 */
struct thpool_* thpool_init(int num_threads) {
    if (num_threads < 0) {
        num_threads = 0;
    }

    if(cppadcg_pool_verbose) {
        fprintf(stdout, "thpool_init(): Thread pool created with %i threads\n", num_threads);
    }

    if(num_threads == 0) {
        cppadcg_pool_disabled = 1; // true
        return NULL;
    }

    /* Make new thread pool */
    thpool_* thpool_p;
    thpool_p = (struct thpool_*) malloc(sizeof(struct thpool_));
    if (thpool_p == NULL) {
        fprintf(stderr, "thpool_init(): Could not allocate memory for thread pool\n");
        return NULL;
    }
    thpool_p->num_threads = num_threads;
    thpool_p->num_threads_alive = 0;
    thpool_p->num_threads_working = 0;
    thpool_p->threads_on_hold = 0;
    thpool_p->threads_keepalive = 1;

    /* Initialize the job queue */
    if (jobqueue_init(thpool_p) == -1) {
        fprintf(stderr, "thpool_init(): Could not allocate memory for job queue\n");
        free(thpool_p);
        return NULL;
    }

    /* Make threads in pool */
    thpool_p->threads = (struct thread**) malloc(num_threads * sizeof(struct thread*));
    if (thpool_p->threads == NULL) {
        fprintf(stderr, "thpool_init(): Could not allocate memory for threads\n");
        jobqueue_destroy(thpool_p);
        free(thpool_p->jobqueue_p);
        free(thpool_p);
        return NULL;
    }

    pthread_mutex_init(&(thpool_p->thcount_lock), NULL);
    pthread_cond_init(&thpool_p->threads_all_idle, NULL);

    /* Thread init */
    int n;
    for (n = 0; n < num_threads; n++) {
        thread_init(thpool_p, &thpool_p->threads[n], n);
    }

    /* Wait for threads to initialize */
    while (thpool_p->num_threads_alive != num_threads) {}

    return thpool_p;
}

/**
 * @brief Add work to the job queue
 *
 * Takes an action and its argument and adds it to the threadpool's job queue.
 * If you want to add to work a function with more than one arguments then
 * a way to implement this is by passing a pointer to a structure.
 *
 * NOTICE: You have to cast both the function and argument to not get warnings.
 *
 * @example
 *
 *    void print_num(int num){
 *       printf("%d\n", num);
 *    }
 *
 *    int main() {
 *       ..
 *       int a = 10;
 *       thpool_add_job(thpool, (void*)print_num, (void*)a);
 *       ..
 *    }
 *
 * @param  threadpool    threadpool to which the work will be added
 * @param  function      pointer to function to add as work
 * @param  arg           pointer to an argument
 * @return 0 on successs, -1 otherwise.
 */
static int thpool_add_job(thpool_* thpool_p,
                          thpool_function_type function,
                          void* arg,
                          float* elapsed) {
    job* newjob;

    newjob = (struct job*) malloc(sizeof(struct job));
    if (newjob == NULL) {
        fprintf(stderr, "thpool_add_job(): Could not allocate memory for new job\n");
        return -1;
    }

    /* add function and argument */
    newjob->function = function;
    newjob->arg = arg;
    newjob->elapsed = elapsed;

    /* add job to queue */
    jobqueue_push(thpool_p, newjob);

    return 0;
}

static int thpool_add_jobs(thpool_* thpool_p,
                           thpool_function_type functions[],
                           void* args[],
                           float elapsed[],
                           int order[],
                           int nJobs) {
    job* newjobs[nJobs];
    int i;
    int j;

    for(i = 0; i < nJobs; ++i) {
        newjobs[i] = (struct job*) malloc(sizeof(struct job));
        if (newjobs[i] == NULL) {
            fprintf(stderr, "thpool_add_jobs(): Could not allocate memory for new jobs\n");
            return -1;
        }

        j = order != NULL ? order[i] : i;
        /* add function and argument */
        newjobs[i]->function = functions[j];
        newjobs[i]->arg = args[j];
        if(elapsed != NULL)
            newjobs[i]->elapsed = &elapsed[j];
        else
            newjobs[i]->elapsed = NULL;
    }

    /* add job to queue */
    jobqueue_multipush(thpool_p, newjobs, nJobs);

    return 0;
}

/**
 * @brief Wait for all queued jobs to finish
 *
 * Will wait for all jobs - both queued and currently running to finish.
 * Once the queue is empty and all work has completed, the calling thread
 * (probably the main program) will continue.
 *
 * Smart polling is used in wait. The polling is initially 0 - meaning that
 * there is virtually no polling at all. If after 1 seconds the threads
 * haven't finished, the polling interval starts growing exponentially
 * untill it reaches max_secs seconds. Then it jumps down to a maximum polling
 * interval assuming that heavy processing is being used in the threadpool.
 *
 * @example
 *
 *    ..
 *    threadpool thpool = thpool_init(4);
 *    ..
 *    // Add a bunch of work
 *    ..
 *    thpool_wait(thpool);
 *    puts("All added work has finished");
 *    ..
 *
 * @param threadpool     the threadpool to wait for
 * @return nothing
 */
static void thpool_wait(thpool_* thpool_p) {
    pthread_mutex_lock(&thpool_p->thcount_lock);
    while (thpool_p->jobqueue_p->len || thpool_p->num_threads_working) {
        pthread_cond_wait(&thpool_p->threads_all_idle, &thpool_p->thcount_lock);
    }
    thpool_p->jobqueue_p->total_time = 0;
    thpool_p->jobqueue_p->highest_expected_return = 0;
    pthread_mutex_unlock(&thpool_p->thcount_lock);
}


/**
 * @brief Destroy the threadpool
 *
 * This will wait for the currently active threads to finish and then 'kill'
 * the whole threadpool to free up memory.
 *
 * @example
 * int main() {
 *    threadpool thpool1 = thpool_init(2);
 *    threadpool thpool2 = thpool_init(2);
 *    ..
 *    thpool_destroy(thpool1);
 *    ..
 *    return 0;
 * }
 *
 * @param threadpool     the threadpool to destroy
 * @return nothing
 */
static void thpool_destroy(thpool_* thpool_p) {
    /* No need to destory if it's NULL */
    if (thpool_p == NULL) return;

    volatile int threads_total = thpool_p->num_threads_alive;

    /* End each thread 's infinite loop */
    thpool_p->threads_keepalive = 0;

    /* Give one second to kill idle threads */
    double TIMEOUT = 1.0;
    time_t start, end;
    double tpassed = 0.0;
    time(&start);
    while (tpassed < TIMEOUT && thpool_p->num_threads_alive) {
        bsem_post_all(thpool_p->jobqueue_p->has_jobs);
        time(&end);
        tpassed = difftime(end, start);
    }

    /* Poll remaining threads */
    while (thpool_p->num_threads_alive) {
        bsem_post_all(thpool_p->jobqueue_p->has_jobs);
        sleep(1);
    }

    /* Job queue cleanup */
    jobqueue_destroy(thpool_p);
    free(thpool_p->jobqueue_p);

    /* Deallocs */
    int n;
    for (n = 0; n < threads_total; n++) {
        thread_destroy(thpool_p->threads[n]);
    }
    free(thpool_p->threads);
    free(thpool_p);
    
    if(cppadcg_pool_verbose) {
        fprintf(stdout, "thpool_destroy(): thread pool destroyed\n");
    }
}

/* ============================ THREAD ============================== */


/* Initialize a thread in the thread pool
 *
 * @param thread        address to the pointer of the thread to be created
 * @param id            id to be given to the thread
 * @return 0 on success, -1 otherwise.
 */
static int thread_init(thpool_* thpool_p,
                       struct thread** thread_p,
                       int id) {

    *thread_p = (struct thread*) malloc(sizeof(struct thread));
    if (thread_p == NULL) {
        fprintf(stderr, "thread_init(): Could not allocate memory for thread\n");
        return -1;
    }

    (*thread_p)->thpool_p = thpool_p;
    (*thread_p)->id = id;

    pthread_create(&(*thread_p)->pthread, NULL, (void*) thread_do, (*thread_p));
    pthread_detach((*thread_p)->pthread);
    return 0;
}

static float get_thread_time(struct timespec* cputime,
                             int* info) {
    *info = clock_gettime(CLOCK_THREAD_CPUTIME_ID, cputime);
    if(*info == 0) {
        return cputime->tv_sec + cputime->tv_nsec * 1e-9f;
    } else {
        fprintf(stderr, "failed clock_gettime()\n");
        return 0;
    }
}

static float get_monotonic_time(struct timespec* time,
                                int* info) {
    *info = clock_gettime(CLOCK_MONOTONIC, time);
    if(*info == 0) {
        return time->tv_sec + time->tv_nsec * 1e-9f;
    } else {
        fprintf(stderr, "failed clock_gettime()\n");
        return 0;
    }
}

/* What each thread is doing
*
* In principle this is an endless loop. The only time this loop gets interrupted is once
* thpool_destroy() is invoked or the program exits.
*
* @param  thread        thread that will run this function
* @return nothing
*/
static void* thread_do(struct thread* thread_p) {
    float elapsed;
    int info;
    struct timespec cputime;
    int i;
    /* Set thread name for profiling and debugging */
    char thread_name[128] = {0};
    sprintf(thread_name, "thread-pool-%d", thread_p->id);

#if defined(__linux__)
    /* Use prctl instead to prevent using _GNU_SOURCE flag and implicit declaration */
    prctl(PR_SET_NAME, thread_name);
#elif defined(__APPLE__) && defined(__MACH__)
    pthread_setname_np(thread_name);
#else
    fprintf(stderr, "thread_do(): pthread_setname_np is not supported on this system");
#endif

    /* Assure all threads have been created before starting serving */
    thpool_* thpool_p = thread_p->thpool_p;

    /* Mark thread as alive (initialized) */
    pthread_mutex_lock(&thpool_p->thcount_lock);
    thpool_p->num_threads_alive += 1;
    pthread_mutex_unlock(&thpool_p->thcount_lock);

    while (thpool_p->threads_keepalive) {

        bsem_wait(thpool_p->jobqueue_p->has_jobs);

        if (thpool_p->threads_keepalive) {

            pthread_mutex_lock(&thpool_p->thcount_lock);
            thpool_p->num_threads_working++;
            pthread_mutex_unlock(&thpool_p->thcount_lock);

            /* Read job from queue and execute it */
            thpool_function_type func_buff;
            void* arg_buff;
            job* job_p;
            job_group* work_group;
            pthread_mutex_lock(&thpool_p->jobqueue_p->rwmutex);
            work_group = jobqueue_pull(thpool_p);
            pthread_mutex_unlock(&thpool_p->jobqueue_p->rwmutex);

            if(cppadcg_pool_verbose) {
                fprintf(stdout, "Thread executing %i jobs\n", work_group->size);
            }
            
            for (i = 0; i < work_group->size; ++i) {
                job_p = &work_group->jobs[i];
                int do_benchmark = job_p->elapsed != NULL && (*job_p->elapsed) == 0;
                if (do_benchmark) {
                    elapsed = -get_thread_time(&cputime, &info);
                }

                /* Execute the job */
                func_buff = job_p->function;
                arg_buff = job_p->arg;
                func_buff(arg_buff);

                if (do_benchmark && info == 0) {
                    elapsed += get_thread_time(&cputime, &info);
                    if (info == 0) {
                        (*job_p->elapsed) = elapsed;
                    }
                }
            }
            free(work_group->jobs);
            free(work_group);

            pthread_mutex_lock(&thpool_p->thcount_lock);
            thpool_p->num_threads_working--;
            if (!thpool_p->num_threads_working) {
                pthread_cond_signal(&thpool_p->threads_all_idle);
            }
            pthread_mutex_unlock(&thpool_p->thcount_lock);

        }
    }
    pthread_mutex_lock(&thpool_p->thcount_lock);
    thpool_p->num_threads_alive--;
    pthread_mutex_unlock(&thpool_p->thcount_lock);

    return NULL;
}


/* Frees a thread  */
static void thread_destroy(thread* thread_p) {
    free(thread_p);
}


/* ============================ JOB QUEUE =========================== */


/* Initialize queue */
static int jobqueue_init(thpool_* thpool_p) {

    thpool_p->jobqueue_p = (struct jobqueue*) malloc(sizeof(struct jobqueue));
    if (thpool_p->jobqueue_p == NULL) {
        return -1;
    }
    thpool_p->jobqueue_p->len = 0;
    thpool_p->jobqueue_p->front = NULL;
    thpool_p->jobqueue_p->rear = NULL;
    thpool_p->jobqueue_p->total_time = 0;
    thpool_p->jobqueue_p->highest_expected_return = 0;

    thpool_p->jobqueue_p->has_jobs = (struct bsem*) malloc(sizeof(struct bsem));
    if (thpool_p->jobqueue_p->has_jobs == NULL) {
        return -1;
    }

    pthread_mutex_init(&(thpool_p->jobqueue_p->rwmutex), NULL);
    bsem_init(thpool_p->jobqueue_p->has_jobs, 0);

    return 0;
}


/* Clear the queue */
static void jobqueue_clear(thpool_* thpool_p) {
    int i;
    struct job_group* group;
    int size;

    do {
        group = jobqueue_pull(thpool_p);
        size = group->size;
        free(group->jobs);
        free(group);
    } while (size > 0);

    thpool_p->jobqueue_p->front = NULL;
    thpool_p->jobqueue_p->rear = NULL;
    bsem_reset(thpool_p->jobqueue_p->has_jobs);
    thpool_p->jobqueue_p->len = 0;
    thpool_p->jobqueue_p->total_time = 0;
    thpool_p->jobqueue_p->highest_expected_return = 0;
}


/**
 * Add (allocated) job to queue without locks (internal function)
 */
static void jobqueue_push_internal(thpool_* thpool_p,
                                   struct job* newjob) {
    newjob->prev = NULL;

    switch (thpool_p->jobqueue_p->len) {

        case 0:  /* if no jobs in queue */
            thpool_p->jobqueue_p->front = newjob;
            thpool_p->jobqueue_p->rear = newjob;
            break;

        default: /* if jobs in queue */
            thpool_p->jobqueue_p->rear->prev = newjob;
            thpool_p->jobqueue_p->rear = newjob;

    }
    if(newjob->elapsed != NULL) {
        thpool_p->jobqueue_p->total_time += *newjob->elapsed;
    }
    thpool_p->jobqueue_p->len++;
}

/**
 * Add (allocated) job to queue
 */
static void jobqueue_push(thpool_* thpool_p,
                          struct job* newjob) {
    pthread_mutex_lock(&thpool_p->jobqueue_p->rwmutex);

    jobqueue_push_internal(thpool_p, newjob);

    bsem_post(thpool_p->jobqueue_p->has_jobs);

    pthread_mutex_unlock(&thpool_p->jobqueue_p->rwmutex);
}


/**
 * Add (allocated) multiple jobs to queue
 */
static void jobqueue_multipush(thpool_* thpool_p,
                               struct job* newjob[],
                               int nJobs) {
    int i;

    pthread_mutex_lock(&thpool_p->jobqueue_p->rwmutex);

    for(i = 0; i < nJobs; ++i) {
        jobqueue_push_internal(thpool_p, newjob[i]);
    }

    bsem_post_all(thpool_p->jobqueue_p->has_jobs);

    pthread_mutex_unlock(&thpool_p->jobqueue_p->rwmutex);
}

static struct job* jobqueue_extract_single(thpool_* thpool_p) {
    struct job* job_p = thpool_p->jobqueue_p->front;

    switch (thpool_p->jobqueue_p->len) {
        case 0:  /* if no jobs in queue */
            return NULL;

        case 1:  /* if one job in queue */
            thpool_p->jobqueue_p->front = NULL;
            thpool_p->jobqueue_p->rear = NULL;
            thpool_p->jobqueue_p->len = 0;
            thpool_p->jobqueue_p->total_time = 0;
            thpool_p->jobqueue_p->highest_expected_return = 0;
            return job_p;

        default: /* if >1 jobs in queue */
            thpool_p->jobqueue_p->front = job_p->prev;
            thpool_p->jobqueue_p->len--;
            if(job_p->elapsed != NULL) {
                thpool_p->jobqueue_p->total_time -= *job_p->elapsed;
            }
            return job_p;
    }
}

static void jobqueue_extract_single_group(thpool_* thpool_p,
                                          job_group* group) {
    struct job* job_p = jobqueue_extract_single(thpool_p);
    if(job_p != NULL) {
        group->size = 1;
        group->jobs = (struct job*) malloc(sizeof(struct job));
        group->jobs[0] = *job_p; // copy
        free(job_p);
    } else {
        group->size = 0;
        group->jobs = NULL;
    }
}

/**
 * Get jobs from the queue(removes them from the queue)
 *
 * Notice: Caller MUST hold a mutex
 */
static struct job_group* jobqueue_pull(thpool_* thpool_p) {

    job_group* group;
    group = (struct job_group*) malloc(sizeof(struct job_group));
    struct job* job_p;
    float current_time;
    float duration, duration_next, min_duration, target_duration;
    struct timespec timeAux;
    int info;
    int i;

    if (group_gen_strategy == SINGLE_JOB || thpool_p->jobqueue_p->len <= 1 || thpool_p->jobqueue_p->total_time <= 0) {
        if (cppadcg_pool_verbose) {
            if(group_gen_strategy == MULTI_JOB) {
                if (thpool_p->jobqueue_p->len <= 1)
                    fprintf(stdout, "jobqueue_pull(): Work group with 1 job\n");
                else if (thpool_p->jobqueue_p->total_time <= 0)
                    fprintf(stdout, "jobqueue_pull(): Using single-job instead of multi-job (no timing information)\n");
            }
        }

        // SINGLE_JOB
        jobqueue_extract_single_group(thpool_p, group);
    } else {
        // MULTI_JOB
        job_p = thpool_p->jobqueue_p->front;

        if (job_p->elapsed == NULL) {
            if (cppadcg_pool_verbose) {
                fprintf(stderr, "jobqueue_pull(): Using single job instead of multi-job (No timing information for current job)\n");
            }
            // cannot use this strategy (something went wrong!)
            jobqueue_extract_single_group(thpool_p, group);

        } else {
            // there are at least 2 jobs in the queue
            group->size = 1;
            duration = *job_p->elapsed;
            duration_next = duration;
            job_p = job_p->prev;
            target_duration = thpool_p->jobqueue_p->total_time * cppadcg_pool_multijob_maxgroupwork / thpool_p->num_threads; // always positive
            current_time = get_monotonic_time(&timeAux, &info);

            if (thpool_p->jobqueue_p->highest_expected_return > 0 && info) {
                min_duration = 0.9f * (thpool_p->jobqueue_p->highest_expected_return - current_time);
                if (target_duration < min_duration) {
                    target_duration = min_duration;
                }
            }

            do {
                if (job_p->elapsed == NULL) {
                    break;
                }
                duration_next += *job_p->elapsed;
                if (duration_next < target_duration) {
                    group->size++;
                    duration = duration_next;
                } else {
                    break;
                }
                job_p = job_p->prev;
            } while (job_p != thpool_p->jobqueue_p->front);

            if (cppadcg_pool_verbose) {
                fprintf(stdout, "jobqueue_pull(): Work group with %i jobs for %e s (target: %e s)\n", group->size, duration, target_duration);
            }

            group->jobs = (struct job*) malloc(group->size * sizeof(struct job));
            for (i = 0; i < group->size; ++i) {
                job_p = jobqueue_extract_single(thpool_p);
                group->jobs[i] = *job_p; // copy
                free(job_p);
            }

            duration_next = current_time + duration; // the time when the current work is expected to end
            if(duration_next > thpool_p->jobqueue_p->highest_expected_return)
                thpool_p->jobqueue_p->highest_expected_return = duration_next;
        }

    }

    /* more than one job in queue -> post it */
    if (thpool_p->jobqueue_p->len > 0) {
        bsem_post(thpool_p->jobqueue_p->has_jobs);
    }

    return group;
}


/* Free all queue resources back to the system */
static void jobqueue_destroy(thpool_* thpool_p) {
    jobqueue_clear(thpool_p);
    free(thpool_p->jobqueue_p->has_jobs);
}





/* ======================== SYNCHRONISATION ========================= */


/* Init semaphore to 1 or 0 */
static void bsem_init(bsem* bsem_p, int value) {
    if (value < 0 || value > 1) {
        fprintf(stderr, "bsem_init(): Binary semaphore can take only values 1 or 0");
        exit(1);
    }
    pthread_mutex_init(&(bsem_p->mutex), NULL);
    pthread_cond_init(&(bsem_p->cond), NULL);
    bsem_p->v = value;
}


/* Reset semaphore to 0 */
static void bsem_reset(bsem* bsem_p) {
    bsem_init(bsem_p, 0);
}


/* Post to at least one thread */
static void bsem_post(bsem* bsem_p) {
    pthread_mutex_lock(&bsem_p->mutex);
    bsem_p->v = 1;
    pthread_cond_signal(&bsem_p->cond);
    pthread_mutex_unlock(&bsem_p->mutex);
}


/* Post to all threads */
static void bsem_post_all(bsem* bsem_p) {
    pthread_mutex_lock(&bsem_p->mutex);
    bsem_p->v = 1;
    pthread_cond_broadcast(&bsem_p->cond);
    pthread_mutex_unlock(&bsem_p->mutex);
}


/* Wait on semaphore until semaphore has value 0 */
static void bsem_wait(bsem* bsem_p) {
    pthread_mutex_lock(&bsem_p->mutex);
    while (bsem_p->v != 1) {
        pthread_cond_wait(&bsem_p->cond, &bsem_p->mutex);
    }
    bsem_p->v = 0;
    pthread_mutex_unlock(&bsem_p->mutex);
}
