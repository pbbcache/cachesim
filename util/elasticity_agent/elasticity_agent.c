#include "schedctl.h"
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
/* For schedctl control */
#include <fcntl.h>
#include <unistd.h> 
#include <sys/mman.h>
#include <stdio.h>
#include <sys/wait.h>
#include <string.h>
#ifndef PAGE_SIZE
#define PAGE_SIZE 4096
#endif
#define SCHEDCTL_PROC_ENTRY "/proc/pmc/schedctl"
//#define DEBUG

schedctl_t* schedctl_data = NULL;
unsigned int current_worker_count = 0;
unsigned int max_worker_count = 0;
int schedctl_fd = -1;

static int test_schedctl_support()
{
    int nr_tries=0;
    static const int max_tries=5;

    /* Just make sure the proc entry exists. If it does, close it.*/
retry:
    int fd = open(SCHEDCTL_PROC_ENTRY, O_RDWR);
    if (fd < 0) {
        perror("schedctl_open()");

        nr_tries++;

        if (nr_tries>=max_tries)
            return 0;
        else 
            goto retry;    
   
        return 0;
    }
    else {
        close(fd);
        return 1;
    }
}


static inline void schedctl_write (int fd, const char* str)
{
    if (fd>=0)
        write(fd,str,strlen(str));
}

static void set_master_thread()
{
    schedctl_write(schedctl_fd, "set master\n");
}

static void unset_master_thread()
{
    schedctl_write(schedctl_fd, "unset master\n");
}


/* Default implementation for LINUX */
static schedctl_t *schedctl_retrieve(int *fd)
{
    static int first_time = 1;
    static int schedctl_support = 0;
    int configfd = -1;
    schedctl_t *schedctl = NULL;

    /* Reset fields ... */
    (*fd) = -1;
    schedctl = NULL;

    if (first_time) {
        schedctl_support = test_schedctl_support();
        first_time = 0;
    }

    if (schedctl_support) {
        configfd = open(SCHEDCTL_PROC_ENTRY, O_RDWR);
        if (configfd < 0) {
            perror("open");
            return NULL;
        }

        schedctl = (schedctl_t *)mmap(NULL, PAGE_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, configfd, 0);

        if (schedctl == MAP_FAILED) {
            perror("mmap");
            return NULL;
        }
#ifdef DEBUG
        printf("Mmap Ok. Address:0x%p\n", schedctl);
        printf("%d %d\n", schedctl->sc_nfc, schedctl->sc_sf);
#endif
    } else {
        /* *
         * Allocate a FAKE schedctl structure.
         * This is just to make sure the runtime works
         * even without kernel extensions.
         */
        schedctl = malloc(sizeof(schedctl_t));
        schedctl->sc_coretype = 1;
        schedctl->sc_prio = 0;
        schedctl->sc_nfc = 0;
        schedctl->sc_sf = 100;
        schedctl->sc_spinning = 0;
        schedctl->sc_num_threads = current_worker_count;
        schedctl->sc_malleable = 0;
        configfd = -1;
//#ifdef DEBUG
        fprintf(stderr, "FAKE SCHEDCTL\n");
//#endif
    }
    (*fd) = configfd;
    return schedctl;
}

void schedctl_release(int fd, schedctl_t* data) {
	if (data)
		munmap(data, PAGE_SIZE);

	if (fd!=-1)
		close(fd);
	else
		free(data); //emulated version	
}

int main(int argc, char *argv[]){

    sigset_t set;
    int sig;
    unsigned char finished = 0;
    pid_t child=0, pid;
#ifdef DEBUG
    unsigned char ascending = 1;
#endif
    if (argc<3) {
        fprintf(stderr, "Usage: %s <initial_worker_count> <max_worker_count> [optional command] \n", argv[0]);
        exit(1);
    }

    if (sscanf(argv[1],"%d",&current_worker_count)!=1) {
        fprintf(stderr, "Argument 1 is not an integer: %s\n", argv[1]);
        exit(2);    
    }

    if (sscanf(argv[2],"%d",&max_worker_count)!=1) {
        fprintf(stderr, "Argument 2 is not an integer: %s\n", argv[2]);
        exit(2);    
    }

    if ((schedctl_data=schedctl_retrieve(&schedctl_fd))==NULL){
        fprintf(stderr, "Cannot obtain reference to schedctl structure\n");
        exit(2);    
    }

    // Create a set of signals to wait for
    sigemptyset(&set);
    sigaddset(&set, SIGUSR1);
    sigaddset(&set, SIGUSR2);
    sigaddset(&set, SIGTERM);
    sigaddset(&set, SIGCHLD);
    sigaddset(&set, SIGINT);

    // Block signals in the set
    sigprocmask(SIG_BLOCK, &set, NULL);


    
    /**
     * Enable malleability as we are ready
     * for action 
     **/
    set_master_thread();
    schedctl_data->sc_malleable=1;

    /* Launch command as workload, to generate certain load in container */
    if (argc>3) {
        child=fork();

        switch (child) {
            case 0: 
                execvp(argv[3],&argv[3]);
                fprintf(stderr, "Exec failed\n");
                exit(1);
                break;
            case -1:
                fprintf(stderr, "Fork failed\n");
                exit(1);
                break;            
            default:
                /* parent does nothing for now */
                break; 
        }
    }

    // Wait for a signal
    // This is an indicator for the other end, so that it is safe to send signals
    printf("Elasticity agent: (PID=%d) waiting for signals...\n", getpid());
    fflush(stdout);
    while (!finished)
    {
        if (sigwait(&set, &sig) != 0)
        {
            perror("sigwait");
            exit(3);
        }

        // Handle the received signal
        switch (sig)
        {
        case SIGUSR1:

#ifdef DEBUG
            switch(ascending){
                case 0:
                        if (schedctl_data->sc_num_threads==1){
                            ascending=1;
                            schedctl_data->sc_num_threads++;
                        } else {
                            schedctl_data->sc_num_threads--;  
                        }
                        break; 
                default:
                        if (schedctl_data->sc_num_threads==max_worker_count){
                            ascending=0;
                            schedctl_data->sc_num_threads--;
                        } else {
                            schedctl_data->sc_num_threads++;  
                        }
                        break;  
            }
            
            current_worker_count=schedctl_data->sc_num_threads;
            printf("%d\n",current_worker_count);
            fflush(stdout);
#else
            if (schedctl_data->sc_num_threads!=current_worker_count){
                current_worker_count=schedctl_data->sc_num_threads;
                /**
                 * Important note: This message should appear on the standard
                 * output, as the other end of the pipe is expecting to receive data
                 * via this stream
                */
                printf("%d\n",current_worker_count);
                fflush(stdout);
            }        
 #endif
            break;
        case SIGUSR2:
            /* Toggle malleability flag */
            schedctl_data->sc_malleable=schedctl_data->sc_malleable?0:1;
            fprintf(stderr, "Elasticity agent: Toggling malleability flag. Now %d\n",schedctl_data->sc_malleable);
            break;
        case SIGINT:
            if (child>0){
                kill(child, SIGKILL);
            }
        case SIGTERM:
        case SIGCHLD:
            if (argc>3) {
                if (child>0){
                    pid=wait(NULL);

                    if (pid==child)
                        fprintf(stderr, "Child terminated\n");
                }
            }
            fprintf(stderr, "Elasticity agent: Terminating\n");
            finished = 1;
            break;
        default:
            fprintf(stderr, "Elasticity agent: Received unknown signal\n");
        }
    }

    // Disable malleability 
    schedctl_data->sc_malleable=0;
    unset_master_thread();

    // Free-up schedctl resources
    schedctl_release(schedctl_fd, schedctl_data);

    return 0;
}




