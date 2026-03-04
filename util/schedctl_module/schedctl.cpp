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
#include <errno.h>
#ifndef PAGE_SIZE
#define PAGE_SIZE 4096
#endif
#define SCHEDCTL_PROC_ENTRY "/proc/pmc/schedctl"
// #define DEBUG

#include <linux/types.h>
#ifndef PAGE_SIZE
#define PAGE_SIZE 4096
#endif

/* The datatype is private */
typedef struct
{
    volatile unsigned char sc_coretype;
    volatile int sc_prio; /* To denote priority among threads */
                          /* Legacy fields */
    volatile unsigned char sc_spinning;
    volatile unsigned int sc_nfc;
    volatile unsigned int sc_sf;
    volatile unsigned int sc_num_threads;
    volatile unsigned char sc_malleable;
} schedctl_t;

schedctl_t *schedctl_data = NULL;
unsigned int current_worker_count = 0;
int schedctl_fd = -1;

static int test_schedctl_support()
{
    int nr_tries = 0;
    static const int max_tries = 5;

    /* Just make sure the proc entry exists. If it does, close it.*/
retry:
    int fd = open(SCHEDCTL_PROC_ENTRY, O_RDWR);
    if (fd < 0)
    {
        perror("schedctl_open()");

        nr_tries++;

        if (nr_tries >= max_tries)
            return 0;
        else
            goto retry;

        return 0;
    }
    else
    {
        close(fd);
        return 1;
    }
}

static inline int schedctl_write(int fd, const char *str)
{
    if (fd == -1)
    {
        errno = ENOENT;
        return -1;
    }

    return write(fd, str, strlen(str));
}

int set_master_thread(void)
{
    return schedctl_write(schedctl_fd, "set master\n");
}

int unset_master_thread(void)
{
    return schedctl_write(schedctl_fd, "unset master\n");
}

/* Default implementation for LINUX */
static schedctl_t *_schedctl_retrieve(int *fd)
{
    static int first_time = 1;
    static int schedctl_support = 0;
    int configfd = -1;
    schedctl_t *schedctl = NULL;

    /* Reset fields ... */
    (*fd) = -1;
    schedctl = NULL;

    if (first_time)
    {
        schedctl_support = test_schedctl_support();
        first_time = 0;
    }

    if (schedctl_support)
    {
        configfd = open(SCHEDCTL_PROC_ENTRY, O_RDWR);
        if (configfd < 0)
        {
            perror("open");
            return NULL;
        }

        schedctl = (schedctl_t *)mmap(NULL, PAGE_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, configfd, 0);

        if (schedctl == MAP_FAILED)
        {
            perror("mmap");
            return NULL;
        }
#ifdef DEBUG
        printf("Mmap Ok. Address:0x%p\n", schedctl);
        printf("%d %d\n", schedctl->sc_nfc, schedctl->sc_sf);
#endif
    }
    else
    {
        /* *
         * Allocate a FAKE schedctl structure.
         * This is just to make sure the runtime works
         * even without kernel extensions.
         */
        schedctl = (schedctl_t*) malloc(sizeof(schedctl_t));
        schedctl->sc_coretype = 1;
        schedctl->sc_prio = 0;
        schedctl->sc_nfc = 0;
        schedctl->sc_sf = 100;
        schedctl->sc_spinning = 0;
        schedctl->sc_num_threads = current_worker_count;
        schedctl->sc_malleable = 0;
        configfd = -1;
        // #ifdef DEBUG
        fprintf(stderr, "FAKE SCHEDCTL\n");
        // #endif
    }
    (*fd) = configfd;
    return schedctl;
}

int schedctl_retrieve(void)
{
    if ((schedctl_data = _schedctl_retrieve(&schedctl_fd)) == NULL)
    {
        fprintf(stderr, "Cannot obtain reference to schedctl structure\n");
        return -1;
    }
    return 0;
}

static void _schedctl_release(int fd, schedctl_t *data)
{
    if (data)
        munmap(data, PAGE_SIZE);

    if (fd != -1)
        close(fd);
    else
        free(data); // emulated version
}

int schedctl_release(void)
{
    if (schedctl_data == NULL || schedctl_fd == -1)
    {
        errno = ENOENT;
        return -1;
    }
    _schedctl_release(schedctl_fd, schedctl_data);
    /* Clear global variables */
    schedctl_fd = -1;
    schedctl_data = NULL;
    return 0;
}

int enable_malleability(void)
{
    if (schedctl_data == NULL)
    {
        errno = ENOENT;
        return -1;
    }

    schedctl_data->sc_malleable = 1;
    return 0;
}


int disable_malleability(void)
{
    if (schedctl_data == NULL)
    {
        errno = ENOENT;
        return -1;
    }

    schedctl_data->sc_malleable = 0;
    return 0;
}

int check_malleability_status(void)
{
    if (schedctl_data == NULL)
    {
        errno = ENOENT;
        return -1;
    }

    return schedctl_data->sc_malleable;
}


int get_num_threads(void)
{
    if (schedctl_data == NULL)
    {
        errno = ENOENT;
        return -1;
    }

    return schedctl_data->sc_num_threads;
}

int set_num_threads(int nr_threads)
{
    if (schedctl_data == NULL)
    {
        errno = ENOENT;
        return -1;
    }

    if (nr_threads <=0)
    {
        errno = EINVAL;
        return -1;
    }

    schedctl_data->sc_num_threads=nr_threads;
    return 0;
}