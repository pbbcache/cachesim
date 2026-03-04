#ifndef PMC_SCHEDCTL_H
#define PMC_SCHEDCTL_H

int set_master_thread(void);
int unset_master_thread(void);
int schedctl_retrieve(void);
int schedctl_release(void);
int enable_malleability(void);
int disable_malleability(void);
int check_malleability_status(void);
int get_num_threads(void);
int set_num_threads(int nr_threads);

#endif