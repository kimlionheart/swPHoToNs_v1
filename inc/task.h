#ifndef TASK_H
#define TASK_H

#define NCORE 64
#define NSEG 10
#define NPACK 16
#define NVECT 4

void fmm_task_parallel_kernel(int gang);
void fmm_task_parallel(int ) ;
void fmm_task_parallel_kernel_athread();


#endif
