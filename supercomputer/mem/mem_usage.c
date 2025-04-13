#include <sys/resource.h>
#include <sys/time.h>
#include <mach/mach.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

long get_footprint(){
    // Get memory regions info
    pid_t pid = getpid();
    task_t task = MACH_PORT_NULL;
    
    if (task_for_pid(mach_task_self(), pid, &task) == KERN_SUCCESS) {
        task_vm_info_data_t vm_info;
        mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
        if (task_info(task, TASK_VM_INFO, (task_info_t)&vm_info, &count) == KERN_SUCCESS) {
            // Returns physical footprint in bytes (includes private and shared memory)
            return vm_info.phys_footprint;
        }
        else {
            return -1;
        }
    }
    else {
        return -1;
    }
}

// long get_maxrss() {
//     struct rusage usage;
//     if (getrusage(RUSAGE_SELF, &usage) == 0) {
//         // All three metrics are in bytes, but measure different aspects of memory usage:
//         // ru_maxrss: Maximum resident set size (historical peak)
//         // foot: Current physical footprint (includes shared memory)
//         long foot = get_footprint();

//         // Print the three memory metrics
//         printf("Memory metrics (bytes): ru_maxrss=%ld, phys_footprint=%ld\n", 
//                usage.ru_maxrss, foot);
               
//         return max_three(usage.ru_maxrss, mach, foot);  // Returns the largest value in bytes
//     } else {
//         return -1;  // Error indicator
//     }
// }

