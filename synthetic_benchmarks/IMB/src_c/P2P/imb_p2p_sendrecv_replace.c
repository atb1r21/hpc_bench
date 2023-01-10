/*
*****************************************************************************
*                                                                           *
* Copyright Intel Corporation.                                              *
*                                                                           *
*****************************************************************************

This code is covered by the Community Source License (CPL), version
1.0 as published by IBM and reproduced in the file "license.txt" in the
"license" subdirectory. Redistribution in source and binary form, with
or without modification, is permitted ONLY within the regulations
contained in above mentioned license.

Use of the name and trademark "Intel(R) MPI Benchmarks" is allowed ONLY
within the regulations of the "License for Use of "Intel(R) MPI
Benchmarks" Name and Trademark" as reproduced in the file
"use-of-trademark-license.txt" in the "license" subdirectory.

THE PROGRAM IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED INCLUDING, WITHOUT
LIMITATION, ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT,
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Each Recipient is
solely responsible for determining the appropriateness of using and
distributing the Program and assumes all risks associated with its
exercise of rights under this Agreement, including but not limited to
the risks and costs of program errors, compliance with applicable
laws, damage to or loss of data, programs or equipment, and
unavailability or interruption of operations.

EXCEPT AS EXPRESSLY SET FORTH IN THIS AGREEMENT, NEITHER RECIPIENT NOR
ANY CONTRIBUTORS SHALL HAVE ANY LIABILITY FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING
WITHOUT LIMITATION LOST PROFITS), HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OR
DISTRIBUTION OF THE PROGRAM OR THE EXERCISE OF ANY RIGHTS GRANTED
HEREUNDER, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

EXPORT LAWS: THIS LICENSE ADDS NO RESTRICTIONS TO THE EXPORT LAWS OF
YOUR JURISDICTION. It is licensee's responsibility to comply with any
export regulations applicable in licensee's jurisdiction. Under
CURRENT U.S. export regulations this software is eligible for export
from the U.S. and can be downloaded by or otherwise exported or
reexported worldwide EXCEPT to U.S. embargoed destinations which
include Cuba, Iraq, Libya, North Korea, Iran, Syria, Sudan,
Afghanistan and any other country to which the U.S. has embargoed
goods and services.
*/

#include "imb_p2p.h"

void imb_p2p_sendrecv_replace() {
    char *buffer;
    size_t msg_size_index;
    int nranks = imb_p2p_config.nranks;
    if (nranks < 2) {
        if (imb_p2p_config.rank == 0) {
            fprintf(unit, "\n");
            fprintf(unit, "# !! Benchmark %s is invalid for %d process !!\n", IMB_P2P_SENDRECV_REPLACE, nranks);
            fflush(unit);
        }
        return;
    }
    if (imb_p2p_config.rank == 0) {
        imb_p2p_print_benchmark_header(IMB_P2P_SENDRECV_REPLACE);
        fprintf(unit, " %12s %12s %12s %12s %12s\n", "#bytes", "#repetitions", "t[usec]", "Mbytes/sec", "Msg/sec");
        fflush(unit);
    }
    buffer = (char *)imb_p2p_alloc_mem(imb_p2p_config.messages.max_size);
    memset(buffer, imb_p2p_config.rank, imb_p2p_config.messages.max_size);
    for (msg_size_index = 0; msg_size_index < imb_p2p_config.messages.length; msg_size_index++) {
        size_t size = imb_p2p_config.messages.array[msg_size_index];
        size_t iteration, number_of_iterations, number_of_warm_up_iterations;
        double time;
        int dst = (imb_p2p_config.rank + 1) % nranks;
        int src = (imb_p2p_config.rank + nranks - 1) % nranks;
        get_iters(size, &number_of_iterations, &number_of_warm_up_iterations);
        imb_p2p_pause();
        imb_p2p_barrier(MPI_COMM_WORLD);
        for (iteration = 0; iteration < number_of_warm_up_iterations; iteration++) {
            touch_send_buff(size, buffer);
            MPI_Sendrecv_replace(buffer, size, MPI_BYTE,
                                 dst, 0, src, 0,
                                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        touch_recv_buff(size, buffer);
        time = MPI_Wtime();
        for (iteration = 0; iteration < number_of_iterations; iteration++) {
            touch_send_buff(size, buffer);
            MPI_Sendrecv_replace(buffer, size, MPI_BYTE,
                                 dst, 0, src, 0,
                                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        touch_recv_buff(size, buffer);
        time = MPI_Wtime() - time;
        imb_p2p_pause();
        imb_p2p_barrier(MPI_COMM_WORLD);
        if (imb_p2p_config.rank) {
            MPI_Send(&time, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        } else {
            const double microseconds = 1000000.0;
            double bandwidth, message_rate, latency;
            double max_time = time;
            double aggregate_time = time;
            int source_rank;
            for (source_rank = 1; source_rank < nranks; source_rank++) {
                double t = 0;
                MPI_Recv(&t, 1, MPI_DOUBLE, source_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                aggregate_time += t;
                if (max_time < t) {
                    max_time = t;
                }
            }
            bandwidth = ((number_of_iterations * nranks * size) / (microseconds * max_time));
            message_rate = ((number_of_iterations * nranks) / max_time);
            latency = (microseconds * aggregate_time) / (number_of_iterations * nranks);
            fprintf(unit, " %12" PRIu64 " %12" PRIu64 " %12.2f %12.2f %12.0f\n", size, number_of_iterations, latency, bandwidth, message_rate);
            fflush(unit);
        }
    }
    imb_p2p_free_mem(buffer);
}
