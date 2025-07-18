(ref_infrastructure)=
# Infrastructure Overview
This page describe the main infrastucutres of the icon model such as the parallelization and I/O.

(ref_infrastructure_parallelization)=
## Parallelization
The ICON model leverages parallelization approaches to enhance computational efficiency and performance. MPI (Message Passing Interface) is employed for distributed memory parallelization, specifically distributing the computational workload horizontally across multiple nodes. This means that the horizontal domain of the model is divided into smaller subdomains, distributed across nodes, allowing for efficient communication and computation on a distributed system. OpenMP (Open Multi-Processing) is used for shared memory parallelization, allowing multiple threads to execute concurrently within a single node, thus speeding up computations. For GPU acceleration, OpenACC compiler directives are utilized to offload compute-intensive tasks to GPUs, significantly reducing runtime by exploiting the massive parallelism offered by modern GPUs. These parallelization strategies collectively ensure that the ICON model can handle the complex and large-scale computations required for accurate and timely weather forecasting.

A more detailed description of the parallelization is provided in the **{term}`ICON Tutorial 2024`** section 8.

(ref_infrastructure_io)=
## Input & Output (I/O)

Detail information to Input and Output can be found in the **{term}`ICON Tutorial 2024`** in respectively section 2 and 7.
