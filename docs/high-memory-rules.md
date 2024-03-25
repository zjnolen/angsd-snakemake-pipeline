# Rules using large amounts of RAM

NOTE: This is a work in progress list. Trying to figure out what

The biggest challenge with using this pipeline with other datasets is ensuring
RAM is properly allocated. Many rules require very little RAM, and so the
default allocations that come on your cluster per thread will likely do fine.
However, some rules require considerably more RAM. These are:

