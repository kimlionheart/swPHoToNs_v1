#!/bin/bash
bsub -o 64.log -J cosmos -q q_sw_share -n 64 -cgsp 64 -share_size 6500 -host_stack 1000   -b ./photoNs-lcdm lcdm.run
