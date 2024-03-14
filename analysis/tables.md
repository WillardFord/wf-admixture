| Variable | Description
|------|--------
| I | Number of Individuals |
| J | Number of Variants |
| K | Number of Populations |

| Matrix Var | Description | Dimension
|------|-------- | ---
| Q | Proportion of individual i's genome from population k | I x K
| F | MAF of variant j in population k | K x J
| G | Minor allele count of variant j in individual i |  I x J


|Tool | Runtime | # Iterations | Time/Iteration
|------|-------- | --- | ----
| wf-admixture | 8 hours, 24 minutes, 29 seconds | 3463 iterations | 8.74 seconds per iteration
| admixture | 21 seconds | 13 iterations | 1.62 seconds per iteration


| Tags | Description |
|------|-------- | 
| -k | number of populations |
| -q | number of threads for multithreading | 
| -t | threshold value for cutoff | 


