# RandomClone

This repo contains methods to generate random subclonal reconstructions from mutation data. There are four variants:

* Stick breaking : Sorts SNVs by their cancer cell fraction (CCF) and breaks the list into a random number of randomly sized chunks, where each chunk is a mutation cluster
* Uniform : Determines the CCF space occupied by SNVs (min CCF and max CCF) and draws a random number of random numbers from the range. Mutations are assigned through binomial assignment
* Informed : Runs the stick breaking procedure multiple times and uses [MutationTimer](https://github.com/gerstung-lab/MutationTime.R) to select the best model. Mutations are then assigned using MutationTimer
* Single : Assigns all SNVs to a single cluster
