netlogo-pea
===========

Pathway Evolution Algorithm (PEA), a minimalist genetic algorithm working on paths in a network and which reproduces some neurobiological phenomenons

## WHAT IS IT?

Pathway Evolution Algorithm (PEA) by Chrisantha Fernando et al. (2011), implemented in NetLogo 5.0.5 by Stephen Larroque (2014).

Also implemented Microbial Genetic Algorithm by I. Harvey (1996) as a comparison (and because PEA is a special case of MGA).

This is a general, minimalist genetic algorithm but with the specificity that it works on pathways in a network as units of evolution instead of individuals.

An MGA algorithm was implemented to compare the phenomenons that appears with PEA.

PEA is NOT meant to be more performant than MGA (even if that's often the case) or any other algorithm, but mostly to demonstrate that some neurobiological phenomenons can naturally emerge from paths in a network, and PEA demonstrate that with a very simple genetic algorithm (since PEA is really just a MGA applied to paths instead of individuals).

## HOW IT WORKS

See the sourcecode for a full algorithm summary in the comments.

## HOW TO USE IT

Select a phenotypes-kind, then a reward-kind, the n-layers (number of layers = number of nodes = length of phenotypes), n-rows (number of distinctive phenotypes on each layer, eg: with phenotypes-kind = Binary, with n-rows = 1 it will initialize with only one row of 0 and 1, if n-rows = 2 there will be two rows with 0 and 1 on each layer).

Then press Setup and Go to see the network expand, contract and learn how to best optimize the reward-kind you chose given the phenotypes-kind.

## THINGS TO NOTICE

There are a few neurobiological phenomenons that can be reproduced with this model, like the expansion/contraction (watch the number of nodes and edges in the chart at the right-side), memory (ability to reuse previously expanded edges and nodes, use reward-kind = Oscillating count and set neuron-drop-after-time = 0) and disequilibrium (ability to maintain several equivalent solutions, use reward-kind = Disequilibrium separate or shift and set neuron-drop-after-time = 0).

## THINGS TO TRY

You can try to play with the parameters (most specifically with synapse-initial-weight-on-mutation, synapse-drop-below-weight and neuron-drop-after-time).

You can also try to enable the (unofficial) extensions of the model:
- redundancy? will break the phenotype unicity constraint (in the base model, a neuron must have a unique phenotype in its layer: no other neuron in this layer must have the same phenotype), and thus allow for more redundancy since several paths can have the same phenotype. This is may be a fundamental property of neurobiological phenomenons (see Claude Berrou's Turbocodes).
- variable-layers will allow the network to increase or reduce the number of layers as needed. This allows to find solutions that are beyond the initial n-layers number. Try this with the reward-kind = "Guess-a-number" or "Min length".

## EXTENDING THE MODEL

- Add more plots.
- Enhance the variable-layers extension: it works well for "Min length" (reducing the network size) but not so well with "Guess-a-number" (increasing the size to reach the required number to guess).

## NETLOGO FEATURES

See the AUX FUNCTIONS part for generic NetLogo functions like ensemble operations, arg-max/arg-min and roulette-wheel.

## CREDITS AND REFERENCES

- "Evolvable Neuronal Paths: A Novel Basis for Information and Search in the Brain", 2011, Fernando C, Vasas V, Szathmáry E, Husbands P. PLoS ONE 6(8): e23534. doi: 10.1371/journal.pone.0023534

Note that another implementation in Objective-C of this algorithm by the authors themselves is provided with the article.

- "Microbial Genetic Algorithm", 1996, I. Harvey

## AUTHOR AND LICENSE

Original Pathway Evolution Algorithm by Chrisantha Fernando et al, implemented in NetLogo by Stephen Larroque, opensource licensed under MIT.
