; Pathway Evolution Algorithm (PEA) by Chrisantha Fernando et al. (2011), implemented in NetLogo 5.0.5 by Stephen Larroque (2014)
; Also implemented Microbial Genetic Algorithm by I. Harvey (1996) as a comparison (and because PEA is a special case of MGA)
; This sourcecode is licensed under FreeBSD License, feel free to do whatever you want with it.
; Ref to the paper: "Evolvable Neuronal Paths: A Novel Basis for Information and Search in the Brain", 2011, Fernando C, Vasas V, Szathmáry E, Husbands P. PLoS ONE 6(8): e23534. doi: 10.1371/journal.pone.0023534

globals[
  ; General
  paths-visu-height layers phenotypes
  layer-width layer-height network-height
  best-phenotypes ; list of best phenotypes to compare to compute the Hamming distance

  ; Oscillating Reward
  reward-oscillating-best-allele reward-oscillating-last-switch

  ; MGA
  mga-pop

  ; Plot
  pea-max-moving pea-max-moving-n pea-avg-moving pea-avg-moving-n pea-moving-last-time ; moving average for pea
  mga-max-moving mga-max-moving-n mga-avg-moving mga-avg-moving-n mga-moving-last-time ; moving average for mga
  ]
breed[neurons neuron]
directed-link-breed[synapses synapse]

synapses-own[weight sy-layer]
neurons-own[start nend my-layer my-row allele activation-last-time] ; layer = column

; k-tournament > 2: le cas qui fache: si deux paths sont winning, mais differents, mais qu'ils partagent quelques aretes. Dans ce cas les chemins winnings (partagés entre winners mais pas avec losers) sont quand meme renforce, et les aretes partagees avec les losers ne le sont pas.
; k-tournament = 2 seems to be optimal



;=======================
;         SETUP
;=======================

; PEA can be initiated any way you like, the minimum required being:
; - a network of nodes and edges.
; - network must start with a start node, and all paths should end to an end node.
; - nodes must store: layer they are in and their allele value
; - links must store: their weight (probability)

to setup
  clear-all
  if not debug [no-display]
  setup-globals
  setup-phenotypes
  setup-network
  setup-neurons-phenotypes
  setup-mga
  if visu? [
    neurons-visu
    synapses-visu
    display
  ]
  reset-ticks
end

to setup-globals
  set paths-visu-height 0.2
  set layers n-layers + 2 ; + 2 for start and end nodes

  set best-phenotypes nobody

  ; Plot
  set pea-max-moving 0
  set pea-max-moving-n 0
  set pea-avg-moving 0
  set pea-avg-moving-n 0
  set pea-moving-last-time 0

  set mga-max-moving 0
  set mga-max-moving-n 0
  set mga-avg-moving 0
  set mga-avg-moving-n 0
  set mga-moving-last-time 0
end

; Setup the phenotypes
; These are represented as a list of possible values
; Modify here if you want to add a new kind of phenotype (integer, float, command control, condition, etc..)
to setup-phenotypes
  ; By default, empty
  set phenotypes []

  ; Binary
  if phenotypes-kind = "Binary" [
    set phenotypes [0 1]
  ]

  ; Integer
  if phenotypes-kind = "Integer" [
    set phenotypes n-values 10 [?]
  ]
end

; Setup neurons/nodes
to setup-neuron [layer row]
  ; Init vars
  set start false
  set nend false
  set heading 90
  set my-layer layer
  set my-row row
  set allele false
  set activation-last-time 0

  ; Init shape
  set shape "circle 2"
  set color orange
end

; Create the initial links between the neurons (start is linked to all first nodes, all last nodes are linked to end, and any other node is linked to the next one in its own path, such that we have parallel paths)
to setup-synapses
  ask neurons with [start] [
    foreach agent-set-to-list neurons with [my-layer = 1] [
      let layer my-layer
      create-synapse-to ? [ set sy-layer layer ]
    ]
  ]
  ask neurons with [my-layer = layers - 2] [
    let layer my-layer
    create-synapse-to one-of neurons with [nend] [ set sy-layer layer ]
  ]
  ask neurons with [not start and not nend and not (my-layer = layers - 2)] [
    let layer my-layer
    let row my-row
    create-synapse-to one-of neurons with [my-layer = layer + 1 and row = my-row] [ set sy-layer layer ]
  ]

  ask synapses [
    set weight 1
    set color grey
  ]
  normalize-all-weights
end

; Initialize the allele values for each neuron
; Guaranties that any neuron has a unique allele in its layer
to setup-neurons-phenotypes
  if phenotypes-kind != "None" [
    ; Check that the number of rows is not above the number of different allele values, to assure that any neuron in a layer has a unique allele
    if n-rows > length phenotypes [
      error "Too many n-rows for this type of phenotypes! Please lower down the number to 2 maximum"
    ]

    ; For each layer/column
    let col n-values (layers - 2) [? + 1]
    foreach col [
      ; Randomize the phenotypes attribution
      let pool shuffle phenotypes ; pool = pool of phenotypes for this layer, thus we can just remove a value to make sure that it won't be assigned again to any other neuron
                                  ; Assign the phenotype to this neuron
      ask neurons with [my-layer = ?] [
        set allele first pool ; Assign phenotype
        set pool but-first pool ; And remove this value from the pool for this layer
      ]
    ]
  ]
end

; Setup the initial network following the description of the paper
; This will also ensures that it's nicely represented visually
; A network is a rooted acyclic graph
; Note: The graph is constrained so that paths are guaranteed to contain exactly one vertex for each distinct parameter.
to setup-network
  ; Init size vars
  set layer-width (world-width - 1) / layers ; width of one layer
  set network-height (world-height - 1) * (1 - paths-visu-height) ; height of the whole network
  set layer-height network-height ; height of one layer
  if n-rows > 1 [ ; if there is only one path, the row will be aligned with start and end neurons, else the paths will be shifted so that we clearly see the bifurcation
    set layer-height network-height / n-rows
  ]

  ; Creating start and end neurons
  let start-posx layer-width / 2
  let end-posx layer-width * (layers - 1) + layer-width / 2
  let se-posy network-height / 2
  create-neurons 1 [ ; start neuron
    setup-neuron 0 -1
    set start true
    setxy start-posx se-posy
  ]
  create-neurons 1 [ ; end neuron
    setup-neuron (layers - 1) -1
    set nend true
    setxy end-posx se-posy
  ]

  ; Creating other neurons
  let column n-values n-layers [? + 1]
  let row n-values n-rows [?]
  foreach column [ ; for each column (layer)
    let layer ?
    let posx layer * layer-width + (layer-width / 2) ; place in the middle of the current layer
    foreach row [ ; for each row (path)
      let posy ? * layer-height + (layer-height / 2) ; place in the middle of the height for the current path

      ; Create one neuron, init its vars and place it adequately
      create-neurons 1 [
        setup-neuron layer ?
        setxy posx posy
      ]
    ]
  ]

  ; Create links (synapses = directed edges) between neurons
  setup-synapses
end



;=======================
;         MAIN
;=======================

; Main procedure
to go
  if not debug [ no-display ] ; disable gui display to speedup processing, the time slider won't influence the setup procedure
  visualization-cleanup
  let reward-func "max"
  if reward-kind = "Min length" or reward-kind = "Guess-a-number" [ set reward-func "min" ]
  pea reward-func choose-reward
  if enable-mga [ mga reward-func choose-reward ]
  visualization-refresh
  tick
  if visu? [ display ] ; reenable gui display
end

; Select the correct reward to evaluate the paths
; Modify here if you want to add a new reward kind for your problem
to-report choose-reward
  if reward-kind = "Random" [
    report (task reward-random)
  ]

  if reward-kind = "Max sum" [
    report (task reward-max-sum)
  ]

  if reward-kind = "Count ones" [
    report (task reward-count-ones)
  ]

  if reward-kind = "Oscillating count" [
    report (task reward-oscillating)
  ]

  if reward-kind = "Disequilibrium separate" [
    report (task reward-disequilibrium-separate)
  ]

  if reward-kind = "Disequilibrium shift" [
    report (task reward-disequilibrium-shift)
  ]

  if reward-kind = "Guess-a-number" [
    report (task reward-guess-a-number)
  ]

  if reward-kind = "Min length" [
    report (task reward-min-length)
  ]
end


;--------
; REWARD
;--------
; Reward functions, modify or add your own depending on the problem you want to solve
; This is also here that you set the best-phenotypes to compute the hamming distance. Note that best-phenotypes must be a list of list (a list of phenotypes), just like the phenotypes, to be able to compare (NOT the value of the reward but the phenotype that should lead to the best reward!). Note: best-phenotypes is a list of best phenotypes, thus you can have multiple optimal phenotypes (or just one if you put one phenotype in a list).


; Random reward, will attribute a random reward for each path
to-report reward-random [Lphenotypes]
  ; First set the best phenotype to compute the hamming distance (which depends on the reward)
  if best-phenotypes = nobody [
    set best-phenotypes (list (n-values n-layers [0]))
  ]

  ; Then compute the reward
  let rewards []
  foreach Lphenotypes [
    set rewards lput (random 10) rewards
  ]
  report rewards
end

; Max sum reward, will attribute the most reward to the path which maximizes the sum of its alleles (eg: in binary, most 1's)
to-report reward-max-sum [Lphenotypes]
  ; First set the best phenotype to compute the hamming distance (which depends on the reward)
  if best-phenotypes = nobody [
    let max-pheno (max phenotypes)
    set best-phenotypes (list (n-values n-layers [max-pheno]))
  ]

  ; Then compute the reward
  let rewards []
  foreach Lphenotypes [
    set rewards lput (sum ?) rewards
  ]
  report rewards
end


; Count alleles = 1 in the path
to-report reward-count-ones [Lphenotypes]
  ; First set the best phenotype to compute the hamming distance (which depends on the reward)
  if best-phenotypes = nobody [
    set best-phenotypes (list (n-values n-layers [1]))
  ]

  ; Then compute the reward
  let rewards []
  foreach Lphenotypes [
    set rewards lput (length filter [? = 1] ?) rewards
  ]
  report rewards
end


; Circular oscillating count of one of the alleles.
; Example with binary: starts by counting alleles = 0 in the path at first, then after a few iterations it will switch and count alleles = 1, then switch again and count alleles = 0, etc.
to-report reward-oscillating [Lphenotypes]
  ; First set the best phenotype to compute the hamming distance (which depends on the reward)
  if best-phenotypes = nobody [
    set best-phenotypes (list (n-values n-layers [first phenotypes]))
    set reward-oscillating-best-allele first phenotypes
    set reward-oscillating-last-switch ticks
  ]

  ; Then compute the reward
  ; switch the counting goal every x ticks (count ones or zeros?)
  if ticks >= reward-oscillating-last-switch + reward-oscillation-time [
    set reward-oscillating-best-allele item (((position reward-oscillating-best-allele phenotypes) + 1) mod (length phenotypes)) phenotypes
    set reward-oscillating-last-switch ticks
    set best-phenotypes (list (n-values n-layers [reward-oscillating-best-allele])) ; change the best-phenotype
  ]
  let rewards []
  foreach Lphenotypes [
    set rewards lput (length filter [? = reward-oscillating-best-allele] ?) rewards
  ]
  report rewards
end

; Linkage disequilibrium demonstration: reward equally any path with identical phenotype (eg: all 0's, or all 1's in binary thus there would be two optimum)
to-report reward-disequilibrium-separate [Lphenotypes]
  ; First set the best phenotype to compute the hamming distance (which depends on the reward)
  if best-phenotypes = nobody [
    set best-phenotypes []
    foreach phenotypes [
      let pheno ?
      set best-phenotypes lput (n-values n-layers [pheno]) best-phenotypes
    ]
  ]

  ; Then compute the reward
  let rewards []
  foreach Lphenotypes [
    let path-pheno ?
    let reward []
    foreach phenotypes [
      let pheno ?
      set reward lput (length filter [? = pheno] path-pheno) reward ; could also use the hamming distance
    ]
    set rewards lput (max reward) rewards
  ]
  report rewards
end

; Linkage disequilibrium second demonstration: generate multiple optimal paths by shifting the alleles, eg in binary and 2 layers: [0 1] and [1 0]
to-report reward-disequilibrium-shift [Lphenotypes]
  ; First set the best phenotype to compute the hamming distance (which depends on the reward)
  if best-phenotypes = nobody [
    set best-phenotypes []
    let index-shift 0
    repeat (length phenotypes) [
      set best-phenotypes lput (n-values n-layers [item ((? + index-shift) mod (length phenotypes)) phenotypes]) best-phenotypes
      set index-shift index-shift + 1
    ]
  ]

  ; Then compute the reward
  let rewards []
  foreach Lphenotypes [
    set rewards lput (n-layers - (hamming-distances best-phenotypes ?)) rewards
  ]
  report rewards
end

to-report reward-guess-a-number [Lphenotypes]
  ; First set the best phenotype to compute the hamming distance (which depends on the reward)
  if best-phenotypes = nobody [
    let gnumber (random 256)
    user-message (word "Will now try to guess number: " gnumber)
    show (word "Decimal number to find: " gnumber " - binary: " (decimal-to-binary gnumber))
    if phenotypes-kind = "Binary" [ set best-phenotypes (list (decimal-to-binary gnumber) ) ]
    if phenotypes-kind = "Integer" [ set best-phenotypes (list (list gnumber)) ]
  ]

  ; Then compute the reward
  let rewards []
  if phenotypes-kind = "Binary" [
    foreach Lphenotypes [
      set rewards lput (hamming-distances-diff best-phenotypes ?) rewards
    ]
  ]
  if phenotypes-kind = "Integer" [
    foreach Lphenotypes [
      set rewards lput ((sum (item 0 best-phenotypes)) - (sum ?)) rewards
    ]
  ]
  report rewards
end

to-report reward-min-length [Lphenotypes]
  ; First set the best phenotype to compute the hamming distance (which depends on the reward)
  if best-phenotypes = nobody [
    set best-phenotypes (list (list 0) (list 1))
  ]

  ; Then compute the reward
  report map [length ?] Lphenotypes
end



;----------------------------------
; MICROBIAL GENETIC ALGORITHM (MGA)
;----------------------------------

to setup-mga
  set mga-pop n-values mga-pop-size [mga-random-phenotype]
end

to-report mga-random-phenotype
  report map [item ? phenotypes] (generate-random-indexes n-layers (length phenotypes))
  ;report n-values n-layers [item (random ((length phenotypes) - 1)) phenotypes]
end

to-report generate-random-indexes [list-length indexes-length]
  report n-values list-length [random indexes-length]
end

; Microbial Genetic Algorithm, a minimalistic genetic algorithm based on steady-state and infection (and a few other tricks) to reduce the size of the algorithm
; Can be done in a single line in C++ (read the article!)
; Ref: Harvey, I. (1996). The Microbial Genetic Algorithm.
to mga [reward-eval reward-callback]
  ; Tournament pair: 2 different at random
  let a-idx random ((length mga-pop) - 1)
  let b-idx a-idx
  while [b-idx = a-idx] [
    set b-idx random ((length mga-pop) - 1)
  ]
  let a item a-idx mga-pop
  let b item b-idx mga-pop
  ; Remove them from the pool (more convenient than tracking the indexes of winner and loser)
  set mga-pop (remove-item a-idx (remove-item b-idx mga-pop))

  ; Selection, with elitism for free
  let rewards (runresult reward-callback (list a b))
  ; Find the winner
  let win-idx 0
  ifelse reward-eval = "min"
  [ set win-idx arg-min rewards ]
  [ set win-idx arg-max rewards ]
  let W a
  let L b
  if win-idx = 1 [
    set W b
    set L a
  ]

  if debug [
    show W
    show L
  ]

  ; Plot fitness and Hamming distance
  if plot-visu? [
    ifelse reward-eval = "min"
    [ plot-fitness "mga" (min rewards) (mean rewards) ]
    [ plot-fitness "mga" (max rewards) (mean rewards) ]
    plot-hamming "mga" W
  ]

  ; Recombination
  ; Uniform crossover
  ifelse mga-crossover-uniform [
    ;let length-crossover mga-crossover * (length L)
    ;let indexes sublist (shuffle n-values (length L) [?]) 0 length-crossover ; optimized version: instead of doing a for loop and call random each time (which is expensive), just shuffle the list and truncate the number of alleles that are likely to be mutated. This reduces drastically the number of calls to random.
    ;foreach indexes [
    ;  set L replace-item ? L (item ? W)
    ;]

    ; Non-optimized version but less biased since even with a very low crossover rate it can still happen
    let length-crossover mga-crossover * (length L)
    let indexes n-values (length L) [?]
    foreach indexes [
      if random-float 1 <= mga-crossover [
        set L replace-item ? L (item ? W)
      ]
    ]
  ]
  ; Two-points crossover
  [
    let length-crossover mga-crossover * (length L)
    let random-start random ((length W) - length-crossover + 1)
    let part1 sublist L 0 random-start
    let part2 sublist W random-start (random-start + length-crossover)
    let part3 sublist L (random-start + length-crossover) (length L)
    set L (sentence part1 part2 part3)
  ]

  ; Mutation
  ;let length-mutation mga-mutation * (length L)
  ;let indexes sublist (shuffle n-values (length L) [?]) 0 length-mutation ; optimized version, same trick as for uniform crossover
  ;foreach indexes [
  ;  set L replace-item ? L (random ((length phenotypes) - 1))
  ;]
  ; Non-optimized version but less biased since even with a very low mutation rate it can still happen
  let length-crossover mga-crossover * (length L)
  let indexes n-values (length L) [?]
  foreach indexes [
    if random-float 1 <= mga-mutation [
      set L replace-item ? L item (random (length phenotypes)) phenotypes
    ]
  ]

  if debug [
    show W
    show L
    show "------"
  ]

  ; Place them back in the pool
  set mga-pop (lput W (lput L mga-pop))
end

to-report mga-cross [W L]
  let length-crossover mga-crossover * (length L)
  let random-start random ((length W) - length-crossover + 1)
  let part1 sublist L 0 random-start
  let part2 sublist W random-start (random-start + length-crossover)
  let part3 sublist L (random-start + length-crossover) (length L)
  set L (sentence part1 part2 part3)
  report L
end



;---------------
; PEA ALGORITHM
;---------------
; Read the article for full details on the algorithm and a neurally plausible implementation with Izhikevich neurons.
; To summary the PEA algorithm (which primary idea is to use paths in a network as individuals):
;= Initialize a network with (allele, layer number) on each node and a weight (probability) on each link (weights are uniform at the start: weight / number-of-out-links-from-current-neuron). Also the network must start with special nodes: start and end.
;= For each generation
; - Traverse the networks and generate two paths:
;   * use roulette wheel or any stochastic algorithm to walk the path randomly relatively to the links frequencies
;   * mutate a neuron you traverse with probability mu:
;     . generate randomly a new allele different from the current neuron from which we mutate
;     . if another neuron in this layer with the same allele exists, connect to it (this ensures that neurons have a unique allele in a layer, no duplicate allowed!)
;     . else create a new neuron in this layer with the generated, mutated allele
;     . link the mutated neuron with the parent (parent-current-neuron -> mutated-neuron with weight = init-weight) and the child of current neuron (mutated-neuron -> child-of-current-neuron with weight = 1), thus creating two links in total.
; - Tournament selection by using (on the phenotypes of the two paths) the reward/fitness function you want to determine which path is winning
; - Learning rule update:
;   * If it's a draw (no winning path, same reward) then do nothing
;   * Else:
;     . break down the links into three categories: winning (exclusive to winning path: on the winning path but not on the losing path); losing (exclusive to losing path); shared (present in both paths)
;     . apply learning rule:
;        > Winning links: weight = (1 + lambda) * weight
;        > Losing links: weight = (1 - lambda) * weight
;        > Shared links:
;           . if you want to maintain diversity (distinctive solutions), penalize them: weight = (1 - lambda) * weight
;           . else do nothing
; - Two-points crossover (once per generation):
;   * create a bridge between a neuron on winning path to a neuron in a greater layer in losing path (layer-losing > layer-winning) with weight = init-weight
;   * and inversely create a second bridge from losing path to winning path (but this bridge must be necessarily after the first bridge) with weight = init-weight
; - Clean up the network:
;   * Drop links if weight < weight-min-threshold
;   * (optional, but it helps a lot the convergence) Drop neurons if not traversed for period period-min-traversal
;= END
; NOTE: don't forget to normalize the weights (outflow/out-links of neurons) after:
; - learning rule update
; - mutation
; - crossover
; The update can either be done globally (naive) or locally (update only the outflow of the neurons where some out-links weights were changed, thus update only where it's required, this is more computationally efficient).
; NOTE2: init-weight is any value you want, a variable you can easily modify.
;

; PEA algorithm
to pea [reward-eval reward-callback] ; reward-eval = "max" or "min", it's how the winning reward will be evaluated ; reward-callback is the function that will be called to evaluate the rewards for the paths
  ; Path traversal on-demand for the tournament
  let paths []
  repeat k-tournament [
    set paths lput path-traversal paths
  ]

  ; Tournament selection and reward evaluation
  let tmp tournament-selection paths reward-eval reward-callback
  ; Continue only if there's a winner and a losing path (else if it's a draw, we skip without doing any change)
  ifelse tmp = nobody [
    ; No winner? then we color edges in grey
    if visu? [
      ask synapses [
        set color grey
      ]
    ]
  ]
  ; Else if there's a winner, we update the synapses' weights
  [
    let winning-paths item 0 tmp
    let lost-paths item 1 tmp

    ; Learning rule
    learning-rule-update winning-paths lost-paths

    ; Two-points cross-over, done once for each generation (according to crossover probability)
    path-cross-over winning-paths lost-paths

    ; Drop unused synapses (with a too low probability)
    ; This is important as it will allow to mutate and reinitiate a new synapse at the same place but with higher probability than if we had left this one in place. This promote continuous exploration.
    ; Particularly important for adaptive/dynamic environment where the reward function may change over time (even if the input the same)!
    drop-low-synapses

    ; Drop neurons that have not been activated since a long time
    ; Note: this is not necessary for the algorithm since creating mutated neurons is equivalent to connecting to unused neurons, thus we can just leave unused neurons alone and it would have the same result
    drop-inactivated-neurons
  ]

  ; Mutation, after a tournament just like the algoritm
;  let Lneurons extract-neurons-from-paths paths
;  let Lneurons-vrac sentence Lneurons
;  let mutate n-values (length Lneurons-vrac) [random-float 1.0]
;
;  set Lneurons-vrac shuffle Lneurons-vrac
;  foreach mutate [
;    ask item random (length Lneurons-vrac) [
;      mutate-neuron self
;    ]
;  ]

end

; Drop unused synapses (with a too low probability)
; This is important as it will allow to mutate and reinitiate a new synapse at the same place but with higher probability than if we had left this one in place. This promote continuous exploration.
; Particularly important for adaptive/dynamic environment where the reward function may change over time (even if the input the same)!
to drop-low-synapses
  if synapse-drop-below-weight > 0 [
    let to-normalize []
    ask synapses with [weight < synapse-drop-below-weight] [
      set to-normalize lput end1 to-normalize
      die
    ]
    if not empty? to-normalize [ normalize-selected to-normalize ]
  ]
end

; Drop neurons that have not been activated since a long time
; Note: this is not necessary for the algorithm since creating mutated neurons is equivalent to connecting to unused neurons, thus we can just leave unused neurons alone and it would have the same result
to drop-inactivated-neurons
  if neuron-drop-after-time > 0 [
    let to-normalize []
    ask (neurons with [not start and not nend]) with [(ticks - activation-last-time) >= neuron-drop-after-time] [ ; avoid start and end neurons
      ; drop all synapses
      ask my-in-synapses [
        set to-normalize lput other-end to-normalize
        die
      ]
      ask my-out-synapses [die]
      die ; drop the neuron
    ]
    if not empty? to-normalize [ normalize-selected to-normalize ]
  ]
end

; Extract utility functions
to-report extract-neurons-from-paths [paths]
  report map [item 0 ?] paths
end

to-report extract-synapses-from-paths [paths]
  report map [item 1 ?] paths
end

to-report extract-phenotypes [Lneurons]
  report remove-item 0 (remove-item ((length Lneurons) - 1) (map [[allele] of ?] Lneurons) ) ; remove start and end neurons
end

; Path traversal
; Generates a path by walking on the graph network, choosing synapses (directed edges) depending on the weights (probability) they have
to-report path-traversal
  let Lnodes [] ; list of neurons for this path
  let Lsynapses [] ; list of synapses for this path
  let currneuron one-of neurons with [start] ; Begin at the start neuron
  set Lnodes lput currneuron Lnodes ; add the start neuron to the list
  let proba 0.0 ; initial proba for the roulette-wheel sampling
  let parent-node nobody ; necessary for mutation
  let child-node nobody ; necessary for mutation

  ; Loop until we have not reached the end neuron
  while [not [nend] of currneuron] [
    ; Hop from the current neuron to the next one, randomly choosing the synapse but relatively to its weight using roulette-wheel
    ask currneuron [

      ; Update last activation time, to remove unused neurons if they aren't activated for a long enough time
      set activation-last-time ticks

      ; Make sure that there is at least one out-synapse, else we create one
      ; This may happen if we drop neurons that are unactivated for some amount of time: a neuron is dropped, but it was the only child of a parent!
      if not any? my-out-synapses [
        ifelse not variable-layers [
          ; Find the neuron at layer L + 1 of current neuron
          create-synapse-to one-of neurons with [my-layer = [my-layer + 1] of myself] [
            set sy-layer [my-layer] of myself
            set weight 1.0
            set color grey
          ]
        ]
        [
          ; Variable layers: find the first neuron in a subsequent layer from the current neuron (including the ending neuron whatever its layer is)
          create-synapse-to one-of neurons with [nend or my-layer > [my-layer] of myself] [
            set sy-layer [my-layer] of myself
            set weight 1.0
            set color grey
          ]
        ]
      ]

      ; List all out synapses from this neuron as a list
      ; We must use a list instead of an agentset to ensure that a fixed order will be used by the roulette-wheel (agentset have no order), else we can't get the post-synaptic neuron (or rather we will be picking a random one if we use an agentset, not the one selected by the roulette-wheel!).
      let outsyn (agent-set-to-list my-out-synapses)

      ; Choose a random synapse (and thus neuron) using roulette-wheel sampling
      ;let tmp roulette-wheel proba ([weight] of my-out-synapses) ; wrong version, the weights will be given to the roulette in a random order
      ;let tmp roulette-wheel proba (map [[weight] of ?] outsyn) ; correct version without exploration-rate (gamma)
      let m count my-out-synapses
      let tmp roulette-wheel proba (map [ [(weight + exploration-rate) / (1 + m * exploration-rate)] of ?] outsyn) ; correct version with exploration-rate (gamma). Note that if gamma = 0, then both correct versions are equivalent
      set proba item 0 tmp ; extract vars, this is the only way in NetLogo: by using temporary list
      let index item 1 tmp

      ; Get the post-synaptic neuron
      ask item index outsyn [
        set currneuron other-end ; Update the current neuron (walk along the path)
        set Lnodes lput currneuron Lnodes ; Add the neuron to the list of the path
        set Lsynapses lput self Lsynapses ; Add the synapse too
        set child-node other-end
      ]

      ; Mutate (only if not starting nor ending neuron)
      if not start and not nend and parent-node != nobody and child-node != nobody [
        if random-float 1.0 <= mutation [
          mutate-neuron parent-node self child-node
        ]
      ]

      set parent-node self ; for mutation. Careful! Set this only after having called Mutate, else before mutation, self is not the parent but the current node!
    ]
  ]

  ; Return the full path consisting of the neurons + synapses we walked
  report (list Lnodes Lsynapses)
end

; Mutate a neuron when walking on the graph (on-demand)
; This will create a new mutated node (with an allele than does not already exists in any neuron in this layer) and link it to the parent (neuron that activated current neuron) and to the child neuron (neuron that was activated by current neuron in the last path-traversal)
; Note: The graph is constrained so that paths are guaranteed to contain exactly one vertex for each distinct parameter.
to mutate-neuron [parent-node node child-node]
  ask node [
    ; Compute the row, will be used to define a new layer if variable-layers is enabled
    let layer my-layer
    if variable-layers [
      set layer layer + (floor random 2)
      ;if layer < 1 [ set layer 1 ] ; cannot be below or equal the start node (layer 0)
    ]
    let pool shuffle (remove allele phenotypes)
    let mutated? false
    foreach pool [
      let mu-node one-of neurons with [my-layer = layer and allele = ?]

      ; Create a new mutated neuron if it does not exists (mutated in the sense that it has a new allele value that was not yet available in any neuron in this layer)
      ifelse redundancy? or mu-node = nobody [
        ; Compute the y coordinate: try to place the mutated node above first, else place it below the current node
        let mu-layer-height min (list 5 (network-height / (length phenotypes))) ; height of one layer

        ; Try to place above if there's not any neuron at this place, else we place below (whether there is already a neuron or not)
        ; TODO: find a better, guaranteed way to place neurons without overlapping
        ;let posx xcor mod world-width
        let posx (layer * layer-width + (layer-width / 2) ) mod world-width
        let posy ycor + mu-layer-height
        if posy > (world-height - 1) or any? neurons with [abs (ycor - posy) <= 2 and abs (xcor - posx) <= 1] [ ; also we limit if it's outside world's boundaries
          set posy ycor - mu-layer-height
          if posy < 0 [ set posy posy + (world-height - 1) ]
        ]

        let row -1
        ; Create the new neuron
        hatch-neurons 1 [
          setup-neuron layer row
          setxy posx posy
          set allele ?
          set mu-node self
        ]

        ; Create a link from parent-node to mutated-node (bypassing current node)
        ask parent-node [
          let sylayer my-layer
          create-synapse-to mu-node [
            set sy-layer sylayer
            set weight synapse-initial-weight-on-mutation
            set color grey
          ]
        ]

        ; Create a link from mutated-node to the child-node of current node
        ask mu-node [
          let sylayer my-layer
          let c-node child-node
          if variable-layers [
            set c-node one-of neurons with [nend or my-layer > layer]
          ]
          create-synapse-to c-node [
            set sy-layer sylayer
            set weight 1
            set color grey
          ]
        ]

        ; Normalize the weights
        normalize-parent parent-node
        ; Stop here the foreach
        set mutated? true
      ]
      ; Else a neuron with this allele already exists
      [
        ; Node already exists but is not linked by the parent-node (may have been forgotten or just never linked from here), we link to it
        if not [out-synapse-neighbor? mu-node] of parent-node [
          ; Create a link from parent-node to mutated-node (bypassing current node)
          ask parent-node [
            let sylayer my-layer
            create-synapse-to mu-node [
              set sy-layer sylayer
              set weight synapse-initial-weight-on-mutation
              set color grey
            ]
            set mutated? true
          ]
        ]

        ; Create a link from mutated-node to the child-node of current node
        let c-node child-node
        ; If variable-layers, we will link to any neuron in a subsequent layer, not just the child-node of current node
        if variable-layers [
          set c-node one-of neurons with [nend or my-layer > layer]
        ]
        if not [out-synapse-neighbor? c-node] of mu-node [ ; only if the mutated-node is not already connected to the child-node
          ask mu-node [
            let sylayer my-layer
            create-synapse-to c-node [
              set sy-layer sylayer
              set weight 1
              set color grey
            ]
          ]
          set mutated? true
        ]

        ; Normalize the weights
        normalize-parent parent-node
        normalize-parent mu-node
      ]
      ; Stop here the foreach if node was mutated
      if mutated? [ stop ]
      ; Else node already exists and is already linked with both the parent-node and the child-node, so we will try another value from the pool
    ]
  ]
end

; Two-points cross-over between the winning path(s) and the lost path(s)
; This will create a bridge from the losing path to the winning path and back to the losing path. This allows to bypass a part of the losing path, such that the losing path may benefit from some parts of the winning path, and the winning path can explore some of the losing path to see if it gains a better reward maybe.
to path-cross-over [winning-paths lost-paths]
  ; Random crossover
  if random-float 1.0 <= crossover [
    ; Extract the neurons from paths, we only need them (not the synapses)
    let winning-nodes extract-neurons-from-paths winning-paths
    let lost-nodes extract-neurons-from-paths lost-paths

    ; Extract one random path of nodes from the list of paths
    let winner one-of winning-nodes
    let lost-full one-of lost-nodes

    ; Slice to remove useless nodes
    let lost sublist lost-full 1 ((length lost-full) - 1)
    let lost-nolast sublist lost-full 1 ((length lost-full) - 2) ; all but the start, end nodes and also the last normal node (because we the lost path will bridge at the layer L to the winner path at layer L+1, thus we must ensure there's at least a L+1 layer after the selected L layer!)
    set winner sublist winner 1 ((length winner) - 1) ; all but the start and end nodes

    if length lost-nolast < 1 or lost-nolast = winner or lost = winner [ stop ] ; stop here if the lost and winner paths are the same except for the last node in the lost path (because we need at least two nodes in the lost path to make a useful crossover, eg: lost = [start (neuron 2) (neuron 4) (neuron 18) end] winner = [start (neuron 2) (neuron 4) end], when we remove the start node and two last nodes from the lost path, the two paths contain the same neurons: [(neuron 2) (neuron 4)], thus it's useless to make a crossover, and anyway it will produce a bug because no lost edges will be found!)

    set lost list-to-agent-set lost
    set lost-nolast list-to-agent-set lost-nolast
    set winner list-to-agent-set winner

    ; Select where the bridges will be placed
    ; Be careful: these are the nodes, not the synapses (bridge) since synapses are what we will create!
    ;let lost-bridge-won-index floor random (length lost)
    ;let won-bridge-lost-index (floor random ((length winner) - (lost-index + 1))) + lost-index + 1
    let lost-bridge-won one-of lost-nolast ; first select the bridge from losing path to winning path
    let won-bridge-lost one-of (winner with [my-layer > [my-layer] of lost-bridge-won]) ; then select the bridge from winning path to losing path so that this bridge is made after the first one

    if lost-bridge-won = nobody or won-bridge-lost = nobody [ stop ] ; if there's no possibility to crossover between the two paths (maybe they are of different lengths and thus it makes no sense to crossover from the longest path to the shortest one or whatever), we stop

    ;show (word (agent-set-to-list lost) " " (agent-set-to-list winner) " // " lost-bridge-won " " won-bridge-lost) ; debug

    ; Create the first bridge: between a lost neuron at layer L to a winning neuron at layer L + 1 (or L + x if variable-layers is enabled)
    ask lost-bridge-won [
      ; Find the winning neuron to which the losing neuron will be linked to
      let win-bridged nobody
      ifelse variable-layers
      [ set win-bridged one-of (winner with [my-layer > [my-layer] of myself and not in-synapse-neighbor? myself]) ] ; variable layers size: we can bridge towards any winning neuron in any layer as long as the layer is greater (after) than the layer of the losing (current) neuron
      [ set win-bridged one-of (winner with [my-layer = [my-layer + 1] of myself and not in-synapse-neighbor? myself]) ] ; fixed layers size: we bridge necessarily to the winning neuron in the directly next layer

      ; If no neuron was found (either because there's no layer after when variable-layers is enabled, or either because this losing neuron is already linked to the winning neurons), then we do nothing
      if win-bridged != nobody [
        ; Create the bridge
        create-synapse-to win-bridged [
          set sy-layer [my-layer] of myself
          set weight synapse-initial-weight-on-mutation
          set color pink
        ]
        ; Renormalize weights
        normalize-parent self
      ]
    ]

    ; Create the second bridge: between a winning neuron at a later layer than the first bridge
    ask won-bridge-lost [
      ; Find the losing neuron to which we will link this winning neuron to
      let lost-bridged nobody
      ifelse variable-layers
      [ set lost-bridged one-of (lost with [my-layer > [my-layer] of myself and not in-synapse-neighbor? myself]) ]
      [ set lost-bridged one-of (lost with [my-layer = [my-layer + 1] of myself and not in-synapse-neighbor? myself]) ]

      ; If no neuron was found (either because there's no layer after when variable-layers is enabled, or either because this winning neuron is already linked to the losing neurons), then we do nothing
      if lost-bridged != nobody [
        ; Create the bridge
        create-synapse-to lost-bridged [
          set sy-layer [my-layer] of myself
          set weight synapse-initial-weight-on-mutation
          set color pink
        ]
        ; Renormalize weights
        normalize-parent self
      ]
    ]
  ]
end

; Selection by tournament competition (tournament: choose a pair of paths ; selection: select the winner)
; In the paper, it was a pairwise tournament (k = 2: comparing only two paths at a time), but here it is generalized to any k > 1
to-report tournament-selection [paths reward-eval reward-callback]
  ; Extract phenotypes for reward and synapses for discriminating paths
  let Lneurons extract-neurons-from-paths paths
  ;let Lsynapses extract-synapses-from-paths paths
  let Lphenotypes map [extract-phenotypes ?] Lneurons

  ; Compute rewards
  let rewards (runresult reward-callback Lphenotypes)

  ; Plot fitness
  if plot-visu? [
    ifelse reward-eval = "min"
    [ plot-fitness "pea" (min rewards) (mean rewards) ]
    [ plot-fitness "pea" (max rewards) (mean rewards) ]
  ]

  ; Find winner(s) (and losers at the same time)
  let winners []
  ifelse reward-eval = "min"
  [ set winners arg-mins rewards ]
  [ set winners arg-maxs rewards ] ; get the paths which maximizes the reward

  ; Plot Hamming distance for one of the winners
  if plot-visu? [ plot-hamming "pea" (item (one-of winners) Lphenotypes) ]

  ; If all paths are equally winners (or losers), it's a draw, we don't do anything
  if length winners = k-tournament [report nobody]

  ; For each winning path, add it into the list of winning-paths and remove it from the lost-paths
  let winning-paths []
  let lost-paths paths
  foreach (sort-by > winners) [ ; must sort so that we remove in lost-paths in the correct order, else if the indexes are shuffled we may incorrectly remove a path because the index has shifted with previous removes!
    set winning-paths lput (item ? paths) winning-paths
    set lost-paths remove-item ? lost-paths
  ]

  report (list winning-paths lost-paths)
end

; Learning rule
; Updates the graph network depending on the winning path(s) and losing path(s) from the last path traversal and tournament selection
to learning-rule-update [winning-paths lost-paths]
  ; Extract the synapses from the paths
  let winning-syn list-to-agent-set (flatten-list extract-synapses-from-paths winning-paths)
  let lost-syn list-to-agent-set (flatten-list extract-synapses-from-paths lost-paths)

  ; Discriminate synapses according to their status: won (on the winning path(s) and not shared by losing path(s)) - shared (between winning and losing path(s)) - lost (else)
  let syn-shared intersection winning-syn lost-syn
  let syn-won difference winning-syn lost-syn
  let syn-lost difference lost-syn winning-syn

  ; Update each synapse only once even if it's in multiple different paths?
  ;set syn-shared list-to-agent-set remove-duplicates (agent-set-to-list syn-shared)
  ;set syn-won list-to-agent-set remove-duplicates (agent-set-to-list syn-won)
  ;set syn-lost list-to-agent-set remove-duplicates (agent-set-to-list syn-lost)

  ; Update weights rule (learning rule)
  let to-normalize [] ; neurons whose synapses' weight were modified, and thus have to be renormalized (the synapses, not the neurons)
  ; Increase weight of won paths (on the winning path(s) and not shared by losing path(s))
  ask syn-won [
    set weight weight * (1 + learning-rate)
    set to-normalize lput end1 to-normalize ; add the pre-synaptic neuron (neuron that is the parent of this synapse) to the list of nodes to which the outflow must be normalized
    if visu? [ set color green ] ; visualization: color the synapses
  ]
  ; Decrease weight of lost paths (on the losing path(s) and not shared by the winning path(s))
  ask syn-lost [
    set weight weight * (1 - learning-rate)
    set to-normalize lput end1 to-normalize
    if visu? [ set color red ] ; visualization: color the synapses
  ]
  ; Maybe decrease weight of shared paths (shared between winning and losing path(s)) if we want to promote diversity
  ask syn-shared [
    set weight weight * (1 - diversity-rate)
    set to-normalize lput end1 to-normalize
    if visu? [ set color blue ] ; visualization: color the synapses
  ]

  ; Normalize the outflow of all outflow synapses for neurons where some weights were changed
  normalize-selected to-normalize
end

; Normalize all weights of all synapses at once, whether they were modified or not.
; This is the most rigorous approach (keep correct probabilities) but it's also the most consuming since all synapses will be updated, not just the one that are in the last modified paths.
to normalize-all-weights
  ask neurons [
    let wsum 0
    ask my-out-synapses [
      set wsum wsum + weight
    ]
    ask my-out-synapses [
      set weight weight / wsum
    ]
  ]
end

; Local normalization to be called each time a synapse's weight is updated, a lot less CPU intensive since it updates only the required weights
; But it can be called multiple times for the same neuron if multiple synapses of the same neuron are updated (eg: one for winning path and another one for losing path).
; Also it's not very rigorous since the first weights to be updated will weight less in the end (probabilities should be updated all at once for both winning and losing path before normalization).
; Anyway this is the implementation given in the paper.
to normalize-parent [parent-neuron]
  ask parent-neuron [
    let wsum 0
    ask my-out-synapses [
      set wsum wsum + weight
    ]
    ask my-out-synapses [
      set weight weight / wsum
    ]
  ]
end

; Local normalization to be called with the list of neurons that had some synapses updated.
; This is the most optimized version, normalizing synapses only once and only when required.
; It's is also rigorous in terms of probabilities.
to normalize-selected [paths-neurons]
  if is-list? paths-neurons [
    set paths-neurons list-to-agent-set paths-neurons
  ]
  ask paths-neurons [
    let wsum 0
    ask my-out-synapses [
      set wsum wsum + weight
    ]
    ask my-out-synapses [
      set weight weight / wsum
    ]
  ]
end



;=======================
;     VISUALIZATION
;=======================

; Cleanup some visualization stuff
to visualization-cleanup
  if visu? [
    ; Need to reset color to grey to not mixup winning paths from different iterations
    ask synapses [
      set color grey
    ]
  ]
end

; Refresh some general visualizations on the network
to visualization-refresh
  if visu? [
    neurons-visu ; refresh neurons visu (allele, etc.)
    synapses-visu ; refresh synapses weights
  ]
end

; Show allele phenotype on the neurons
to neurons-visu
  ask neurons [
    ifelse allele = false [
      if start [
        set label "S"
      ]
      if nend [
        set label "F"
      ]
    ]
    [
      set label allele
    ]
  ]
end

; Show weights on the synapses
to synapses-visu
  ask synapses [
    set label (precision weight 4)
  ]
end

to plot-fitness [pen f-max f-avg]
  if pen = "pea" [
    set pea-max-moving (pea-max-moving * pea-max-moving-n + f-max) / (pea-max-moving-n + 1)
    set pea-max-moving-n pea-max-moving-n + 1

    set pea-avg-moving (pea-avg-moving * pea-avg-moving-n + f-avg) / (pea-avg-moving-n + 1)
    set pea-avg-moving-n pea-avg-moving-n + 1

    if ticks >= (pea-moving-last-time + moving-avg-time) [
      set-current-plot "Fitness"

      set-current-plot-pen (word pen "-max")
      plot pea-max-moving

      set-current-plot-pen (word pen "-avg")
      plot pea-avg-moving

      set pea-max-moving 0
      set pea-max-moving-n 0
      set pea-avg-moving 0
      set pea-avg-moving-n 0
      set pea-moving-last-time ticks
    ]
  ]

  if pen = "mga" [
    set mga-max-moving (mga-max-moving * mga-max-moving-n + f-max) / (mga-max-moving-n + 1)
    set mga-max-moving-n mga-max-moving-n + 1

    set mga-avg-moving (mga-avg-moving * mga-avg-moving-n + f-avg) / (mga-avg-moving-n + 1)
    set mga-avg-moving-n mga-avg-moving-n + 1

    if ticks >= (mga-moving-last-time + moving-avg-time) [
      set-current-plot "Fitness"

      set-current-plot-pen (word pen "-max")
      plot mga-max-moving

      set-current-plot-pen (word pen "-avg")
      plot mga-avg-moving

      set mga-max-moving 0
      set mga-max-moving-n 0
      set mga-avg-moving 0
      set mga-avg-moving-n 0
      set mga-moving-last-time ticks
    ]
  ]
end

to plot-hamming [pen pheno]
  let hamming-dist (hamming-distances-diff best-phenotypes pheno)

  set-current-plot "Hamming distance"
  set-current-plot-pen pen
  plotxy ticks hamming-dist
end

; Compute the hamming distance between a list L and a target list Lbest
to-report hamming-distance [Lbest L]
  report sum (map [ifelse-value (?1 = ?2) [0] [1]] L Lbest) ; +1 if character is different, else 0. At the end we do the sum.
end

; Compute the hamming distance between a list L and a target list Lbest
; This version works on strings of different length by comparing the characters of the substrings with length = smallest string length, and then the rest of the longer string is added (just like if the added characters are addition to the smallest string)
to-report hamming-distance-diff [Lbest L]
  ifelse length L = length Lbest
  [ report hamming-distance Lbest L ]
  [
    let minlength (min (list (length L) (length Lbest)))
    let maxlength (max (list (length L) (length Lbest)))
    report (maxlength - minlength) + (hamming-distance (sublist Lbest 0 minlength) (sublist L 0 minlength))
  ]
end

; Compute the hamming distance between a list L against a megalist megaLbest containing several target lists Lbest. The minimum (best) hamming distance will be returned.
to-report hamming-distances [megaLbest L]
  report min (map [hamming-distance ? L] megaLbest)
end

to-report hamming-distances-diff [megaLbest L]
  report min (map [hamming-distance-diff ? L] megaLbest)
end


;=======================
;     AUX FUNCTIONS
;=======================

; Recursively print a mega list (recursive list of listes or in fact any other element)
; Should be used with print, eg: print list:pretty-print-recursive (list (list 1 2) (list 3 4))
to-report list:pretty-print-recursive [megalist]
  report list:pretty-print-recursive-aux megalist 1
end

to-report list:pretty-print-recursive-aux [megalist level] ; megalist is a list of listes
  ifelse is-list? megalist and (is-list? item 0 megalist or is-agent? item 0 megalist) [
    let n (length megalist)
    let i 0
    let text ""
    let levelstring (word (reduce [word ?1 ?2] (n-values level ["-"])) ">")
    while [i < n] [
      set text (word text "\n" (word levelstring i ": " (list:pretty-print-recursive-aux (item i megalist) (level + 1))))
      set i (i + 1)
    ]
    report text
  ]
  [
    if megalist = nobody [
      report "nobody"
    ]
    report megalist
  ]
end

; Retourne l'ensemble des patchs entre une cible et soi
; beam-width permet de definir la largeur du trait de patchs qu'on va checker
to-report in-between [objects target beam-width]
  let dist distance target
  report objects with [distance target + distance myself < dist + beam-width]
end
to-report patches-in-between [target beam-width]
  report in-between patches target beam-width
end

; Faster implementation of in-radius-nowrap than NetLogo native
to-report fast-in-radius-nowrap [objects radius]
  report objects with [distance-nowrap myself < radius]
end

to-report intersection [set1 set2]
  if is-agent? set2 [ set set2 (list set2) ] ; if it's only a single agent and not an agent-set, we convert to a list (not to an agent-set because we don't know if it's a turtle-set, patch-set or link-set)
  report set1 with [member? self set2]
end

to-report union [set1 set2]
  report (patch-set set1 set2)
end

to-report difference [set1 set2]
  if is-agent? set2 [ set set2 (list set2) ] ; if it's only a single agent and not an agent-set, we convert to a list (not to an agent-set because we don't know if it's a turtle-set, patch-set or link-set)
  report set1 with [not member? self set2]
end

to-report exclusive [set1 set2]
  if is-agent? set1 [ set set1 (list set1) ] ; if it's only a single agent and not an agent-set, we convert to a list (not to an agent-set because we don't know if it's a turtle-set, patch-set or link-set)
  if is-agent? set2 [ set set2 (list set2) ]
  report (patch-set (set1 with [not member? self set2]) (set2 with [not member? self set1]))
end

to-report agent-set-to-list [agent-set]
  ;report sort agent-set ;; list sorted by default agent order (e.g. by who-number)
  ;report sort-by [ criteria ] agent-set ;; list sorted by user-defined order
  if not any? agent-set [ report [] ]
  report [ self ] of agent-set ;; list, in randomized order
end

to-report list-to-agent-set [L]
  if is-list? L and empty? L [ report no-turtles ]
  ifelse is-turtle? (item 0 L)
  [ report turtles with [member? self L] ]
  [
    ifelse is-link? (item 0 L)
    [ report links with [member? self L] ]
    [ report patches with [member? self L] ]
    ]
end

; Get all max values in a list
to-report maxs [L]
  ;set L sort-by > L
  let rmax max L
  report filter [? = rmax] L
end

; Get the index of the first max value in a list
to-report arg-max [L]
  report position (max L) L
end

; Get the indexes of all max values in a list
to-report arg-maxs [L]
  let rmax max L ; get the max value
  let indexes n-values (length L) [?] ; precache the list of indexes of the list
  report filter [? >= 0] (map [ifelse-value (?1 = rmax) [?2][-1]] L indexes) ; for each item in the list that has the max value, we return the index (map index if value = max, else -1 ; then filter out indexes that are -1 since they are not possible)
end

; Get all min values in a list
to-report mins [L]
  ;set L sort-by < L
  let rmin min L
  report filter [? = rmin] L
end

; Get the index of the first min value in a list
to-report arg-min [L]
  report position (min L) L
end

; Get the indexes of all min values in a list
to-report arg-mins [L]
  let rmin min L
  let indexes n-values (length L) [?]
  report filter [? >= 0] (map [ifelse-value (?1 = rmin) [?2][-1]] L indexes)
end

; Check if a value is unique in a list
to-report unique-value? [value L]
  report ifelse-value ((length filter [? = value] L) > 1) [false][true]
end

; Roulette-wheel random selection aka stochastic sampling with replacement [Bak87]
to-report roulette-wheel [proba L]
  ;let pmax max L
  ;set proba proba + random-float 1.0 * 2 * pmax
  set proba proba + random-float 1.0
  let Lsize (length L)
  let index random Lsize
  while [proba > (item index L)] [
    set proba proba - (item index L)
    set index (index + 1) mod Lsize
  ]
  report (list proba index)
end

; Recursively extract all items inside a list and flatten the list onto one level
to-report flatten-list [L]
  let ret []
  foreach L [
    ifelse is-list? ? [
      set ret (sentence ret (flatten-list ?))
    ]
    [
      set ret lput ? ret
    ]
  ]
  report ret
  ;report map [ifelse-value (is-list? ?) [flatten-list ?] [?]] L
end

to print-agentset [agentset]
  print [self] of agentset
end

to-report binary-to-decimal [L]
  let powers n-values (length L) [?]
  report sum (map [(2 ^ ?2) * ?1] L powers)
end

to-report decimal-to-binary [number]
  let L []
  while [number > 0] [
    set L lput (number mod 2) L
    set number (floor (number / 2))
  ]
;  ifelse number = 1
;  [ set L lput 1 L ]
;  [ set L lput 0 L ]
  report L
end
@#$#@#$#@
GRAPHICS-WINDOW
350
10
789
470
-1
-1
13.0
1
14
1
1
1
0
0
0
1
0
32
0
32
0
0
1
ticks
30.0

SLIDER
5
120
177
153
learning-rate
learning-rate
0
1
0.1
0.1
1
NIL
HORIZONTAL

SLIDER
5
150
177
183
mutation
mutation
0
1
0.01
0.01
1
NIL
HORIZONTAL

SLIDER
5
180
177
213
crossover
crossover
0
1
0.1
0.1
1
NIL
HORIZONTAL

INPUTBOX
130
50
180
110
n-rows
2
1
0
Number

SWITCH
5
400
95
433
visu?
visu?
0
1
-1000

SLIDER
5
210
177
243
epsilon
epsilon
0
1
0
0.1
1
NIL
HORIZONTAL

BUTTON
285
35
345
68
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
285
105
345
138
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

INPUTBOX
80
50
130
110
n-layers
2
1
0
Number

SWITCH
5
430
95
463
debug
debug
1
1
-1000

CHOOSER
190
50
282
95
phenotypes-kind
phenotypes-kind
"None" "Binary" "Integer"
1

SLIDER
5
240
177
273
diversity-rate
diversity-rate
0
1
0.1
0.1
1
NIL
HORIZONTAL

BUTTON
285
70
345
103
1-step
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
190
95
282
140
reward-kind
reward-kind
"Random" "Max sum" "Count ones" "Oscillating count" "Disequilibrium separate" "Disequilibrium shift" "Guess-a-number" "Min length"
1

SLIDER
5
305
220
338
synapse-initial-weight-on-mutation
synapse-initial-weight-on-mutation
0
1
0.1
0.05
1
NIL
HORIZONTAL

INPUTBOX
5
340
145
400
synapse-drop-below-weight
0.05
1
0
Number

INPUTBOX
145
340
265
400
neuron-drop-after-time
0
1
0
Number

INPUTBOX
5
50
80
110
k-tournament
2
1
0
Number

SWITCH
15
505
135
538
variable-layers
variable-layers
1
1
-1000

PLOT
790
75
1120
270
Fitness
Generation
Reward
0.0
5.0
0.0
1.0
true
true
"" ""
PENS
"pea-max" 1.0 0 -2674135 true "" ""
"pea-avg" 1.0 0 -13840069 true "" ""
"mga-max" 1.0 0 -955883 true "" ""
"mga-avg" 1.0 0 -11221820 true "" ""

INPUTBOX
885
10
975
70
moving-avg-time
50
1
0
Number

PLOT
790
270
1120
435
Hamming distance
Generation
Distance
0.0
5.0
0.0
1.0
true
true
"" ""
PENS
"pea" 1.0 0 -2674135 true "" ""
"mga" 1.0 0 -13345367 true "" ""

INPUTBOX
190
140
310
200
reward-oscillation-time
200
1
0
Number

SWITCH
230
240
340
273
enable-mga
enable-mga
0
1
-1000

SWITCH
790
10
885
43
plot-visu?
plot-visu?
0
1
-1000

PLOT
790
435
1120
585
Network size
NIL
NIL
0.0
2.0
0.0
1.0
true
true
"" ""
PENS
"edges" 1.0 0 -2674135 true "" "plot count synapses"
"nodes" 1.0 0 -16777216 true "" "plot count neurons with [not start and not nend]"

SLIDER
5
270
177
303
exploration-rate
exploration-rate
0
1
0
0.1
1
NIL
HORIZONTAL

SLIDER
230
270
340
303
mga-crossover
mga-crossover
0
1
0.5
0.1
1
NIL
HORIZONTAL

SLIDER
230
300
340
333
mga-mutation
mga-mutation
0
1
0.1
0.1
1
NIL
HORIZONTAL

TEXTBOX
230
205
345
226
Microbial GA
20
0.0
1

TEXTBOX
205
10
325
35
Global param
20
0.0
1

TEXTBOX
45
10
150
31
PEA param
20
0.0
1

TEXTBOX
1010
20
1080
65
Plots & stats
20
0.0
1

INPUTBOX
265
330
340
390
mga-pop-size
50
1
0
Number

SWITCH
190
405
350
438
mga-crossover-uniform
mga-crossover-uniform
1
1
-1000

SWITCH
15
550
125
583
redundancy?
redundancy?
1
1
-1000

TEXTBOX
15
475
295
521
PEA extensions (experimental)
20
0.0
1

TEXTBOX
135
550
375
605
<- removes the unique phenotype constraint (eg: you can have multiple neurons with phenotype 0 on the same layer, creating multiple \"redundant\" paths which will produce the same phenotype)
11
0.0
1

TEXTBOX
140
500
435
556
<- allows neurons and synapses to create shorter or longer paths than n-layers (by adding or shortcircuiting layers). Use with reward-kind \"Guess-a-number\" or \"Min-length\".
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

Pathway Evolution Algorithm (PEA) by Chrisantha Fernando et al. (2011), implemented in NetLogo 5.0.5 by Stephen Larroque (2014) and updated to work with NetLogo v5.3.1 (2017).

PEA is a general, minimalist genetic algorithm but with the specificity that the units of evolution are pathways in a network instead of individuals of a population.

Also implemented Microbial Genetic Algorithm by I. Harvey (1996) as a comparison (and because PEA is a special case of MGA). MGA is a minimalist genetic algorithm that implements the basic tenets of evolutionary algorithm and can be programmed in one line of C code.

The MGA algorithm was implemented to compare the phenomenons that appears with PEA.

PEA is NOT meant to be more performant than MGA (even if that's often the case) or any other algorithm, but mostly to demonstrate that some neurobiological phenomenons can naturally emerge from paths in a network, and PEA demonstrate that with a very simple genetic algorithm (since PEA is really just a MGA applied to paths instead of individuals).

## HOW IT WORKS

Basically it's a minimalistic genetic algorithm (mutation, crossovers) but applied on paths instead of individuals: each path is considered as an individual, so any combination of paths in a network from the Start node to the End node is an individual.

At each iteration, a tournament is done to compare 2 randomly chosen paths, and the one that provides the best flow (ie, the best fit error) is reinforced, and the other one's weights are lowered.

The change of weights is a bit complicated to explain, but we change each link's weight of the winner path by comparing to the loser path: this allows to accelerate the best fit by also keeping what was good in the loser path. See the crossover rules from the code comments or the reference below.

See the sourcecode for a full algorithm summary in the comments.

## HOW TO USE IT

Please select:
* a phenotypes-kind = range of possible values for each node (None, Binary, Integer)
* a reward-kind = what is the objective to optimize (eg, max sum)
* n-layers = number of layers = number of nodes per path = length of phenotypes
* n-rows = number of distinctive phenotypes on each layer = number of "rows" of nodes

Example: with phenotypes-kind = Binary, with n-rows = 1 it will initialize with only one row of 0 and 1, if n-rows = 2 there will be two rows with 0 and 1 on each layer).

Then press Setup and Go to see the network expand, contract and learn how to best optimize the reward-kind you chose given the phenotypes-kind.

## THINGS TO NOTICE

There are a few neurobiological phenomenons that can be reproduced with this model, like the expansion/contraction (watch the number of nodes and edges in the chart at the right-side), memory (ability to reuse previously expanded edges and nodes, use reward-kind = "Oscillating count" and set neuron-drop-after-time = 0) and disequilibrium (ability to maintain several equivalent solutions, use reward-kind = Disequilibrium separate or shift and set neuron-drop-after-time = 0).

## THINGS TO TRY

You can try to play with the parameters (most specifically with synapse-initial-weight-on-mutation, synapse-drop-below-weight and neuron-drop-after-time).

You can also try to enable the (unofficial) extensions of the model:

* redundancy? will break the phenotype unicity constraint (in the base model, a neuron must have a unique phenotype in its layer: no other neuron in this layer must have the same phenotype), and thus allow for more redundancy since several paths can have the same phenotype. This is may be a fundamental property of neurobiological phenomenons (see Claude Berrou's Turbocodes and Gripon-Berrou Neural Networks, an extension of associative Hopfield memory model).
* variable-layers will allow the network to increase or reduce the number of layers as needed. This allows to find solutions that are beyond the initial n-layers number. Try this with the reward-kind = "Guess-a-number" or "Min length".

Variable-layers modify the original PEA algorithm's mutation process by adding the following rules:

1. the mutated neuron can be placed randomly on either the current neuron's layer L, or on L + 1 (this allows to increase the path size).

2. the mutated neuron can link to any neuron in a subsequent layer or to the end node directly (thus it doesn't just link to child-neuron at L+1 but any neuron > L). This allows to reduce the path size.

This strategy works well for "Min length" so 2. seems to be correct, but it does not fit well for "Guess-a-number" so the 1. needs another strategy.

## EXTENDING THE MODEL

* Add more plots.
* Enhance the variable-layers extension: it works well for "Min length" (reducing the network size) but not so well with "Guess-a-number" (increasing the size to reach the required number to guess).

## NETLOGO FEATURES

See the AUX FUNCTIONS part for generic NetLogo functions like ensemble operations, arg-max/arg-min and roulette-wheel.

## CREDITS AND REFERENCES

* Pathway Evolution Algorithm from: "Evolvable Neuronal Paths: A Novel Basis for Information and Search in the Brain", 2011, Fernando C, Vasas V, Szathmáry E, Husbands P. PLoS ONE 6(8): e23534. doi: 10.1371/journal.pone.0023534

Note that another implementation in Objective-C of this algorithm by the authors themselves is provided with the article.

* Microbial Genetic Algorithm from: "Microbial Genetic Algorithm", 1996, I. Harvey
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.3.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
1
@#$#@#$#@
