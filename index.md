---
created: 2025-05-18T01:07
---
# Predator-Prey Cellular Automaton Model

Based on: Cattaneo, Dennunzio, and Farina (2006) - "A Full Cellular Automaton to Simulate Predator-Prey Systems"

This implementation extends the original model to support:
- 1 prey species
- Up to 2 predator species (user-configurable)
- Fully local cellular automaton rules
- Configurable initial populations

---

## Table of Contents

1. [Introduction](#introduction)
2. [State Space](#state-space)
3. [Parameters](#parameters)
4. [Neighborhoods](#neighborhoods)
5. [Evolution Rules](#evolution-rules)
   - [Phase 1A: Attack Sub-Phase](#phase-1a-attack-sub-phase)
   - [Phase 1B: Reproduction Sub-Phase](#phase-1b-reproduction-sub-phase)
   - [Phase 2: Movement](#phase-2-movement)
6. [Enhanced Model](#enhanced-model)
7. [Complete Time Step Algorithm](#complete-time-step-algorithm)
8. [Implementation Notes](#implementation-notes)

---


## Introduction

<img src = "./assets/predator_prey_animation_20260218_193621.gif">Fig 1: Example image of the predator-prey CA model

The predator-prey model, describing the competition between predator and prey in an ecological system, is one of the most popular and widely studied biological models. One of the most fundamental ways to study such a system is by using the Lotka-Volterra equations, a set of non-linear ordinary partial differential equations. It assumes that, for a set of fixed positive constants A (the growth rate of prey), B (the rate at which predators destroy prey), C (the death rate of predators), and D (the rate at which predators increase by consuming prey), the following conditions hold :

1. A prey population x increases at a rate $dx=Axdt$ (proportional to the number of prey) but is simultaneously destroyed by predators at a rate $dx=-Bxydt$ (proportional to the product of the numbers of prey and predators).

2. A predator population y decreases at a rate $dy=-Cydt$ (proportional to the number of predators), but increases at a rate $dy=Dxydt$ (again proportional to the product of the numbers of prey and predators). 

 This gives the coupled differential equations

$$
\frac{dx}{dt} = Ax-Bxy
$$
and,
$$
\frac{dy}{dt} = -Cy+Dxy
$$

In this sort of model, the prey curve always lead the predator curve.

Critical points occur when $dx/dt=dy/dt=0$, such that the sole stationary point is therefore located at $(x,y)=(C/D,A/B)$

Traditionally, such systems have often been studied in a time-discretized manner using a lattice model evolved using Monte-Carlo methods to model the spatial evolution of the system. However, due to problems with Monte-Carlo evolution of such a system 
(such as non-local movement rules, lack exclusivity at lattice sites), CA (Cellular Automaton) models have become an increasingly popular choice to model such systems accurately, with local movement rules and lattice-site exclusivity. Here, we model such a predator-prey system using a CA framework. This also allows us to model an instance where we have 2 predators in a system hunting for the same kind of prey, something which is difficult to execute using a Monte Carlo framework.

The cellular automaton operates on a 2D rectangular lattice with **periodic boundary conditions** (toroidal topology). Each cell can contain at most one organism (**exclusion principle**).

At each discrete time step, the automaton applies two phases **synchronously** to all cells:
1. **Reaction Phase** (Attack → Reproduction)
2. **Movement Phase**

This is a **fully local** model - all interactions occur only between neighboring cells. There are no random "jumps" or Monte Carlo steps.

---

## State Space

Each cell (x, y) can be in one of the following states:

### Primary States (visible)
| State | Symbol | Description |
|-------|--------|-------------|
| **Empty** | `0a` | No organism present |
| **Prey** | `1` | Contains a prey organism |
| **Predator 1** | `2a` | Contains predator type 1 |
| **Predator 2** | `3a` | Contains predator type 2 (optional) |

### Temporary States (computational)
| State | Symbol | Description |
|-------|--------|-------------|
| **Empty after attack** | `0b` | Cell just emptied (prey was killed) |
| **Predator 1 fed** | `2b` | Predator 1 just successfully hunted |
| **Predator 2 fed** | `3b` | Predator 2 just successfully hunted |

**Note:** Temporary states only exist during computation and are converted back to primary states at the end of each time step.

---

## Parameters

### Prey Parameters
| Symbol | Name | Range | Description |
|--------|------|-------|-------------|
| **bp** | Prey birth probability | [0, 1] | Probability that one neighboring prey reproduces into an empty cell |
| **dp** | Prey death probability per predator | [0, 1] | Probability that one neighboring predator kills the prey |

### Predator Parameters
| Symbol | Name | Range | Description |
|--------|------|-------|-------------|
| **bh** | Predator birth probability | [0, 1] | Probability that a predator reproduces after successfully hunting |
| **dh** | Predator death probability | [0, 1] | Probability that a predator dies naturally each time step |

### Movement Parameter
| Symbol | Name | Range | Description |
|--------|------|-------|-------------|
| **r** | Movement sensing radius | [1, 3] | Radius of Moore neighborhood used for sensing (typical: 2) |

<!-- ### Typical Parameter Values
- **bp** = 0.6
- **dp** = 0.7
- **bh** = 0.3
- **dh** = 0.4
- **r** = 2

--- -->

## Neighborhoods

### Von Neumann Neighborhood
The **4 nearest neighbors** (N, S, E, W):
```
      N
      ↑
W ← (x,y) → E
      ↓
      S
```

Positions: `{(x, y-1), (x, y+1), (x-1, y), (x+1, y)}`

**Used for:** Attack phase, reproduction phase, movement execution

---

### Moore Neighborhood (Radius r)
All cells within distance r (using Chebyshev distance):

For radius r=2:
```
□ □ □ □ □
□ ■ ■ ■ □
□ ■ X ■ □
□ ■ ■ ■ □
□ □ □ □ □
```

**Used for:** Movement sensing, overcrowding detection

#### Moore Neighborhood Quadrants

The Moore neighborhood is divided into 4 quadrants for movement decisions:
```
     NORTH
       ↑
WEST ← X → EAST
       ↓
     SOUTH
```

**Quadrant formulas:**

For cell at position (x, y) with radius r:

- **North:** All cells `(x+j, y+i)` where `i ∈ [1, r]` and `j ∈ [-i, i]`
- **South:** All cells `(x+j, y-i)` where `i ∈ [1, r]` and `j ∈ [-i, i]`
- **East:** All cells `(x+i, y+j)` where `i ∈ [1, r]` and `j ∈ [-i, i]`
- **West:** All cells `(x-i, y+j)` where `i ∈ [1, r]` and `j ∈ [-i, i]`

---

## Evolution Rules

### Phase 1A: Attack Sub-Phase

This phase determines predator-prey interactions and immediate deaths.

---

#### Rule A1: PREY Cell Processing

**Input:** Cell in state `1` (PREY) at position (x, y)

**Step 1:** Count predators in Von Neumann neighborhood
```
npt(x, y) = number of cells in Von Neumann neighborhood with state ∈ {2a, 3a}
```

**Step 2:** Apply attack logic

**Case 1: No predators present** (`npt = 0`)
- If **basic model**: Cell remains `1` (PREY)
- If **enhanced model**: Apply overcrowding check (see Enhanced Model section)

**Case 2: Predators present** (`npt > 0`)

Calculate survival probability:
```
P(survival) = (1 - dp)^npt
```

**Interpretation:** 
- With 1 predator: survive with probability `(1 - dp)`
- With 2 predators: survive with probability `(1 - dp)²`
- With n predators: survive with probability `(1 - dp)^n`

**Generate random number** `U ~ Uniform(0, 1)`
```
If U < P(survival):
    new_state = 1 (PREY survives)
    If enhanced model: Apply overcrowding check
Else:
    new_state = 0b (PREY dies, cell becomes EMPTY_AFTER_ATTACK)
```

**Example:**
- dp = 0.7 (70% death probability per predator)
- npt = 2 predators
- P(survival) = (1 - 0.7)² = 0.09 = 9%
- P(death) = 91%

---

#### Rule A2: PREDATOR Cell Processing

**Input:** Cell in state `2a` or `3a` (PREDATOR) at position (x, y)

**Step 1:** Count prey in Von Neumann neighborhood
```
npr(x, y) = number of cells in Von Neumann neighborhood with state = 1
```

**Step 2:** Apply hunting logic

**Case 1: No prey present** (`npr = 0`)
```
new_state = 2a (PREDATOR, no change - hunt fails)
```

**Case 2: Prey present** (`npr > 0`)

Calculate hunt success probability:
```
P(hunt_success) = 1 - (1 - dp)^npr
```

**Interpretation:**
- With 1 prey: succeed with probability `dp`
- With 2 prey: succeed with probability `1 - (1-dp)²`
- With n prey: succeed with probability `1 - (1-dp)^n`

**Generate random number** `U ~ Uniform(0, 1)`
```
If U < P(hunt_success):
    new_state = 2b (PREDATOR_FED - successful hunt)
Else:
    new_state = 2a (PREDATOR - hunt failed)
```

**Example:**
- dp = 0.7 (70% success per prey)
- npr = 2 prey nearby
- P(hunt_success) = 1 - (1 - 0.7)² = 1 - 0.09 = 91%

---

#### Rule A3: EMPTY Cell Processing

**Input:** Cell in state `0a` (EMPTY)

**Action:** No change during attack phase
```
new_state = 0a (remains EMPTY)
```

---

### Phase 1B: Reproduction Sub-Phase

This phase handles births and natural deaths.

---

#### Rule R1: PREY Cell Processing

**Input:** Cell in state `1` (PREY)

**Action:** No change
```
new_state = 1 (remains PREY)
```

**Rationale:** Prey don't die during reproduction phase (only during attack). Prey reproduction happens in empty cells, not in occupied prey cells.

---

#### Rule R2: PREDATOR Cell Processing

**Input:** Cell in state `2a`, `2b`, `3a`, or `3b` (PREDATOR or PREDATOR_FED)

**Step 1:** Apply natural death

**Generate random number** `U ~ Uniform(0, 1)`
```
If U < dh:
    new_state = 0a (PREDATOR dies, cell becomes EMPTY)
Else:
    new_state = 2a (PREDATOR survives, normalized from 2a/2b to 2a)
```

**Note:** Both fed and unfed predators have the same death probability `dh`.

**Example:**
- dh = 0.4 (40% death probability)
- Each predator has 40% chance of dying, regardless of whether it ate

---

#### Rule R3: EMPTY Cell (Regular) Processing

**Input:** Cell in state `0a` (EMPTY) - was empty before attack phase

**Step 1:** Count prey and predators in Von Neumann neighborhood
```
npr(x, y) = number of cells in Von Neumann neighborhood with state = 1
npt(x, y) = number of cells in Von Neumann neighborhood with state ∈ {2a, 3a}
```

**Step 2:** Apply prey birth logic

**Case 1: Predators present OR no prey present**
```
If npt(x, y) > 0 OR npr(x, y) = 0:
    new_state = 0a (remains EMPTY)
```

**Rationale:** Prey don't reproduce near predators, and can't reproduce if no prey nearby.

**Case 2: Only prey present** (`npt = 0` AND `npr > 0`)

Calculate prey birth probability:
```
P(prey_birth) = 1 - (1 - bp)^npr
```

**Interpretation:**
- With 1 prey neighbor: birth with probability `bp`
- With 2 prey neighbors: birth with probability `1 - (1-bp)²`
- With n prey neighbors: birth with probability `1 - (1-bp)^n`

**Generate random number** `U ~ Uniform(0, 1)`
```
If U < P(prey_birth):
    new_state = 1 (PREY is born)
Else:
    new_state = 0a (remains EMPTY)
```

**Example:**
- bp = 0.6 (60% birth probability per neighbor)
- npr = 2 prey neighbors
- P(prey_birth) = 1 - (1 - 0.6)² = 1 - 0.16 = 84%

---

#### Rule R4: EMPTY_AFTER_ATTACK Cell Processing

**Input:** Cell in state `0b` (EMPTY_AFTER_ATTACK) - prey was just killed here

**Step 1:** Count fed predators in Von Neumann neighborhood
```
npt2(x, y) = number of cells in Von Neumann neighborhood with state ∈ {2b, 3b}
```

**Important:** Only count predators that successfully hunted (state 2b/3b), not all predators!

**Step 2:** Apply predator birth logic

**Case 1: No fed predators** (`npt2 = 0`)
```
new_state = 0a (becomes regular EMPTY)
```

**Case 2: Fed predators present** (`npt2 > 0`)

Calculate predator birth probability:
```
P(predator_birth) = 1 - (1 - bh)^npt2
```

**Interpretation:**
- With 1 fed predator: birth with probability `bh`
- With 2 fed predators: birth with probability `1 - (1-bh)²`
- With n fed predators: birth with probability `1 - (1-bh)^n`

**Generate random number** `U ~ Uniform(0, 1)`
```
If U < P(predator_birth):
    new_state = 2a (PREDATOR is born)
Else:
    new_state = 0a (becomes regular EMPTY)
```

**Example:**
- bh = 0.3 (30% birth probability per fed predator)
- npt2 = 2 fed predators nearby
- P(predator_birth) = 1 - (1 - 0.3)² = 1 - 0.49 = 51%

**Note:** If multiple predator types are present, priority can be given to one type (e.g., Predator 1 has priority over Predator 2).

---

### Phase 2: Movement

Movement is **fully local** - organisms only move to immediately adjacent cells (Von Neumann neighborhood).

---

#### Step M1: Determine Movement Intentions

For each cell containing an organism (PREY or PREDATOR):

**Substep M1.1:** Count organisms in Moore neighborhood quadrants

For cell at (x, y) with sensing radius r:
```
For each quadrant Q ∈ {North, South, East, West}:
    Count relevant organisms in quadrant Q
```

**Substep M1.2:** Determine preferred direction

**For PREY cells:**
```
Goal: Move away from predators

For each quadrant Q:
    predator_count[Q] = number of PREDATOR cells in quadrant Q

min_predators = minimum(predator_count values)

If min_predators = 0 AND all counts are 0:
    intention = NONE (no movement)
Else:
    candidate_quadrants = all quadrants with predator_count = min_predators
    chosen_quadrant = random_choice(candidate_quadrants)
    intention = Von Neumann neighbor in direction of chosen_quadrant
```

**For PREDATOR cells:**
```
Goal: Move toward prey

For each quadrant Q:
    prey_count[Q] = number of PREY cells in quadrant Q

max_prey = maximum(prey_count values)

If max_prey = 0:
    intention = NONE (no movement)
Else:
    candidate_quadrants = all quadrants with prey_count = max_prey
    chosen_quadrant = random_choice(candidate_quadrants)
    intention = Von Neumann neighbor in direction of chosen_quadrant
```

**Quadrant to direction mapping:**
```
North quadrant → Move to cell (x, y+1)
South quadrant → Move to cell (x, y-1)
East quadrant  → Move to cell (x+1, y)
West quadrant  → Move to cell (x-1, y)
```

**Example (PREY):**
```
Quadrant predator counts:
  North: 5
  South: 2  ← minimum
  East: 7
  West: 2   ← minimum (tied)

Candidates: {South, West}
Chosen: South (random)
Intention: Move to (x, y-1)
```

---

#### Step M2: Execute Movements with Conflict Resolution

**Substep M2.1:** Group intentions by target cell
```
For each organism with intention to move to target cell T:
    Add organism to list: target_cell_map[T]
```

**Substep M2.2:** Resolve conflicts
```
For each target cell T that is EMPTY:
    candidates = target_cell_map[T]
    
    If len(candidates) = 0:
        # No one wants this cell
        T remains EMPTY
    
    Else If len(candidates) = 1:
        # Only one organism wants this cell
        Move that organism to T
    
    Else:
        # Multiple organisms want this cell - CONFLICT
        chosen_organism = random_choice(candidates)
        Move chosen_organism to T
        All other candidates stay in their original positions
```

**Substep M2.3:** Handle non-moving organisms
```
For each organism that didn't move (either had no intention or lost conflict):
    Organism stays in its current cell
```

**Example:**
```
Cell (5, 5) is EMPTY

Intentions pointing to (5, 5):
  - Prey at (5, 6) wants to move to (5, 5)
  - Prey at (4, 5) wants to move to (5, 5)
  - Predator at (5, 4) wants to move to (5, 5)

Conflict! 3 organisms want same cell.

Resolution:
  random_choice selects Predator at (5, 4)
  
Result:
  - Predator moves from (5, 4) to (5, 5)
  - Prey at (5, 6) stays at (5, 6)
  - Prey at (4, 5) stays at (4, 5)
  - Cell (5, 4) becomes EMPTY
```

---

## Enhanced Model

The enhanced model adds **density-dependent mortality** for prey to create oscillatory dynamics similar to the logistic equation.

### Enhancement Rule: Prey Overcrowding

**When to apply:** 

This rule applies to PREY cells in two scenarios:
1. **No predators in Von Neumann neighborhood** (during attack phase)
2. **Prey survived attack** but predators were present (during attack phase)

**Rationale:** Even without predation, prey can die from overcrowding (competition for resources).

---

**Step E1:** Count prey in Moore neighborhood of radius r
```
For cell at (x, y):
    total_prey = 0
    For all cells (x', y') in Moore neighborhood of radius r:
        If cell (x', y') contains PREY:
            total_prey += 1
    
    total_cells = (2*r + 1)²
```

**Step E2:** Calculate enhancement function

Choose one of two functions:

**Cosine function:**
```
f(bp) = 1 - cos(π/2 × bp)
```

**Exponential function:**
```
f(bp) = 1 - exp(-e × bp)
```

**Properties:**
- Both functions map bp ∈ [0, 1] to output ∈ [0, 1]
- Both are monotonically increasing
- Cosine: smoother, more gradual
- Exponential: steeper, more pronounced effect

**Step E3:** Calculate overcrowding death probability
```
P(death_from_overcrowding) = (total_prey × f(bp)) / total_cells
```

**Interpretation:**
- Higher prey density → higher death probability
- Birth parameter bp influences mortality (higher reproduction potential → higher competition)
- Normalized by total cells in neighborhood

**Step E4:** Apply overcrowding death

**Generate random number** `U ~ Uniform(0, 1)`
```
If U < P(death_from_overcrowding):
    new_state = 0b (PREY dies from overcrowding)
Else:
    new_state = 1 (PREY survives)
```

---

**Example:**

Parameters:
- r = 2 (Moore radius)
- bp = 0.6
- Enhancement function: cosine

Scenario:
- Cell (x, y) contains PREY
- No predators in Von Neumann neighborhood
- Moore neighborhood (5×5 = 25 cells) contains 15 prey

Calculation:
```
total_prey = 15
total_cells = 25
f(0.6) = 1 - cos(π/2 × 0.6) = 1 - cos(0.942) = 1 - 0.588 = 0.412

P(death) = (15 × 0.412) / 25 = 6.18 / 25 = 0.247 = 24.7%
```

Result: 24.7% chance prey dies from overcrowding.

---

### Effect on Dynamics

Without enhancement:
- Prey population grows exponentially (limited only by grid size)
- No intrinsic carrying capacity
- Difficult to match Lotka-Volterra oscillations

With enhancement:
- Prey population follows logistic growth
- Carrying capacity emerges from overcrowding
- Oscillations similar to Lotka-Volterra equations
- Better match to differential equation models

---

## Complete Time Step Algorithm

### Input
- Current grid state at time t: `grid[t]`
- All parameters: bp, dp, bh, dh, r

### Output
- New grid state at time t+1: `grid[t+1]`

---

### Algorithm
```
FUNCTION evolve_one_step(grid[t]):

    // PHASE 1A: ATTACK
    grid_after_attack = empty grid
    
    FOR each cell (x, y) in grid[t]:
        state = grid[t][x, y]
        
        IF state = PREY:
            new_state = apply_rule_A1(x, y, grid[t])
        
        ELSE IF state = PREDATOR:
            new_state = apply_rule_A2(x, y, grid[t])
        
        ELSE IF state = EMPTY:
            new_state = apply_rule_A3(x, y, grid[t])
        
        grid_after_attack[x, y] = new_state
    END FOR
    
    
    // PHASE 1B: REPRODUCTION
    grid_after_reproduction = empty grid
    
    FOR each cell (x, y) in grid_after_attack:
        state = grid_after_attack[x, y]
        
        IF state = PREY:
            new_state = apply_rule_R1(x, y, grid_after_attack)
        
        ELSE IF state ∈ {PREDATOR, PREDATOR_FED}:
            new_state = apply_rule_R2(x, y, grid_after_attack)
        
        ELSE IF state = EMPTY:
            new_state = apply_rule_R3(x, y, grid_after_attack)
        
        ELSE IF state = EMPTY_AFTER_ATTACK:
            new_state = apply_rule_R4(x, y, grid_after_attack)
        
        grid_after_reproduction[x, y] = new_state
    END FOR
    
    
    // PHASE 2: MOVEMENT
    
    // Step M1: Determine intentions
    intentions = empty map
    
    FOR each cell (x, y) in grid_after_reproduction:
        state = grid_after_reproduction[x, y]
        
        IF state ∈ {PREY, PREDATOR}:
            target = determine_movement_intention(x, y, state, grid_after_reproduction)
            IF target ≠ NONE:
                intentions[(x, y)] = target
    END FOR
    
    
    // Step M2: Execute movements
    grid[t+1] = execute_movements(grid_after_reproduction, intentions)
    
    
    RETURN grid[t+1]

END FUNCTION
```

---

## Implementation Notes


### 1. Parameter Sensitivity

**Key relationships for stable dynamics:**

For coexistence (neither species goes extinct):
  - bp should be moderately high (0.5-0.7)
  - dp should be moderate (0.6-0.8)
  - bh should be lower than bp (0.2-0.4)
  - dh should be moderate to high (0.3-0.5)

For oscillations:
  - Use enhanced model
  - bp > dp × (1 - dh/bh) approximately
  
For prey-only survival:
  - Set bh very low or dh very high
  - Predators will go extinct

For predator-only survival:
  - Not possible (predators need prey)
  - Will always go extinct without prey


<!-- ---

### 6. Validation Against Lotka-Volterra

The CA should produce dynamics similar to discrete Lotka-Volterra equations:
```
x₁(t+1) = x₁(t)[(1 + k₁) - a₁₁x₁(t)] - a₁₂x₁(t)x₂(t)
x₂(t+1) = x₂(t)[(1 + k₂) - a₂₂x₂(t)] + a₂₁x₁(t)x₂(t)
```

**Approximate parameter mapping:**
- k₁ (prey growth) ≈ bp × 2
- k₂ (predator growth) ≈ -dh + bh × dp
- a₁₂ (predation rate) ≈ dp × 10
- a₂₁ (conversion) ≈ bh × dp × 5

These are rough approximations; exact calibration requires fitting.

--- -->

### 2. Common Issues and Solutions

**Issue:** Predators go extinct immediately
- **Cause:** bh too low or dh too high
- **Solution:** Increase bh or decrease dh

**Issue:** Prey go extinct immediately
- **Cause:** dp too high, bp too low, or too many initial predators
- **Solution:** Decrease dp, increase bp, or reduce predator density

**Issue:** No oscillations (populations stabilize)
- **Cause:** Enhanced model not enabled
- **Solution:** Enable enhanced model with overcrowding

**Issue:** Populations grow without bound
- **Cause:** Grid too large, no density dependence
- **Solution:** Enable enhanced model or reduce grid size

**Issue:** Movement seems random/uncoordinated
- **Cause:** Sensing radius r too small
- **Solution:** Increase r to 2 or 3

---

## Summary

This cellular automaton implements a **fully local, spatially explicit** predator-prey model with:

1. **Local interactions:** All processes occur between immediate neighbors
2. **Exclusion principle:** One organism per cell maximum
3. **Stochastic rules:** Probabilistic birth, death, and hunting
4. **Intelligent movement:** Based on local sensing of environment
5. **Density dependence:** Optional overcrowding mortality for realistic dynamics
6. **2 predators** can be used with different species parameters

The model produces population dynamics qualitatively similar to classical Lotka-Volterra differential equations while also providing **spatial information** about organism distribution and movement patterns.

---