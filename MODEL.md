# Process-Based Simulation Model of Konservat-Lagerstätten Formation

## Abstract

This model tests whether the preservation window in Burgess Shale–type Lagerstätten can emerge endogenously from coupled ecological–geochemical dynamics. A spatially explicit reaction–diffusion model of decay, oxygen consumption, and anaerobic biofilm growth demonstrates that intermediate oxygen resupply produces a self-organised anoxia field and realistic tissue selectivity, without prescribing the timing of preservation. The results support a conceptual shift from externally imposed preservation conditions toward **emergent kinetic arrest**.

---

## 1. Motivation

Mechanistic explanations for Konservat-Lagerstätten formation remain largely qualitative or reconstructed post hoc from preserved geochemical signatures. The framework of Gaines et al. (2012), characterised by oxidant limitation suppressing early decay, provides the most rigorous current account of Burgess Shale–type preservation. However, it does not explicitly model:

- how these conditions arise dynamically  
- why they develop at the correct rate relative to carcass decay  

This model addresses a complementary question:

> **Can the preservation window arise as an emergent property of coupled ecological–geochemical dynamics, rather than as a prescribed boundary condition?**

---

## 2. Core Hypothesis

Exceptional preservation requires **kinetic arrest within a critical intermediate window**:

- **Too early** → insufficient decay structuring  
- **Too late** → loss of morphological fidelity  

The central claim is that this window is **endogenously generated** through coupling between:

- organic decay  
- oxygen drawdown  
- anaerobic biofilm growth  

> Preservation fidelity depends on the *rate* of anoxia development relative to the decay trajectory, and this rate is itself an emergent system property.

---

## 3. Model Structure

### 3.1 Spatial Framework

The model operates on a 2D 100 × 100 grid with:

- Dirichlet boundary (O₂ = 1) at the sediment–water interface  
- Neumann (zero-flux) boundaries laterally and at depth  

Each grid cell tracks:

- **C** — organic concentration  
- **O₂** — porewater oxygen  
- **B** — anaerobic biofilm  
- **φ** — clay heterogeneity scalar  

A *Hallucigenia sparsa* body plan is embedded near the sediment–water interface  
(`hallucigenia_bodyplan.R`), selected due to well-documented taphonomic bias in preservation.

---

### 3.2 Tissue Classes

Four tissue types are defined with distinct decay constants:

| Tissue | Decay rate (k) |
|--------|----------------|
| Spines | 0.005 |
| Trunk | 0.015 |
| Lobopods | 0.030 |
| Soft tissue | 0.060 |

These values are plausible but currently unconstrained by experiment.

---

### 3.3 Emergent Kinetic Arrest

Anoxia is defined as a sigmoid function of oxygen concentration:


```
anox = 1 / (1 + exp(20 * (O₂ - 0.2)))
```

This field simultaneously:

- suppresses effective diffusion (pore sealing analogue)  
- inhibits aerobic decay  
- promotes anaerobic biofilm growth  

Critically, **anoxia is not imposed** — it arises from internal system dynamics.

---

### 3.4 Transport and Reaction Processes

Implemented in `run_model()` (see main R script):

- **Organic matter (C)** diffuses under zero-flux boundaries  
- **Oxygen (O₂)** diffuses with sub-stepping (numerical stability)  
- **Biofilm (B)** grows as a function of C, anoxia, and φ  

Oxygen consumption includes:

- aerobic respiration: `k_resp_aero * C * (1 - anox)`  
- anaerobic metabolism: `k_resp_anaer * B * anox`  
- **background sediment demand**: `k_sed`  

#### v10 Update: Sediment Oxygen Demand

A constant background sink is included:


```
k_sed = 0.002
```


This produces realistic downcore O₂ gradients independent of carcass processes and removes artefacts present in earlier versions.

---

### 3.5 Numerical Stability

Oxygen diffusion is sub-stepped to enforce the von Neumann stability condition:


```
D * dt ≤ 0.25
```

``


This prevents artificial oscillations and spurious anoxia at high diffusivity.

---

## 4. Experimental Design

The primary control parameter is **oxygen diffusivity (D_O₂)**:

| Condition     | D_O2 | Interpretation        |
|---------------|------|----------------------|
| Fast arrest | 0.05 | Rapid anoxia |
| Intermediate | 0.10 | Balanced regime |
| Slow arrest | 0.20 | Delayed anoxia |
| Oxic flush | 0.50 | Continuous oxygenation |

- 20 Monte Carlo replicates per condition  
- φ initialised as random heterogeneity  

### Fossil Preservation Index (FPI)


```
FPI = mean(C_final * anox_final) / mean(C_initial)
```

interpreted as the **fraction of initial preservation potential realised under anoxia**.

---

## 5. Results

### 5.1 Preservation Regimes

The model robustly produces:


Fast arrest > Intermediate > Slow arrest > Oxic flush

in both mean FPI and probability of successful preservation.

---

### 5.2 Tissue Selectivity

Under intermediate conditions:

- Spines show highest preservation  
- Soft tissues decay rapidly  
- Lobopods and trunk occupy intermediate space  

This reproduces observed Burgess Shale–type patterns.

---

### 5.3 Anoxia Dynamics

Anoxia develops as a **self-organised wave**:

- initiates at the carcass  
- propagates outward via coupled reaction–diffusion  

This is consistent with localised redox zonation inferred in Cambrian Lagerstätten.

---

### 5.4 Effect of Sediment O₂ Demand

Including `k_sed`:

- produces realistic vertical gradients  
- removes deep oxic artefacts  
- slightly elevates background anoxia  

Importantly, it **does not alter regime ordering**, indicating robustness.

---

## 6. Limitations

### 6.1 Clay Heterogeneity

φ is random rather than sedimentologically structured.

### 6.2 Biofilm Representation

Biofilm growth is phenomenological and does not distinguish metabolic pathways.

### 6.3 Tissue Decay Parameters

Decay constants are unconstrained by experimental datasets.

### 6.4 Dimensionality

Model is 2D and does not explicitly resolve vertical porewater structure.

### 6.5 Mineralisation

No explicit authigenic mineral formation (e.g. pyrite, phosphate) is modelled.

---

## 7. Next Steps

Priority developments:

1. Constrain decay parameters using experimental taphonomy datasets  
2. Introduce structured sediment fields (spatially correlated φ)  
3. Couple to mineralisation pathways (starting with pyritisation)  
4. Extend to 3D to resolve vertical redox gradients  

---

## 8. Conceptual Contribution

This model demonstrates that:

> **A preservation window and selective fidelity can emerge from internal dynamics without prescribing anoxia timing.**

The central proposal is **emergent kinetic arrest** as the mechanism controlling exceptional preservation.

---

## 9. Key Reference

Gaines et al. (2012)  
*Mechanism for Burgess Shale-type preservation*  
PNAS 109(14): 5180–5184

---

## 10. Code Structure

- `hallucigenia_bodyplan.R` — body plan generator  
- `run_model()` — core simulation  
- visualisation scripts — plotting outputs and diagnostics  

---

## 11. Contact

Billy Hunter  
AFBI Environment and Marine Science Division  

Code and further details available on request.


