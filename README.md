# Burgess Shale Taphonomic ABM

## Overview
This repository contains a spatially explicit reaction–diffusion model simulating microbial and geochemical controls on Burgess Shale–style soft-tissue preservation.

### Key Concept
The unique aspect of this model is its explicit coupling of **microbial processes with sedimentary geochemistry**. Rather than treating preservation as purely chemical or purely biological, the model captures the **feedbacks between microbial activity and the evolving geochemical environment** that collectively determine preservation outcomes.

In particular, the model represents how:
- Microbial metabolism reshapes local oxygen conditions
- Oxygen availability regulates microbial pathways and decay rates
- These feedbacks create emergent anoxic microenvironments within and around carcasses

This microbial–geological interaction is hypothesised to be a **primary control on exceptional fossil preservation**.

---

## Experiment: Oxygen Resupply and Preservation Potential

This experiment tests whether **rapid geochemical “lock-in” into anoxic conditions is necessary for exceptional preservation**, using a stylised *Hallucigenia* body plan embedded in sediment.

Four environmental regimes are compared:

- **Fast arrest**
- **Intermediate**
- **Slow arrest**
- **Oxic flush**

Each condition varies:
- Oxygen diffusivity (`D_O2`)
- Aerobic respiration rate (`k_resp_aero`)

to simulate different burial and ventilation scenarios.

---

## Microbial–Geochemical Coupling

The model explicitly links:

### Microbial processes
- Aerobic respiration (oxygen consumption tied to organic matter availability)
- Anaerobic biofilm growth under anoxic conditions
- Spatially heterogeneous microbial activity emerging from local conditions

### Geological / geochemical processes
- Oxygen diffusion from the water column
- Background sediment oxygen demand (`k_sed`)
- Sediment boundary conditions (oxic interface, zero-flux base)

### Feedbacks
- Microbial respiration drives oxygen depletion  
- Oxygen depletion promotes anoxia  
- Anoxia enhances biofilm growth and inhibits decay  
- These processes reinforce each other, generating **self-organised preservation niches**

**This feedback structure is the central innovation of the model and underpins all emergent preservation behaviour.**

---

## Key Model Features

- **Background sediment oxygen demand (`k_sed`)**
  - Constant zero-order sink representing ambient microbial respiration  
  - Produces realistic downcore O₂ gradients  

- **Sub-stepped O₂ diffusion**
  - Ensures numerical stability across regimes  

- **Normalised Fossil Preservation Index (FPI)**

```text
FPI = mean(C_final × anoxia_final) / mean(C_initial)
FPI = mean(C_final × anoxia_final) / mean(C_initial)

```

---
``


## Hallucigenia Body Plan Representation

The model uses a spatially explicit representation of a *Hallucigenia*-like organism to simulate anatomically structured decay and preservation.

The body plan is constructed on a 100 × 100 grid using the `hallucigenia_bodyplan.R` module, which generates:

- An initial organic carbon distribution (`C`)
- A tissue-specific decay rate matrix (`k_tissue`)
- A categorical tissue map for analysis and visualisation

### Morphology

The synthetic organism includes:

- **Trunk** — elongated, slightly curved elliptical body  
- **Head bulb** — anterior expansion  
- **Dorsal spines** — 7 pairs, elongated and outward-angled  
- **Ventral lobopods** — 7 pairs, shorter and broader appendages  
- **Internal soft tissue core** — central tract and pharyngeal region  

### Tissue Types and Decay Rates

Four tissue classes are defined, each with distinct decay constants:

| Tissue type | Description                        | Decay rate (k) |
|-------------|----------------------------------|----------------|
| Spine       | Sclerotised structures           | 0.005          |
| Trunk       | Cuticularised body wall          | 0.015          |
| Lobopod     | Appendages                       | 0.030          |
| Soft tissue | Internal organs / pharynx        | 0.060          |

This hierarchy reflects increasing lability from structural to soft tissues.

### Functional Role in the Model

The body plan is not just visual—it directly controls system dynamics:

- **Spatial heterogeneity in decay**  
  Different tissues degrade at different rates  

- **Variable oxygen consumption**  
  Organic matter distribution shapes local respiration  

- **Emergent preservation patterns**  
  Differential decay produces realistic fossilisation signatures  

Because the organism is explicitly resolved in space, the model captures how **anatomy interacts with microbial–geochemical processes**, leading to non-uniform preservation across the body.

### Usage

```r
source("hallucigenia_bodyplan.R")

bp <- make_hallucigenia(N = 100)

C        <- bp$C
k_tissue <- bp$k_tissue
```

---

## License

This project is licensed under the MIT License.  
Copyright (c) 2026 Billy Hunter

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
