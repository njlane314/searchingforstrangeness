# Inclusive Strangeness (post-FSI)

An event is counted if the generator final-state list (particles exiting the
nucleus, pre-detector) contains one or more hadrons with non-zero net
strangeness (K±, K⁰/K̄⁰, K_S, K_L, Λ/Σ/Ξ/Ω), excluding

1. kaons whose immediate mother is a final-state φ(1020), and
2. strange hadrons with a charm-hadron or τ ancestor when tracing ancestry
   up to the nearest final-state ancestor.

This definition is intended to be generator-agnostic and free of switches or
variants. In LArSoft/GENIE the final-state set corresponds to particles with
`StatusCode()==1` and `Process()=="primary"`. The TruthAnalysis tool contains a
compact implementation that mirrors the following pseudocode:

```python
for h in final_state_particles:
    if h is not open strangeness: continue
    if h is kaon and mother is final-state phi: veto
    walk ancestry to nearest final-state ancestor;
    if any ancestor is charm hadron or tau: veto
    return True  # surviving open strangeness
return False
```

## PDG Codes

- Kaons: 321, −321, 311, −311, 310 (K_S), 130 (K_L)
- Strange baryons: 3122, 3112, 3212, 3222, 3312, 3322, 3334
- φ(1020) (hidden strangeness): 333
- Heavy-feeddown veto parents: open-charm hadrons (4xx, 4xxx excluding 441/443/445) and τ (15)

The TruthAnalysis tool applies the two vetoes above and returns `true` if at
least one surviving open-strangeness hadron remains in the final-state set. For
downstream studies, it also fills an `InclusiveStrangenessResult` structure with
diagnostics:

- `has_inclusive_strangeness`: event passes the definition above
- `n_phi`: number of final-state φ(1020)
- `has_phi_daughter`: at least one kaon removed by the φ-daughter veto
- `has_heavy_feeddown`: open-strangeness carrier removed due to a charm/τ ancestor

These flags let analyses independently track hidden-strangeness and heavy-flavour
feed-through channels while still using a common inclusive-strangeness tag.
