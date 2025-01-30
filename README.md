# Stochastic Modeling of Partially Stirred Reactor (PaSR)
The _PaSR4Comb_ code is used to Turbulence-Chemistry Interaction for combustion system. 


## Numerical Implementation
### Through-flow
```math
N_\text{replaced} = N_p \cdot \frac{\Delta t}{\tau_\text{res}}
```
### Molecular mixing process 
+ IEM
+ MCM
+ MMC

> [!NOTE]
> The MMC mixing model will be updated soon.The block for MMC is still not be tested intensively.

### Chemical reaction process

```math
\frac{\text{d}\Vec{\phi}^{(k)}}{\text{d}t}=\mathcal{R}(\vec{\phi}^{(k)}),
```
