# Stochastic Modeling of Partially Stirred Reactor (PaSR)
The _PaSR4Comb_ code is used to Turbulence-Chemistry Interaction for combustion system. 


## Numerical Implementation
#### Through-flow
```math
N_\text{replaced} = N_p \cdot \frac{\Delta t}{\tau_\text{res}}
```
#### Molecular mixing process 
+ IEM
```math
\frac{\text{d} \phi^{(k)}}{\text{d} t} = - \frac{1}{2} C_\phi \omega_\text{mix} \cdot ( \phi^{(k)} - \langle \phi \rangle )
```

```math
\phi^{(k)}_\text{new} = \langle \phi \rangle + \left( \phi^{(k)} - \langle \phi \rangle \right) \cdot \exp{\left(-\frac{1}{2} C_\phi \omega_\text{mix} \Delta t \right) }
```


+ MCM
+ MMC

> [!NOTE]
> The MMC mixing model will be updated soon.The block for MMC is still not be tested intensively.

#### Chemical reaction process

```math
\frac{\text{d}\phi^{(k)}}{\text{d}t}=\mathcal{R}(\phi^{(k)}),
```

## Other models in _PaSR4Comb_ Code
#
