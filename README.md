# Stochastic Modeling of Partially Stirred Reactor (_PaSR4Comb_)
The _PaSR4Comb_ code is used to Turbulence-Chemistry Interaction for combustion system. 

## Stochastic modeling

```math
    \frac{\partial \widetilde{f}_{\phi}(\psi,t)}{\partial t}  =  
    +  \frac{1}{\tau_\text{res}} \bigl\{ \widetilde{f}_{\phi,\text{inlet}}(\psi,t) - \widetilde{f}_{\phi}(\psi,t) \bigr\}
    - \sum_{\alpha=1}^{n_s+2} \frac{\partial}{\partial \psi_\alpha} \bigl\{ R_\alpha (\Vec{\psi}) \cdot \widetilde{f}_{\phi}(\psi,t) \bigr\} - \sum_{\alpha=1, \beta=1}^{n_s+2} \frac{\partial^2}{\partial \psi_\alpha \psi_\beta} \bigl\{ \langle \varepsilon_{\alpha \beta} | \phi = \psi \rangle \widetilde{f}_{\phi}(\psi,t) \bigr\}
```

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

#### Chemical reaction process (detailed chemistry)

```math
\frac{\text{d}\phi^{(k)}}{\text{d}t}=\mathcal{R}(\phi^{(k)}),
```

## Chemistry models
#### Quasi-steady state assumption (_QSSA_)

The simulation based on the GQL reduced chemistry is to solve the DAE system [2]:

```math
\textbf{Q}_s \frac{\text{d}\phi}{\text{d}t}=R(\phi),
```
where $\textbf{M}_s$ is the mass matrix defined as:

```math
\textbf{Q}_s = \begin{pmatrix}
                  1 ~ 0 \\
                  0 ~ 1
                  \end{pmatrix},
```

#### Virtual chemistry (_VC_)
> [!WARNING]
> in progress
