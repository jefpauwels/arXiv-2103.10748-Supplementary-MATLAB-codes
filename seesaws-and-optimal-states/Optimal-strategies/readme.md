.mat files with explicit strategies have a data format consistent with their respective see-saw optimization codes.

**Information-constrained communication**

- `Gallego-1bit.mat`
     
  `rho(:,:,x)` is the density matrix corresponding to $x$th preparation of Alice.
  `B(:,:,y,b)` is the measurement operator corresponding to outcome $b$ of Bob's $y$th measurement.
  
**Entanglement-assisted classical communication**
  
  - `Gallego-EntC2.mat`
  - `Gallego-EntC3.mat`
  - `Gallego-EntC4.mat`
  - `RAC-EntC3.mat`


 `rho(:,:,x,c)` is the density matrix corresponding to conditional marginal state of Bob, when Alice got outcome $c$ for measurement $x$.
 `E(:,:,y,c)` is the measurement observable of Bob when he receives the message $c$ from Alice and input $y$.


- `tritRAC-EntC3.mat`

 `rho(:,:,x,c)` is the density matrix corresponding to conditional marginal state of Bob, when Alice got outcome $c$ for measurement $x$.
  `B(:,:,y,c,b)` is the measurement operator corresponding to outcome $b$ of the measurement performed by Bob when he receives the classical message $c$ from Alice, and the classical input $y$.

**Quantum communication**

- `Gallego-Q2.mat`
- `Gallego-Q3.mat`
- `Gallego-Q4-EntQ2.mat` (note that because this strategy is real, it is equivalent to a 'generalised dense coding' protocol)

 `rho{x}` is the density matrix corresponding to $x$th preparation of Alice.
 `E{y}`is the measurement observable corresponding to Bob's $y$th measurement
