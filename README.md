# Efficient-direct-quantum-state-tomography-using-fan-out-couplings

Code for reproducing results in Efficient direct quantum state tomography using fan-out couplings.

## Installation
```bash
pip install -r requirements.txt
```
## Usage
### 1. Density matrix reconstruction
`DQST_full_reconstruction_data.csv` contains experimental data obtained from IBM Quantum (backend: IBM_aachen) for three 4-qubit states: the GHZ state, the computational basis state |0000⟩, and the |++++⟩ state.

Each state is reconstructed using 31 measurement circuits. The circuit order is given by:
- $U_{\mathrm{ES}} = \mathrm{IIII}$ (diagonal measurement),
- the set {XXXX, XIII, IXII, IIXI, IIIX, XXII, XIXI, XIIX, IXXI, IXIX, IIXX, XXXI, XXIX, XIXX, IXXX} with X-basis measurements,
- the same set with Y-basis measurements.

### 2. GHZ-state fidelity estimation
`GHZ_fidelity_estimation_data.csv` contains experimental data obtained from IBM Quantum (backend: IBM_aachen) for system sizes n = 4, 5, 6, 7, 8, 9, 10, 15, and 20, where n denotes the number of system qubits.

For each value of n, the GHZ fidelity is estimated using a single circuit configuration given by
$U_{\mathrm{ES}}^{\mathbf{1}} = X^{\otimes n}$.

### 3. Confusion matrices
`confusion_full_reconstruction.pkl` contains the readout confusion matrices for each qubit used in the density matrix reconstruction experiments.

`confusion_GHZ_fidest.pkl` contains the readout confusion matrices for each qubit used in the n-qubit GHZ state fidelity estimation experiments.
