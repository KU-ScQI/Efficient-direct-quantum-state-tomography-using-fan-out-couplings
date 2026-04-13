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
#### Data format
Datas can be loaded by the following code.
```bash
df = pd.read_csv("DQST_full_reconstruction_data.csv")
df['GHZ'] = df['GHZ'].apply(ast.literal_eval)
df['0state'] = df['0state'].apply(ast.literal_eval)
df['+state'] = df['+state'].apply(ast.literal_eval)
```
The dataset is stored as a pandas DataFrame with three columns:
- `GHZ`
- `0state`
- `+state`

Each column corresponds to measurement data for a specific 4-qubit state.

Each entry in these columns is a list of length 31, where each element is a dictionary representing the measurement counts obtained from a single circuit execution.

- The indices `[0]` to `[30]` correspond to the 31 circuit configurations used in the DQST protocol, in the following order:
  - `[0]`: $U_{\mathrm{ES}} = \mathrm{IIII}$ (diagonal measurement)
  - `[1]`–`[15]`: X-basis measurement circuits
  - `[16]`–`[30]`: Y-basis measurement circuits

Each dictionary has the form:

```bash
{
    '00000': count,
    '00001': count,
    ...
    '11111': count
}```

### 2. GHZ-state fidelity estimation
`GHZ_fidelity_estimation_data.csv` contains experimental data obtained from IBM Quantum (backend: IBM_aachen) for system sizes n = 4, 5, 6, 7, 8, 9, 10, 15, and 20, where n denotes the number of system qubits.

For each value of n, the GHZ fidelity is estimated using a single circuit configuration given by
$U_{\mathrm{ES}}^{\mathbf{1}} = X^{\otimes n}$.



### 3. Confusion matrices
`confusion_full_reconstruction.pkl` contains the readout confusion matrices for each qubit used in the density matrix reconstruction experiments.

`confusion_GHZ_fidest.pkl` contains the readout confusion matrices for each qubit used in the n-qubit GHZ state fidelity estimation experiments.

### 4. Demo.ipynb
`Demo.ipynb` demonstrates how the data from `DQST_full_reconstruction_data.csv` and `GHZ_fidelity_estimation_data.csv` are processed to reproduce Fig. 2, Table 1 (DQST results), and Fig. 3b and 3c.
