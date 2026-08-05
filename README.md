# Efficient-direct-quantum-state-tomography-using-fan-out-couplings

Code for reproducing results in Efficient direct quantum state tomography using fan-out couplings.

## Installation
```bash
pip install -r requirements.txt
```
## Usage
### 1. Density matrix reconstruction
`DQST_full_reconstruction_data.csv` contains experimental data obtained from IBM Quantum (backend: IBM_aachen) for three 4-qubit states: the GHZ state, the computational basis state |0000⟩, and the |++++⟩ state.

#### Data format
Datas can be loaded by the following code.
For datas without readout mitigation:
```bash
df = pd.read_csv("DQST_full_reconstruction_data.csv")
df['GHZ'] = df['GHZ'].apply(ast.literal_eval)
df['0state'] = df['0state'].apply(ast.literal_eval)
df['+state'] = df['+state'].apply(ast.literal_eval)
```

For datas with readout mitigation applied:
```bash
df_QREM = pd.read_csv("DQST_full_reconstruction_data_QREM.csv")
df_QREM['GHZ'] = df_QREM['GHZ'].apply(ast.literal_eval)
df_QREM['0state'] = df_QREM['0state'].apply(ast.literal_eval)
df_QREM['+state'] = df_QREM['+state'].apply(ast.literal_eval)
```
The dataset is stored as a pandas DataFrame with three columns:
- `GHZ`
- `0state`
- `+state`

Each column corresponds to measurement data for a specific 4-qubit state.

Each entry in these columns is a list of length 31, where each element is a dictionary representing the measurement counts obtained from a single circuit execution. The counts have been preprocessed to correct the bitstring ordering of IBM Quantum outputs.

- The indices `[0]` to `[30]` correspond to the 31 circuit configurations used in the DQST protocol, in the following order:
  - `[0]`: $U_{\mathrm{ES}}^{\mathbf{0}} = \mathrm{IIII}$ (diagonal measurement)
  - `[1]`–`[15]`: the set $U_{\mathrm{ES}}^{\mathbf{k}}=$ {XXXX, XIII, IXII, IIXI, IIIX, XXII, XIXI, XIIX, IXXI, IXIX, IIXX, XXXI, XXIX, XIXX, IXXX} with the meter qubit X-basis measurement circuits
  - `[16]`–`[30]`: the same set with the meter qubit Y-basis measurement circuits

Each dictionary has the form:

```python
{
    '00000': count,
    '00001': count,
    ...
    '11111': count
}
```

### 2. GHZ-state fidelity estimation
`GHZ_fidelity_estimation_data_raw.csv` contains experimental data obtained from IBM Quantum (backend: IBM_aachen) for system sizes n = 4, 5, 6, 7, 8, 9, 10, 15, and 20, where n denotes the number of system qubits. 

For each value of n, the GHZ fidelity is estimated using a single circuit configuration given by
$U_{\mathrm{ES}}^{\mathbf{1}} = X^{\otimes n}$.
#### Data format
Datas can be loaded by the following code.
For datas without readout mitigation:
```bash
df_GHZ_raw = pd.read_csv("GHZ_fidelity_estimation_data_raw.csv", header=[0,1],index_col=0)
df_GHZ_raw = df_GHZ_raw.applymap(ast.literal_eval)
```

For datas with readout mitigation applied:
```bash
df_GHZ_mit = pd.read_csv("GHZ_fidelity_estimation_data_QREM.csv", header=[0,1],index_col=0)
df_GHZ_mit = df_GHZ_mit.applymap(ast.literal_eval)
```

The GHZ fidelity estimation data are stored in a pandas DataFrame indexed by system size:
- `n=4, 5, 6, 7, 8, 9, 10, 15, 20`

Each entry (e.g., `df_GHZ_raw['n=4']`) contains results for different noise amplification factors:
- `zne=1`, `zne=3`, `zne=5`

These correspond to repeating the $U_{\mathrm{ES}}^{\mathbf{1}}$ gate 1, 3, and 5 times, respectively.

For each `(n, zne)` pair, the data consist of a list of length 100:
- Each element corresponds to a randomly sampled Pauli twirling instance,
- A total of 100 random Pauli sets were created.

Each element in the list is a dictionary of measurement counts. For the GHZ-state fidelity estimation, only four outcomes are required: 00...00, 00...01, 11...10, 11...11. Therefore, for every system size n and each ZNE noise-scaling factor, we report only the counts corresponding to these four bit strings. Both the raw measurement counts and the counts after QREM are provided in this format.


### 4. Demo.ipynb
`Demo.ipynb` demonstrates how the data from `DQST_full_reconstruction_data.csv` and `GHZ_fidelity_estimation_data.csv` are processed to reproduce Fig. 2, Table 1 (DQST results), and Fig. 3b and 3c, which constitute the main results of this work.

`DQST_full_tomography.py` contains the code used to reproduce Fig. 2 and Table 1 (DQST results).

`GHZ_fidelity_estimation.py` contains the code used to reproduce Fig. 3b and 3c.

`helper.py` contains additional utility functions used throughout the codebase.


