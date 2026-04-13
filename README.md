<img width="581" height="331" alt="image" src="https://github.com/user-attachments/assets/293c03fc-2041-4383-a657-c275bdc1d245" /># Efficient-direct-quantum-state-tomography-using-fan-out-couplings

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

Each entry in these columns is a list of length 31, where each element is a dictionary representing the measurement counts obtained from a single circuit execution. The counts have been preprocessed to correct the bitstring ordering of IBM Quantum outputs.

- The indices `[0]` to `[30]` correspond to the 31 circuit configurations used in the DQST protocol, in the following order:
  - `[0]`: $U_{\mathrm{ES}} = \mathrm{IIII}$ (diagonal measurement)
  - `[1]`–`[15]`: X-basis measurement circuits
  - `[16]`–`[30]`: Y-basis measurement circuits

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
`GHZ_fidelity_estimation_data.csv` contains experimental data obtained from IBM Quantum (backend: IBM_aachen) for system sizes n = 4, 5, 6, 7, 8, 9, 10, 15, and 20, where n denotes the number of system qubits.

For each value of n, the GHZ fidelity is estimated using a single circuit configuration given by
$U_{\mathrm{ES}}^{\mathbf{1}} = X^{\otimes n}$.
#### Data format
Datas can be loaded by the following code.
```bash
df_GHZ_fidest = pd.read_csv("GHZ_fidelity_estimation_data.csv", header=[0,1],index_col=0)
df_GHZ_fidest = df_GHZ_fidest.applymap(ast.literal_eval)
```

The GHZ fidelity estimation data are stored in a pandas DataFrame indexed by system size:
- `n=4, 5, 6, 7, 8, 9, 10, 15, 20`

Each entry (e.g., `df_GHZ_fidest['n=4']`) contains results for different noise amplification factors:
- `zne=1`, `zne=3`, `zne=5`

These correspond to repeating the $U_{\mathrm{ES}}$ gate 1, 3, and 5 times, respectively.

For each `(n, zne)` pair, the data consist of a list of length 100:
- Each element corresponds to a randomly sampled Pauli twirling instance,
- A total of 100 random Pauli sets were created.

Each element in the list is a dictionary of measurement counts.


### 3. Confusion matrices
`confusion_full_reconstruction.pkl` contains the readout confusion matrices for each qubit used in the density matrix reconstruction experiments. 

`confusion_GHZ_fidest.pkl` contains the readout confusion matrices for each qubit used in the n-qubit GHZ state fidelity estimation experiments.

#### Data format
Datas can be loaded by the following code.

```bash
confu_4q=pd.read_pickle("confusion_full_reconstruction.pkl") #for density matrix reconstruction experiment
df_confu = pd.read_pickle("confusion_GHZ_fidest.pkl") # for GHZ-state fidelity estimation experiment
```
##### 1) Data format for confusion matrix datas in density matrix reconstruction experiment
For the density matrix reconstruction experiment, the confusion matrix is stored as a full-system matrix.

`confu_4q` is a $2^5 \times 2^5$ NumPy array obtained by taking the tensor product of the single-qubit confusion matrices.

##### 2) Data format for confusion matrix datas in GHZ-state fidelity estimation experiment
The confusion matrices are stored in a dictionary-like structure indexed by system size:
- `n=4, 5, 6, 7, 8, 9, 10, 15, 20`

Each entry (e.g., `df_confu['n=4']`) is a list of 2×2 readout confusion matrices, one for each qubit used in the experiment.

For example:

```python
df_confu['n=4']
```
returns a list of length 5 (4 system qubits + 1 ancilla), where each element is a 2$\times 2$ matrix of the form 

```python
[[P(0|0), P(0|1)],
 [P(1|0), P(1|1)]]
```
Here, $P(i|j)$ dentoes the probability of measuring state i given that the true state is j.

### 4. Demo.ipynb
`Demo.ipynb` demonstrates how the data from `DQST_full_reconstruction_data.csv` and `GHZ_fidelity_estimation_data.csv` are processed to reproduce Fig. 2, Table 1 (DQST results), and Fig. 3b and 3c.
