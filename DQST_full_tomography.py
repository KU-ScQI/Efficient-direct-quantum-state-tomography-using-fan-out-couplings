from helper import *

def obtain_physical_dens(rho):

    rho = (rho + rho.conj().T) / 2
    n = rho.shape[0]
    rho_physical = cp.Variable((n, n), complex=True)
    objective = cp.Minimize(cp.norm(rho - rho_physical, "fro"))
    constraints = [
        rho_physical >> 0, 
        cp.trace(rho_physical) == 1, 
        rho_physical == rho_physical.H 
    ]
    problem = cp.Problem(objective, constraints)
    problem.solve()
    
    return rho_physical.value

def state_space_projection(density_matrix):

    rho_unphysical = density_matrix.full()
    rho_physical = obtain_physical_dens(rho_unphysical)
    physical_density_matrix=qutip.Qobj(rho_physical)
        
    return physical_density_matrix


def process_measurement(data_dict, total_shot, order, op='diag'):
    grouped = defaultdict(lambda: [0, 0])

    for key, val in data_dict.items():
        prefix, bit = key[:4], int(key[4])
        grouped[prefix][bit] = val

    if op == 'diag':
        result = {k: (v[0] + v[1]) / total_shot for k, v in grouped.items()}
    else:  # op == 'off_diag'
        result = {k: (v[1] - v[0]) / total_shot for k, v in grouped.items()}

    return [result[k] for k in order]

def DQST_full_tomography(count_list,total_shot):
    diag_data=count_list[0] #first datas are count data for diagonal elements

    off_diag_real_data=[]
    off_diag_imag_data=[]
    for i in range(1,16):
        off_diag_real_data.append(count_list[i]) #count_list[1]~[15] are count datas acquired for off diagonal real (meter X measurement)

    for i in range(16,31):
        off_diag_imag_data.append(count_list[i]) #count_list[16]~[30] are count datas acquired for off diagonal imaginary (meter Y measurement)

    order = ['0000', '0001', '0010', '0011', '0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111']

    diag_result = []
    diag_result.append(process_measurement(diag_data, total_shot, order, op='diag'))

    real_result = []
    for data in off_diag_real_data:
        real_result.append(process_measurement(data, total_shot, order, op='off_diag'))

    imag_result = []
    for data in off_diag_imag_data:
        imag_result.append(process_measurement(data, total_shot, order, op='off_diag'))

    circuit_index=['XXXX', 'XIII', 'IXII','IIXI','IIIX','XXII','XIXI','XIIX','IXXI','IXIX','IIXX','XXXI','XXIX','XIXX','IXXX'] #order depends on the order of circuit we executed
    
    density_matrix_element = ['00000000', '00010001', '00100010', '00110011', '01000100', '01010101', '01100110', '01110111',
                '10001000', '10011001', '10101010', '10111011', '11001100', '11011101', '11101110', '11111111']

    same_circuit=[]
    for m in circuit_index:
        for st in density_matrix_element:
            d=kb(modify_last_4_bits(st,m))
            same_circuit.append(d)

    dens = np.array(same_circuit).reshape(15, 16)
    result_off_diagonal = dens * (np.array(real_result) + 1j * np.array(imag_result))

    diagonal_term = np.array([kb(s) for s in density_matrix_element], dtype=object)
    result_diagonal = np.array(diag_result)[0] * diagonal_term

    density_matrix = np.sum(result_off_diagonal) + np.sum(result_diagonal)

    return density_matrix, state_space_projection(density_matrix)


def QREM_mitigated_counts(count_list,num_qubit,confusion_matrix):
    count_data_mit=[]
    for i in range (len(count_list)):
        count_data_mit.append(renormalize_counts(Counts_Readout_Dic(count_list[i],int(num_qubit+1),confusion_matrix)))
    return count_data_mit


def visualize_DQST_results(dens,target):
    density = dens.full()
    data_real = np.real(density)
    data_imag = np.imag(density)

    n = data_real.shape[0]
    _x = np.arange(n)
    _y = np.arange(n)
    _xx, _yy = np.meshgrid(_x, _y)

    x = _xx.ravel() - 0.7
    y = _yy.ravel() - 0.7
    z = np.zeros_like(x)

    dx = dy = 0.7

    labels = ['0000','0001','0010','0011','0100','0101','0110','0111',
              '1000','1001','1010','1011','1100','1101','1110','1111']
    important_indices = [0, 15]

    fig = plt.figure(figsize=(8, 2), dpi=300)
    ax1 = fig.add_subplot(121, projection='3d')
    dz_real = data_real.ravel()

    norm_r = Normalize(vmin=0, vmax=np.max(dz_real))
    cmap_r = cm.get_cmap('Blues')
    colors_r = cmap_r(norm_r(dz_real))

    ax1.bar3d(x, y, z, dx, dy, dz_real, color=colors_r, edgecolor='black', linewidth=0.1)

    ax1.set_xticks(important_indices)
    ax1.set_xticklabels([labels[i] for i in important_indices], fontsize=7)
    ax1.set_yticks(important_indices)
    ax1.set_yticklabels([labels[i] for i in important_indices], fontsize=7)
    ax1.set_zlim(0, 1)
    ax1.view_init(elev=30, azim=35)

    mappable_r = ScalarMappable(cmap=cmap_r, norm=norm_r)
    mappable_r.set_array(dz_real)
    cbar_r = plt.colorbar(mappable_r, ax=ax1, orientation='vertical', pad=0.1, shrink=0.6)
    cbar_r.set_label(r'$\mathrm{Re}(\rho)$', fontsize=8)

    ax2 = fig.add_subplot(122, projection='3d')
    dz_imag = data_imag.ravel()

    norm_i = Normalize(vmin=np.min(dz_imag), vmax=np.max(dz_imag))
    cmap_i = cm.get_cmap('Reds')
    colors_i = cmap_i(norm_i(dz_imag))

    ax2.bar3d(x, y, z, dx, dy, dz_imag, color=colors_i, edgecolor='black', linewidth=0.1)

    ax2.set_xticks(important_indices)
    ax2.set_xticklabels([labels[i] for i in important_indices], fontsize=7)
    ax2.set_yticks(important_indices)
    ax2.set_yticklabels([labels[i] for i in important_indices], fontsize=7)
    ax2.set_zlim(0, 1)
    ax2.view_init(elev=30, azim=35)

    mappable_i = ScalarMappable(cmap=cmap_i, norm=norm_i)
    mappable_i.set_array(dz_imag)
    cbar_i = plt.colorbar(mappable_i, ax=ax2, orientation='vertical', pad=0.1, shrink=0.6)
    cbar_i.set_label(r'$\mathrm{Im}(\rho)$', fontsize=8)

    plt.tight_layout()
    plt.show()


def Monte_carlo_list(data):
    n_qubits = 5
    all_outcomes = [''.join(bits) for bits in product('01', repeat=n_qubits)]
    count_data=data
    monte_DQST = []
    for i in range(len(count_data)):
        original_counts = count_data[i]
        total_shots = sum(original_counts.values())

        outcomes = list(original_counts.keys())
        probs = [original_counts[o] / total_shots for o in outcomes]

        all_resampled_counts = []
        for _ in range(500):
            resampled = np.random.choice(outcomes, size=10000, p=probs)
            resampled_counts = dict(Counter(resampled))
            full_counts = {outcome: resampled_counts.get(outcome, 0) for outcome in all_outcomes}
            all_resampled_counts.append(full_counts)
        monte_DQST.append(all_resampled_counts)
        
    return monte_DQST

def Monte_carlo_list_mit(monte_DQST,confu):
    monte_DQST_mitigated=[]
    for i in range (len(monte_DQST)):
        iter_circ=monte_DQST[i]
        mit_per_circ=[]
        for j in range (500):
            mit_per_circ.append(renormalize_counts(Counts_Readout_Dic(iter_circ[j],5,confu),10000))
        monte_DQST_mitigated.append(mit_per_circ)

    return monte_DQST_mitigated

def statistical_analysis(data,confu,target,shot):

    monte_nomit=Monte_carlo_list(data)
    monte_mit=Monte_carlo_list_mit(monte_nomit,confu)

    fid_list_nomit=[]
    for j in range (500):
        final_p=[]
        for i in range(len(monte_nomit)):
            final_p.append(monte_nomit[i][j])
        fid_list_nomit.append(qutip.fidelity(DQST_full_tomography(final_p,shot)[1],target))

    fid_list_mit=[]
    for j in range (500):
        final_p=[]
        for i in range(len(monte_mit)):
            final_p.append(monte_mit[i][j])
        fid_list_mit.append(qutip.fidelity(DQST_full_tomography(final_p,shot)[1],target))
        

    mean_fid_nomit = np.mean(np.array(fid_list_nomit))
    std_fid_nomit = np.std(np.array(fid_list_nomit))
    mean_fid_mit = np.mean(np.array(fid_list_mit))
    std_fid_mit = np.std(np.array(fid_list_mit))

    m1, s1 = mean_fid_nomit*100,   std_fid_nomit*100
    m2, s2 = mean_fid_mit*100, std_fid_mit*100


    print(f'1) Without QREM:      {m1:.2f} ± {s1:.2f} %')
    print(f'2) With QREM:    {m2:.2f} ± {s2:.2f} %')

    return mean_fid_nomit, std_fid_nomit, mean_fid_mit, std_fid_nomit