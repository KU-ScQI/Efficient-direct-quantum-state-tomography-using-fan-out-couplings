from helper import *

def GHZ_fidelity(data_list):
    GHZ_fidelity = []
    for counts in data_list:
        total_count = 100000
        num_qubits = len(next(iter(counts))) 

        leading_zero_count = sum(v for k, v in counts.items() if k[:-1] == '0' * (num_qubits - 1))

        all_one_except_last_0 = sum(v for k, v in counts.items() if k[:-1] == '1' * (num_qubits - 1) and k[-1] == '0')
        all_one_except_last_1 = sum(v for k, v in counts.items() if k[:-1] == '1' * (num_qubits - 1) and k[-1] == '1')
        leading_one_count = (all_one_except_last_1 - all_one_except_last_0)

        prob = (leading_zero_count + leading_one_count) / total_count if total_count > 0 else 0
        GHZ_fidelity.append(prob)

    return GHZ_fidelity

def indexing(confu):
    confus = np.linalg.inv(confu)
    a = confus[0][0]
    b = confus[0][1]
    c = confus[1][0]
    d = confus[1][1]

    A_0 = {'0': a, '1': b}
    A_1 = {'0': c, '1': d}
    return A_0, A_1

def build_correction_matrices(confu_list):

    return [indexing(confu) for confu in confu_list]

def correct_counts_selected(raw_counts, correction_matrices, targets):
    n = len(correction_matrices)
    corrected_counts = {}

    for corrected_bits in targets:
        total = 0.0

        for measured_key, count in raw_counts.items():
            prob = 1.0
            for i in range(n):
                A_0, A_1 = correction_matrices[i]
                meas_bit = measured_key[i]
                corr_bit = corrected_bits[i]
                prob *= A_0[corr_bit] if meas_bit == '0' else A_1[corr_bit]
            total += prob * count

        corrected_counts[corrected_bits] = total

    return corrected_counts

def RO_mitigation(data,confu_list,total_counts):
    correction_matrices = build_correction_matrices(confu_list)
    n = len(correction_matrices)
    target_keys = [
        '0' * n,
        '0' * (n - 1) + '1',
        '1' * (n - 1) + '0',
        '1' * n]

    total_corrected_counts = defaultdict(float)

    for raw_counts in data[0]:
        corrected_counts = correct_counts_selected(raw_counts, correction_matrices, target_keys)
        for key in target_keys:
            total_corrected_counts[key] += corrected_counts.get(key, 0.0)

    result_counts = abs(total_corrected_counts['0' * n] + total_corrected_counts['0' * (n - 1) + '1']) \
            + abs(total_corrected_counts['1' * (n - 1) + '0'] - total_corrected_counts['1' * n])

    result = result_counts / total_counts
    return result

def Bootstrapped_data(data, sampling_number):
    zne_1 = []
    zne_3 = []
    zne_5 = []

    for i in range(1):
        zne_1.extend(data['zne=1']) 
        zne_3.extend(data['zne=3'])
        zne_5.extend(data['zne=5'])

    Bootdata = []
    for _ in range(sampling_number):
        zne_1_boot = random.choices(zne_1, k=100)
        zne_3_boot = random.choices(zne_3, k=100)
        zne_5_boot = random.choices(zne_5, k=100)
        Bootdata.append([zne_1_boot, zne_3_boot, zne_5_boot])

    return Bootdata  # shape: [sampling_number][3][100]


def sum_counts(job_results):
    aggregated = Counter()
    for result in job_results:
        aggregated.update(result)

    return dict(aggregated)


def elaborate_zne_from_counts(count_list):

    sums = [sum_counts(c) for c in count_list]
    zne_vals = [GHZ_fidelity([s]) for s in sums]   # [zne_1, zne_3, zne_5]

    x = np.array([1, 3, 5], dtype=float)
    y = np.array(zne_vals, dtype=float)

    slope, intercept = np.polyfit(x, y, 1)

    return zne_vals[0], zne_vals[1], zne_vals[2], intercept

def Elaborate_GHZ_fidest(data, confu_lists, num_qubit):
    
    # 1) No QREM (with / without ZNE)
    Boot_list_with_zne = []
    Boot_list_no_zne = []

    for d in data:
        zne_vals = elaborate_zne_from_counts(d)
        Boot_list_no_zne.append(zne_vals[0][0])
        Boot_list_with_zne.append(zne_vals[-1][0])

    mean_nozne_noQREM   = np.mean(Boot_list_no_zne)
    std_nozne_noQREM    = np.std(Boot_list_no_zne)
    mean_withzne_noQREM = np.mean(Boot_list_with_zne)
    std_withzne_noQREM  = np.std(Boot_list_with_zne)

    # 2) With QREM (with / without ZNE)
    q_values = {1: [], 3: [], 5: []}
    confu = confu_lists[f'n={num_qubit}']

    for d in data:
        q_values[1].append(RO_mitigation([d[0]], confu, 100000))
        q_values[3].append(RO_mitigation([d[1]], confu, 100000))
        q_values[5].append(RO_mitigation([d[2]], confu, 100000))

    x = np.array([1, 3, 5])
    y = np.array([q_values[1], q_values[3], q_values[5]])

    slope, intercept = np.polyfit(x, y, 1)
    intercepts = [intercept]

    mean_withzne_withQREM = np.mean(intercepts)
    std_withzne_withQREM  = np.std(intercepts[0])

    mean_nozne_withQREM = np.mean(q_values[1])
    std_nozne_withQREM  = np.std(q_values[1])

    m1, s1 = mean_nozne_noQREM*100,   std_nozne_noQREM*100
    m2, s2 = mean_nozne_withQREM*100, std_nozne_withQREM*100
    m3, s3 = mean_withzne_noQREM*100, std_withzne_noQREM*100
    m4, s4 = mean_withzne_withQREM*100, std_withzne_withQREM*100

    print(f'---------------- Results for GHZ fidelity estimation for n={num_qubit} ---------------- ')

    print(f'1) no ZNE, no QREM:      {m1:.1f} ± {s1:.1f} %')
    print(f'2) no ZNE, with QREM:    {m2:.1f} ± {s2:.1f} %')
    print(f'3) with ZNE, no QREM:    {m3:.1f} ± {s3:.1f} %')
    print(f'4) with ZNE, with QREM:  {m4:.1f} ± {s4:.1f} %')


    return (
        mean_nozne_noQREM,  std_nozne_noQREM,
        mean_nozne_withQREM, std_nozne_withQREM,
        mean_withzne_noQREM, std_withzne_noQREM,
        mean_withzne_withQREM, std_withzne_withQREM
    )


def visualize_GHZ_fidest_result(boot_data,confu,n_vals):
    all_results = [Elaborate_GHZ_fidest(boot_data[f'n={n}'],confu, n) for n in n_vals]

    col1, col2, col3, col4, col5, col6, col7, col8 = zip(*all_results)
    y1, err1 = col1, col2
    y2, err2 = col3, col4
    y3, err3 = col5, col6
    y4, err4 = col7, col8


    plt.figure(figsize=(16, 10))

    # 1) W/O(ZNE), W/O(QREM) - black
    plt.errorbar(n_vals, y1, yerr=err1, marker='o', linestyle='none',
                color='black',markersize=18, label="None")

    # 2) W/O(ZNE), W(QREM) - gray
    plt.errorbar(n_vals, y2, yerr=err2, marker='^', linestyle='none',
                color='gray',markersize=18, label="QREM")

    # 3) W(ZNE), W/O(QREM) - red transparent
    plt.errorbar(n_vals, y3, yerr=err3, marker='d', linestyle='none',
                color='red', alpha=0.4,markersize=18, label="ZNE")

    # 4) W(ZNE), W(QREM) - red solid
    plt.errorbar(n_vals, y4, yerr=err4, marker='s', linestyle='none',
                color='red',markersize=18, label="QREM, ZNE")


    plt.xlabel("n", fontsize=20)
    plt.ylabel("GHZ Fidelity (%)", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(fontsize=20)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
