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


def Bootstrapped_data(data1, data2, sampling_number):
    zne_1 = []
    zne_3 = []
    zne_5 = []

    zne_1_mit = []
    zne_3_mit = []
    zne_5_mit = []

    for i in range(1):
        zne_1.extend(data1['zne=1'])
        zne_3.extend(data1['zne=3'])
        zne_5.extend(data1['zne=5'])

        zne_1_mit.extend(data2['zne=1'])
        zne_3_mit.extend(data2['zne=3'])
        zne_5_mit.extend(data2['zne=5'])

    Bootdata = []
    Bootdata_mit = []

    for _ in range(sampling_number):

        idx1 = random.choices(range(len(zne_1)), k=100)
        idx3 = random.choices(range(len(zne_3)), k=100)
        idx5 = random.choices(range(len(zne_5)), k=100)

        zne_1_boot = [zne_1[i] for i in idx1]
        zne_3_boot = [zne_3[i] for i in idx3]
        zne_5_boot = [zne_5[i] for i in idx5]

        zne_1_boot_mit = [zne_1_mit[i] for i in idx1]
        zne_3_boot_mit = [zne_3_mit[i] for i in idx3]
        zne_5_boot_mit = [zne_5_mit[i] for i in idx5]

        Bootdata.append([zne_1_boot, zne_3_boot, zne_5_boot])
        Bootdata_mit.append([zne_1_boot_mit, zne_3_boot_mit, zne_5_boot_mit])

    return Bootdata, Bootdata_mit


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

def Elaborate_GHZ_fidest(data1, num_qubit):
    
    # 1) No QREM (with / without ZNE)
    Boot_list_with_zne = []
    Boot_list_no_zne = []

    for d in data1[0]:
        zne_vals = elaborate_zne_from_counts(d)
        Boot_list_no_zne.append(zne_vals[0][0])
        Boot_list_with_zne.append(zne_vals[-1][0])

    mean_nozne_noQREM   = np.mean(Boot_list_no_zne)
    std_nozne_noQREM    = np.std(Boot_list_no_zne)
    mean_withzne_noQREM = np.mean(Boot_list_with_zne)
    std_withzne_noQREM  = np.std(Boot_list_with_zne)

    # 2) With QREM (with / without ZNE)
    QREM_Boot_list_with_zne = []
    QREM_Boot_list_no_zne = []
    for d in data1[1]:
        QREM_zne_vals = elaborate_zne_from_counts(d)
        QREM_Boot_list_no_zne.append(QREM_zne_vals[0][0])
        QREM_Boot_list_with_zne.append(QREM_zne_vals[-1][0])

    mean_nozne_withQREM   = np.mean(QREM_Boot_list_no_zne)
    std_nozne_withQREM    = np.std(QREM_Boot_list_no_zne)
    mean_withzne_withQREM = np.mean(QREM_Boot_list_with_zne)
    std_withzne_withQREM  = np.std(QREM_Boot_list_with_zne)

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


def visualize_GHZ_fidest_result(boot_data,n_vals):
    all_results = [Elaborate_GHZ_fidest(boot_data[f'n={n}'], n) for n in n_vals]

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
    plt.ylabel("GHZ Fidelity", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(fontsize=20)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

