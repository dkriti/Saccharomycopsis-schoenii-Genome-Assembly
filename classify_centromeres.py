def classify_chromosome(chr_name, chr_length, cen_start, cen_end):
    arm1 = cen_start
    arm2 = chr_length - cen_end

    q_arm = max(arm1, arm2)
    p_arm = min(arm1, arm2)

    # Arm ratio (r)
    if p_arm == 0:
        ratio = float('inf')
    else:
        ratio = q_arm / p_arm

    # Based on Levan et al. (1964) rules
    if 1.0 <= ratio < 1.7:
        c_type = "Metacentric"
    elif 1.7 <= ratio < 3.0:
        c_type = "Submetacentric"
    elif 3.0 <= ratio < 7.0:
        c_type = "Acrocentric"
    else:
        c_type = "Telocentric"

    print(f"{chr_name}:")
    print(f"  Short Arm (p): {p_arm:,} bp | Long Arm (q): {q_arm:,} bp")
    print(f"  Arm Ratio (r): {ratio:.2f} -> Type: {c_type}\n")

classify_chromosome("Chromosome I", 1831595, 1280001, 1290000)
classify_chromosome("Chromosome II", 3497966, 518406, 528405)
classify_chromosome("Chromosome III", 3172219, 2720440, 2730439)
classify_chromosome("Chromosome IV", 2335549, 498221, 508220)
classify_chromosome("Chromosome V", 2250200, 2222672, 2227671)
classify_chromosome("Chromosome VI", 1178080, 297472, 302471)
