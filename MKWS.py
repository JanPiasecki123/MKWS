import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

lambdas = np.linspace(0.8, 1.3, 6)
egr_rates = np.arange(0.0, 0.35, 0.025)

stability_lambda = []
stability_max_egr = []
tabela_txt = []
quench_points = []

tabela_txt.append(f"{'Lambda (L)':<12} | {'EGR [%]':<10} | {'Velocity S_L [cm/s]':<22} | {'Max Temperature [K]':<20} | {'Flame Status':<18}")
tabela_txt.append("-" * 92)

fig1, ax1 = plt.subplots(figsize=(7, 5))
fig2, ax2 = plt.subplots(figsize=(7, 5))
fig3, ax3 = plt.subplots(figsize=(7, 5))

kolory = plt.cm.plasma(np.linspace(0, 0.85, len(lambdas)))

print("Starting 1D FreeFlame simulation with physical quenching limit...")

def solve_flame(phi, egr, prev_flame=None):
 
    gas = ct.Solution('gri30.yaml')
    o2_frac = 0.21 * (1.0 - egr)
    n2_frac = 0.79 * (1.0 - egr)
    oxidizer_comp = f'O2:{o2_frac}, N2:{n2_frac}, CO2:{egr}'
    gas.set_equivalence_ratio(phi, 'CH4', oxidizer_comp)
    gas.TP = 300.0, ct.one_atm

    flame = ct.FreeFlame(gas, width=0.03)
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)

    if prev_flame is not None:
        prev_arr = prev_flame.to_array()
        flame.set_initial_guess(data=prev_arr)

    flame.solve(loglevel=0, auto=True)
    S_L = flame.velocity[0] * 100.0
    T_max = float(max(flame.T))
    return flame, S_L, T_max


def bisect_quench_limit(phi, egr_lo, egr_hi, prev_flame, tol=0.002):
  
    print(f"     [~] Bisecting quench limit between {egr_lo*100:.2f}% and {egr_hi*100:.2f}%...")
    best_stable_egr = egr_lo
    best_flame = prev_flame

    while (egr_hi - egr_lo) > tol:
        egr_mid = (egr_lo + egr_hi) / 2.0
        try:
            fl, S_L, _ = solve_flame(phi, egr_mid, prev_flame=best_flame)
            if S_L < 3.5:
                egr_hi = egr_mid
            else:
                egr_lo = egr_mid
                best_stable_egr = egr_mid
                best_flame = fl
        except Exception:
            egr_hi = egr_mid

    return best_stable_egr, best_flame


for idx, l in enumerate(lambdas):
    phi = 1.0 / l
    speeds_for_this_lambda = []
    egr_for_this_lambda = []
    temp_for_this_lambda = []
    max_stable_egr = 0.0
    prev_flame = None   # warm-start carrier

    print(f"---> Analyzing mixture for Lambda = {l:.2f}")

    quench_egr_coarse = None  # first EGR step that triggered quench/failure

    for egr in egr_rates:
        try:
            flame, S_L, T_max = solve_flame(phi, egr, prev_flame=prev_flame)

            if S_L < 3.5:
                print(f"     [!] Coarse quench detected at EGR = {egr*100:.1f}% (S_L = {S_L:.2f} cm/s)")
                quench_egr_coarse = egr
                break

            # Stable point — record and carry flame forward as warm start
            speeds_for_this_lambda.append(S_L)
            temp_for_this_lambda.append(T_max)
            egr_for_this_lambda.append(egr * 100)
            max_stable_egr = egr * 100
            prev_flame = flame

            tabela_txt.append(f"{l:<12.2f} | {egr*100:<10.1f} | {S_L:<22.2f} | {T_max:<20.1f} | {'Stable':<18}")
            print(f"     EGR = {egr*100:>4.1f}% | S_L = {S_L:>5.2f} cm/s | T_max = {T_max:>6.1f} K")

        except Exception as exc:
            print(f"     [!] Solver failed at EGR = {egr*100:.1f}% — {type(exc).__name__}: {exc}")
            quench_egr_coarse = egr
            break

    if quench_egr_coarse is not None and prev_flame is not None:
        egr_lo = max_stable_egr / 100.0
        egr_hi = quench_egr_coarse
        precise_egr, precise_flame = bisect_quench_limit(phi, egr_lo, egr_hi, prev_flame)
        max_stable_egr = precise_egr * 100.0

        try:
            _, S_L_q, T_max_q = solve_flame(phi, precise_egr, prev_flame=precise_flame)
        except Exception:
            S_L_q, T_max_q = 0.0, 300.0

        tabela_txt.append(f"{l:<12.2f} | {precise_egr*100:<10.3f}    | {S_L_q:<22.2f} | {T_max_q:<20.1f} | {'Quenched (Limit)':<18}")
        quench_points.append((l, precise_egr * 100.0, S_L_q, T_max_q))
        print(f"     [✓] Precise quench limit: EGR = {max_stable_egr:.3f}%")

    ax1.plot(egr_for_this_lambda, speeds_for_this_lambda, marker='o', markersize=4, linestyle='-', linewidth=1.5, color=kolory[idx], label=f'L={l:.2f}')
    ax3.plot(egr_for_this_lambda, temp_for_this_lambda, marker='x', markersize=4, linestyle='-', linewidth=1.5, color=kolory[idx], label=f'L={l:.2f}')

    stability_lambda.append(l)
    stability_max_egr.append(max_stable_egr)

quench_txt = []
quench_txt.append(f"{'Lambda (L)':<12} | {'Precise EGR Limit [%]':<22} | {'S_L at Limit [cm/s]':<22} | {'T_max at Limit [K]':<20}")
quench_txt.append('-' * 84)
for lam, egr_pct, sl, tmax in quench_points:
    quench_txt.append(f"{lam:<12.2f} | {egr_pct:<22.3f} | {sl:<22.2f} | {tmax:<20.1f}")

with open('Tabela_Granic_Gaszenia.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(quench_txt))

print("Saved Tabela_Granic_Gaszenia.txt")
print("\nGenerating plots...")

ax1.set_title('Laminar Burning Velocity vs EGR', fontsize=12)
ax1.set_xlabel('EGR Fraction [%]')
ax1.set_ylabel('Burning Velocity S_L [cm/s]')
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend(title='Lambda (L)', loc='upper right', fontsize=9)
fig1.tight_layout()
fig1.savefig('Wykres_1_Predkosc.png', dpi=300, bbox_inches='tight')

ax2.plot(stability_lambda, stability_max_egr, color='#C0392B', linewidth=3, marker='s', markersize=6)
ax2.fill_between(stability_lambda, 0, stability_max_egr, color='#2ECC71', alpha=0.15, label='Flammability Region')
EGR_SCAN_MAX_PCT = (egr_rates[-1] + egr_rates[1] - egr_rates[0]) * 100  # upper bound of scan range
ax2.fill_between(stability_lambda, stability_max_egr, EGR_SCAN_MAX_PCT, color='#E74C3C', alpha=0.1, label='Quenching Region')
ax2.set_title('Stability Limits Map', fontsize=12)
ax2.set_xlabel('Air Excess Ratio (Lambda)')
ax2.set_ylabel('Max Allowable EGR [%]')
ax2.grid(True, linestyle='--', alpha=0.5)
ax2.legend(loc='lower center', fontsize=9)
fig2.tight_layout()
fig2.savefig('Wykres_2_Mapa.png', dpi=300, bbox_inches='tight')

ax3.set_title('Max Flame Temperature vs EGR', fontsize=12)
ax3.set_xlabel('EGR Fraction [%]')
ax3.set_ylabel('Temperature [K]')
ax3.grid(True, linestyle='--', alpha=0.5)
ax3.legend(title='Lambda (L)', loc='upper right', fontsize=9)
fig3.tight_layout()
fig3.savefig('Wykres_3_Temperatura.png', dpi=300, bbox_inches='tight')

with open('Tabela_Wynikow.txt', 'w', encoding='utf-8') as file:
    file.write('\n'.join(tabela_txt))

print("=== SUCCESS ===")
print("Saved 3 separate plot files, Tabela_Wynikow.txt, and Tabela_Granic_Gaszenia.txt")
plt.show()
