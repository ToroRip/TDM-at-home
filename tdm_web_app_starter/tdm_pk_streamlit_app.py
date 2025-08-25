
import streamlit as st
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

st.set_page_config(page_title="TDM PK Calculator", layout="wide")

st.title("🧪 TDM At home")

with st.expander("📘 Instructions", expanded=False):
    st.markdown(
        """
        **What this does**
        - Converts your Excel/VBA macro into a web UI built with Python (Streamlit).
        - Estimates PK parameters (Ke, Vd, Cl, Half-life), suggests safe infusion time, suggests dosing interval,
          and calculates a new dose to target a desired trough (default 15 mg/L).
        - Plots steady-state concentration–time curves for the current vs. adjusted regimen.

        **Notes**
        - Units assumed: dose in **mg**, weight in **kg**, times in **hours**, concentrations in **mg/L**.
        - Model options: **Bolus** or **Infusion** (intermittent infusion).
        - If **Cpeak** is blank, we infer it using Vd = 0.7 * weight (same assumption as your VBA).
        """
    )

# --- Sidebar Inputs ---
st.sidebar.header("Inputs")

# Required inputs
realCt = st.sidebar.number_input("Ctrough (mg/L) [required]", min_value=0.0, value=10.0, step=0.1, format="%.2f")
dose = st.sidebar.number_input("Dose per administration (mg)", min_value=1.0, value=1000.0, step=10.0, format="%.0f")
tau = st.sidebar.number_input("Dosing interval τ (hr)", min_value=0.1, value=12.0, step=0.5, format="%.1f")
tinf = st.sidebar.number_input("Infusion time (hr) [use 0 for Bolus]", min_value=0.0, value=1.0, step=0.25, format="%.2f")
modelUsed = st.sidebar.selectbox("Model", options=["Bolus", "Infusion"], index=1)
weight = st.sidebar.number_input("Weight (kg)", min_value=1.0, value=70.0, step=0.5, format="%.1f")

# Optional inputs
realCp_opt = st.sidebar.text_input("Cpeak (mg/L) [optional; blank to infer]")
Ct_target = st.sidebar.number_input("Target trough (mg/L)", min_value=0.1, value=15.0, step=0.5, format="%.1f")
tStart = st.sidebar.number_input("Plot start time tStart (hr)", min_value=0.0, value=0.0, step=0.5, format="%.1f")
custom_suggested_tau = st.sidebar.text_input("Override Suggested Interval (hr) [optional]")

# Validate Cpeak input
realCp = None
if realCp_opt.strip() != "":
    try:
        realCp = float(realCp_opt)
        if realCp <= 0:
            st.sidebar.error("Cpeak must be positive if provided.")
            st.stop()
    except Exception:
        st.sidebar.error("Invalid Cpeak: must be numeric or blank.")
        st.stop()

# --- Helper functions ---
def ceiling_precise(x, step):
    # Mimics Excel CEILING.PRECISE for positive x
    if step == 0:
        return x
    m = math.ceil(x / step)
    return m * step

def estimate_params(realCt, realCp, dose, tau, tinf, modelUsed, weight):
    # Returns dict with ke, vd, Cl, realCp_used, halfLife; raises ValueError on invalid
    if realCt <= 0:
        raise ValueError("Ctrough must be > 0.")
    if dose <= 0 or tau <= 0 or weight <= 0:
        raise ValueError("Dose, tau, and weight must be > 0.")
    # Infer Cp if missing
    if realCp is None:
        vd = 0.7 * weight
        realCp_used = (dose / vd) + realCt
        # VBA uses ln(Cp/Ct)/(tau - 0.5)
        denom = (tau - 0.5)
        if denom <= 0:
            raise ValueError("tau must be > 0.5 hr when Cpeak is inferred.")
        ke = math.log(realCp_used / realCt) / denom
        # In VBA, Bolus and Infusion branches both compute Cl = ke * vd when Cp is inferred
        Cl = ke * vd
    else:
        realCp_used = realCp
        # VBA uses a fixed 2.5 hr between peak and trough sample when Cp provided
        ke = math.log(realCp_used / realCt) / 2.5
        if modelUsed == "Bolus":
            # vd from trough at time (tau - 0.5) hr post-dose at steady state
            # Rearranged from Ct = (dose/vd)*exp(-ke*(tau-0.5)) / (1 - exp(-ke*tau))
            num = dose * math.exp(-ke * (tau - 0.5))
            den = realCt * (1 - math.exp(-ke * tau))
            vd = num / den
            Cl = ke * vd
        else:
            # Intermittent infusion: use VBA rearrangement for Cl first, then vd = Cl/ke
            # Ct = [dose * exp(-ke * (tau - tinf - 0.5)) * (1 - exp(-ke * tinf))] / [Cl * tinf * (1 - exp(-ke * tau))]
            num = dose * math.exp(-ke * (tau - tinf - 0.5)) * (1 - math.exp(-ke * tinf))
            den = realCt * tinf * (1 - math.exp(-ke * tau))
            Cl = num / den
            vd = Cl / ke

    if ke <= 0 or vd <= 0 or Cl <= 0 or not np.isfinite(ke*vd*Cl):
        raise ValueError("Invalid PK parameters calculated. Please check inputs.")

    halfLife = math.log(2) / ke
    return dict(ke=ke, vd=vd, Cl=Cl, realCp_used=realCp_used, halfLife=halfLife)

def suggest_tinf_safe(dose):
    # Max rate 10 mg/min => min infusion minutes = dose/10
    min_inf_min = dose / 10.0
    min_inf_hr = min_inf_min / 60.0
    rounded = ceiling_precise(min_inf_hr, 0.5)
    return min_inf_hr, rounded

def suggest_interval(tau, halfLife, override_text):
    if override_text.strip():
        try:
            return float(override_text)
        except Exception:
            return None
    if tau < halfLife:
        return round(halfLife, 1)
    elif tau < 2 * halfLife:
        return round(2 * halfLife, 1)
    else:
        return tau

def steady_state_profile_bolus(dose, tau, ke, vd, t_grid):
    # C_ss(t) = (dose/vd) * exp(-ke * t) / (1 - exp(-ke * tau)), repeating every tau
    acc = 1.0 / (1.0 - math.exp(-ke * tau))
    c = (dose / vd) * np.exp(-ke * (t_grid % tau)) * acc
    return c

def steady_state_profile_infusion(dose, tau, tinf, ke, Cl, t_grid):
    # Intermittent infusion with rate R0 = dose / tinf
    # Compute steady-state Cmax and Cmin, then piecewise within each cycle
    if tinf <= 0:
        tinf = 1e-9  # avoid zero division; will behave close to bolus for tiny tinf
    R0 = dose / tinf
    # At steady state:
    Cmax = (R0 / Cl) * (1 - math.exp(-ke * tinf)) / (1 - math.exp(-ke * tau))
    Cmin = Cmax * math.exp(-ke * (tau - tinf))

    t_mod = t_grid % tau
    c = np.zeros_like(t_grid, dtype=float)
    # During infusion: 0 <= t <= tinf
    mask_infusion = t_mod <= tinf
    # C(t) = Cmin * e^{-ke*t} + (R0/Cl) * (1 - e^{-ke t})
    c[mask_infusion] = Cmin * np.exp(-ke * t_mod[mask_infusion]) + (R0 / Cl) * (1 - np.exp(-ke * t_mod[mask_infusion]))
    # Post-infusion: tinf < t <= tau
    mask_post = ~mask_infusion
    c[mask_post] = Cmax * np.exp(-ke * (t_mod[mask_post] - tinf))
    return c, Cmin, Cmax

def model_predicted_Ctrough(dose, tau, tinf, ke, vd, Cl, model):
    if model == "Bolus":
        # Ct at time tau - 0.5 after dose (as in VBA), using steady-state trough formula
        num = (dose / vd) * math.exp(-ke * (tau - 0.5))
        den = (1 - math.exp(-ke * tau))
        return num / den
    else:
        # Ct measured at 0.5 hr before next dose => time = tau - 0.5
        t = tau - 0.5
        if t <= tinf:
            # If 0.5 hr before end of cycle falls within infusion window (rare), use during-infusion formula
            R0 = dose / max(tinf, 1e-9)
            # Need Cmin
            _, Cmin, Cmax = steady_state_profile_infusion(dose, tau, tinf, ke, Cl, np.array([0.0,]))
            return float(Cmin * math.exp(-ke * t) + (R0 / Cl) * (1 - math.exp(-ke * t)))
        else:
            # Post-infusion decay from Cmax
            _, Cmin, Cmax = steady_state_profile_infusion(dose, tau, tinf, ke, Cl, np.array([0.0,]))
            return float(Cmax * math.exp(-ke * (t - tinf)))

def predict_Cpeak_adjusted(newDose, suggestedTau, tinf, ke, vd, Cl, model):
    if model == "Bolus":
        # VBA uses sampling at 0.5 hr after dose
        return (newDose / vd) * math.exp(-ke * 0.5) / (1 - math.exp(-ke * suggestedTau))
    else:
        # During infusion at t = 0.5 hr (within infusion)
        R0 = newDose / max(tinf, 1e-9)
        # Need Cmin at steady state for the adjusted regimen
        Cmax = (R0 / Cl) * (1 - math.exp(-ke * tinf)) / (1 - math.exp(-ke * suggestedTau))
        Cmin = Cmax * math.exp(-ke * (suggestedTau - tinf))
        t = 0.5
        return Cmin * math.exp(-ke * t) + (R0 / Cl) * (1 - math.exp(-ke * t))

# --- Core calculation ---
errs = []
try:
    # Safety infusion time
    min_inf_hr, rounded_tinf = suggest_tinf_safe(dose)
    if modelUsed == "Infusion" and tinf < rounded_tinf:
        st.warning(f"Current infusion time ({tinf:.2f} hr) is too short for safe rate (≤ 10 mg/min). Suggested minimum: **{rounded_tinf:.2f} hr**")

    # Estimate PK parameters
    res = estimate_params(realCt, realCp, dose, tau, tinf, modelUsed, weight)
    ke = res["ke"]; vd = res["vd"]; Cl = res["Cl"]; realCp_used = res["realCp_used"]; halfLife = res["halfLife"]

    # Suggested interval
    suggestedTau = suggest_interval(tau, halfLife, custom_suggested_tau)
    if suggestedTau is None or suggestedTau <= 0:
        st.error("Invalid override for Suggested Interval. Please enter a positive number.")
        st.stop()

    # Model-predicted trough at suggestedTau (with current dose)
    modelCt2 = None
    if modelUsed == "Bolus":
        modelCt2 = (dose / vd) * math.exp(-ke * (suggestedTau - 0.5)) / (1 - math.exp(-ke * suggestedTau))
    else:
        modelCt2 = (dose * (1 - math.exp(-ke * tinf)) * math.exp(-ke * (suggestedTau - tinf - 0.5))) / (Cl * tinf * (1 - math.exp(-ke * suggestedTau)))

    if modelCt2 <= 0 or not np.isfinite(modelCt2):
        raise ValueError("Error calculating model-predicted trough (modelCt2).")

    # Adjusted dose to reach Ct_target
    newDose = dose * (Ct_target / modelCt2)
    predictedCpeak = predict_Cpeak_adjusted(newDose, suggestedTau, tinf, ke, vd, Cl, modelUsed)

    # Prepare outputs
    col1, col2, col3 = st.columns(3)
    with col1:
        st.subheader("Estimated PK")
        st.metric("Ke (/hr)", f"{ke:.4f}")
        st.metric("Half-life (hr)", f"{halfLife:.2f}")
    with col2:
        st.subheader("Disposition")
        st.metric("Vd (L)", f"{vd:.2f}")
        st.metric("Cl (L/hr)", f"{Cl:.2f}")
    with col3:
        st.subheader("Dosing")
        st.metric("Suggested Interval (hr)", f"{suggestedTau:.2f}")
        st.metric("Adjusted Dose (mg) → Ct = {:.1f}".format(Ct_target), f"{newDose:.0f}")
        st.metric("Predicted Adjusted Cpeak (mg/L)", f"{predictedCpeak:.2f}")

    st.info(f"Suggested **tinf** to keep ≤ 10 mg/min: **{rounded_tinf:.2f} hr** (raw minimum: {min_inf_hr:.2f} hr)")

    # --- Simulation & Plot ---
    t_end = max(suggestedTau, tau) * 3.0
    t_grid = np.linspace(tStart, t_end, 600)

    if modelUsed == "Bolus":
        c_current = steady_state_profile_bolus(dose, tau, ke, vd, t_grid)
        c_adjusted = steady_state_profile_bolus(newDose, suggestedTau, ke, vd, t_grid)
    else:
        c_current, cmin_cur, cmax_cur = steady_state_profile_infusion(dose, tau, tinf, ke, Cl, t_grid)
        c_adjusted, cmin_adj, cmax_adj = steady_state_profile_infusion(newDose, suggestedTau, tinf, ke, Cl, t_grid)

    fig1 = plt.figure()
    plt.plot(t_grid, c_current, label="Current regimen")
    plt.plot(t_grid, c_adjusted, label="Adjusted regimen")
    plt.axhline(Ct_target, linestyle="--", label=f"Target Ct = {Ct_target:.1f}")
    plt.xlabel("Time (hr)")
    plt.ylabel("Concentration (mg/L)")
    plt.title("Steady-state Concentration–Time Profile")
    plt.legend()
    st.pyplot(fig1)

    # Table of key values
    df = pd.DataFrame({
        "Parameter": ["Ke (/hr)", "Half-life (hr)", "Vd (L)", "Cl (L/hr)", "Suggested τ (hr)", "Adjusted Dose (mg)",
                      "Predicted Adjusted Cpeak (mg/L)", "Inferred/Used Cpeak (mg/L)"],
        "Value": [round(ke,4), round(halfLife,2), round(vd,2), round(Cl,2), round(suggestedTau,2), round(newDose,0),
                  round(predictedCpeak,2), round(realCp_used,2)]
    })
    st.subheader("Summary")
    st.dataframe(df, use_container_width=True)

except Exception as e:
    st.error(f"Error: {str(e)}")
