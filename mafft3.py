import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.stats.power as sm_power

# Read the data from the CSV file
data = pd.read_csv("~/Desktop/tf_trial/analysis/analysis2/combined_sequences.csv")

# Separate the data for MinION 1 and MinION 2
minion1_data = data[data["MinION1"] == "MinION1"]
minion2_data = data[data["MinION2"] == "MinION2"]


# Calculate the error rates and pairwise identities for MinION 1
minion1_error_rates = (minion1_data["Mismatches"] + minion1_data["Indels"]) / minion1_data["Sequence Length"] * 100
minion1_pairwise_identities = (minion1_data["Sequence Length"] - minion1_data["Mismatches"] - minion1_data["Indels"]) / minion1_data["Sequence Length"] * 100

# Calculate the error rates and pairwise identities for MinION 2
minion1_error_rates = (minion1_data["Mismatches"] + minion1_data["Indels"]) / minion1_data["Sequence Length"] * 100
minion1_pairwise_identities = (minion1_data["Sequence Length"] - minion1_data["Mismatches"] - minion1_data["Indels"]) / minion1_data["Sequence Length"] * 100

minion2_error_rates = (minion2_data["Mismatches"] + minion2_data["Indels"]) / minion2_data["Sequence Length"] * 100
minion2_pairwise_identities = (minion2_data["Sequence Length"] - minion2_data["Mismatches"] - minion2_data["Indels"]) / minion2_data["Sequence Length"] * 100

# Perform t-tests for error rates and pairwise identities
t_statistic_error_rates, p_value_error_rates = stats.ttest_ind(minion1_error_rates, minion2_error_rates)
t_statistic_pairwise_identities, p_value_pairwise_identities = stats.ttest_ind(minion1_pairwise_identities, minion2_pairwise_identities)

# Calculate the confidence intervals for error rates and pairwise identities
confidence_interval_minion1_error_rates = stats.t.interval(0.95, len(minion1_error_rates) - 1, loc=np.mean(minion1_error_rates), scale=stats.sem(minion1_error_rates))
confidence_interval_minion2_error_rates = stats.t.interval(0.95, len(minion2_error_rates) - 1, loc=np.mean(minion2_error_rates), scale=stats.sem(minion2_error_rates))
confidence_interval_minion1_pairwise_identities = stats.t.interval(0.95, len(minion1_pairwise_identities) - 1, loc=np.mean(minion1_pairwise_identities), scale=stats.sem(minion1_pairwise_identities))
confidence_interval_minion2_pairwise_identities = stats.t.interval(0.95, len(minion2_pairwise_identities) - 1, loc=np.mean(minion2_pairwise_identities), scale=stats.sem(minion2_pairwise_identities))

# Perform power analysis for error rates
effect_size_error_rates = np.abs(np.mean(minion1_error_rates) - np.mean(minion2_error_rates)) / np.std(minion1_error_rates + minion2_error_rates)
power_error_rates = sm_power.tt_ind_solve_power(effect_size=effect_size_error_rates, nobs1=len(minion1_error_rates), alpha=0.05, power=None, ratio=len(minion2_error_rates) / len(minion1_error_rates), alternative="two-sided")

# Perform power analysis for pairwise identities
effect_size_pairwise_identities = np.abs(np.mean(minion1_pairwise_identities) - np.mean(minion2_pairwise_identities)) / np.std(minion1_pairwise_identities + minion2_pairwise_identities)
power_pairwise_identities = sm_power.tt_ind_solve_power(effect_size=effect_size_pairwise_identities, nobs1=len(minion1_pairwise_identities), alpha=0.05, power=None, ratio=len(minion2_pairwise_identities) / len(minion1_pairwise_identities), alternative="two-sided")

# Print the results
print("Comparative Analysis:")
print(f"t-test for Error Rates: t-statistic = {t_statistic_error_rates}, p-value = {p_value_error_rates}")
print(f"t-test for Pairwise Identities: t-statistic = {t_statistic_pairwise_identities}, p-value = {p_value_pairwise_identities}")
print()

print("Confidence Intervals:")
print(f"MinION 1 Error Rate Confidence Interval: {confidence_interval_minion1_error_rates}")
print(f"MinION 2 Error Rate Confidence Interval: {confidence_interval_minion2_error_rates}")
print(f"MinION 1 Pairwise Identity Confidence Interval: {confidence_interval_minion1_pairwise_identities}")
print(f"MinION 2 Pairwise Identity Confidence Interval: {confidence_interval_minion2_pairwise_identities}")
print()

print("Power Analysis:")
print(f"Power for Error Rates: {power_error_rates}")
print(f"Power for Pairwise Identities: {power_pairwise_identities}")
print()

# Visualization
sns.set(style="whitegrid")

plt.figure(figsize=(10, 6))

plt.subplot(2, 2, 1)
sns.boxplot(data=data, x="Sample_id", y="MinION1", palette="Set3")
plt.xlabel("MinION Dataset")
plt.ylabel("Error Rate (%)")
plt.title("Error Rates Comparison")

plt.subplot(2, 2, 2)
sns.violinplot(data=data, x="Sample_id", y="MinION1", palette="Set3")
plt.xlabel("MinION Dataset")
plt.ylabel("Error Rate (%)")
plt.title("Error Rates Comparison")

plt.subplot(2, 2, 3)
sns.boxplot(data=data, x="Sample_id", y="MinION2", palette="Set3")
plt.xlabel("MinION Dataset")
plt.ylabel("Pairwise Identity (%)")
plt.title("Pairwise Identities Comparison")

plt.subplot(2, 2, 4)
sns.violinplot(data=data, x="Sample_id", y="MinION2", palette="Set3")
plt.xlabel("MinION Dataset")
plt.ylabel("Pairwise Identity (%)")
plt.title("Pairwise Identities Comparison")

plt.tight_layout()
plt.show()
