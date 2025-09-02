import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO

# Provided dataset as a string
data = """host.id RF_pred_linear cohort
phyloscanner_2023000916 0.330020136 0.265753425
phyloscanner_2023001299 0.061225545 0.249315068
phyloscanner_2023002072 2.406356125 0.224657534
phyloscanner_2023002418 2.884574939 0.21369863
phyloscanner_2023003610 0.084756606 0.22739726
phyloscanner_2023005569 0.084670333 0.252054795
phyloscanner_2023006138 0.031311851 0.205479452
phyloscanner_469 0.058212157 1.582191781
phyloscanner_577 3.036113238 2
phyloscanner_579 0.058716883 1.498630137
phyloscanner_GHM03156 1.301533279 0.194520548
phyloscanner_GHM03227 0.061112264 0.315068493
phyloscanner_GHM03828 0.050140459 0.257534247
phyloscanner_GHM04028 0.062948033 0.268493151
phyloscanner_GHM10101 0.246185186 0.284931507
phyloscanner_GHM11068 0.203587332 0.189041096
phyloscanner_GHM12290 0.891155607 0.104109589
phyloscanner_GHM13862 0.862900822 0.328767123
phyloscanner_GHM14754 0.567350006 0.216438356
phyloscanner_GHM14891 0.685489858 0.136986301
phyloscanner_GHM15575 0.712061481 0.273972603
phyloscanner_OOB12152913 0.697757089 0.128767123
phyloscanner_OOB12152949 0.552119949 0.224657534
phyloscanner_OOC194444 0.170198027 0.136986301
phyloscanner_OOC194532 0.242261978 0.161643836
phyloscanner_OOC270001 0.318178038 0.161643836
phyloscanner_OOC270535 0.057972835 0.2
"""

# Prediction intervals
pred_intervals = """RF_pred_min_linear RF_pred_max_linear
0 3.780453206
0 1.615188033
0.311603335 6.473322089
0.265923723 8.300906408
0 1.7548661
0 1.694608543
0 0.324535882
0 1.643951031
0.042174437 10.75528471
0 1.646041615
0.121286242 3.738165042
0 1.699681642
0 1.700114452
0 1.69484505
0.003129597 0.876841666
0.00199742 0.73568455
0.290720427 1.819358029
0 3.772428703
0.141854425 1.276486745
0.234573813 1.372548358
0.259333986 1.388688443
0.230586731 1.417153676
0.140525767 1.23482727
0 2.384736715
0 3.356017422
0 2.393753443
0 1.652617112
"""

# Read into pandas
df = pd.read_csv(StringIO(data), sep=" ")
intervals = pd.read_csv(StringIO(pred_intervals), sep=" ")

# Merge
df = pd.concat([df, intervals], axis=1)

# Plot
plt.figure(figsize=(12, 6))

# Main lines
plt.plot(df["host.id"], df["RF_pred_linear"], marker='o', markersize=3, label="RF_pred_linear")
plt.plot(df["host.id"], df["cohort"], marker='s', markersize=3, label="cohort")

# Prediction interval shading
plt.fill_between(df["host.id"], df["RF_pred_min_linear"], df["RF_pred_max_linear"],
                  alpha=0.2, label="Prediction Interval")

plt.xticks(rotation=90)
plt.xlabel("Host ID")
plt.ylabel("Estimated Duration (years)")
plt.title("Comparison of HIV-Phylotsi and UVRI TSI estimates")
plt.legend()
plt.tight_layout()
plt.show()
