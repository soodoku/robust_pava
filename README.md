## Robust PAVA

Isotonic regression generally produces few unique values because all violations are pooled. This property can be consequential if we want to use the outcome for anything other than ranking. To increase uniqueness, we can amend PAVA to ignore small violations. The trade-off is a slightly lower rank-ordered correlation. The repository contains code for the same along with a simulation that illustrates its impact.



