"""
Global settings
:global s_base: The system complex power value to use as a basis for the power system in per unit.
:type s_base: int
:global mm_max: The maximum amount of mismatch on power convergence in power analysis.
:type mm_max: float
:global v_max: A bus maximum voltage limit to check for compliance.
:type v_max: float
:global v_min: A bus minimum voltage limit to check for compliance.
:type v_min: float
"""
# Complex Power Base
s_base = 100 # MVA

# Maximum Mismatch
mm_max = 0.1 / s_base

# Voltage Limits
v_max = 1.05 # p.u.
v_min = 0.95 # p.u.