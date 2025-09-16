#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def canonicalize_param_label(label: str) -> str:
    """
    Map UI-facing labels (metric or imperial) to internal canonical names.
    Return one of:
        'L2', 'L1', 'grad', 'D', 'Tinj', 'k'
    """
    if not label:
        return ""

    lbl = label.strip().lower()

    # horizontal extent
    if "horizontal" in lbl:
        return "L2"
    # vertical extent
    if "vertical" in lbl:
        return "L1"
    # geothermal gradient
    if "gradient" in lbl:
        return "grad"
    # diameter / borehole diameter
    if "diameter" in lbl:
        return "D"
    # injection temperature
    if "injection" in lbl and "temp" in lbl:
        return "Tinj"
    # thermal conductivity / rock thermal conductivity
    if "conductivity" in lbl:
        return "k"

    # fallback: nothing recognized
    return ""
