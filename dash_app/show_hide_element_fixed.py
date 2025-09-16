def show_hide_element(visibility_state, tab, fluid, end_use, model):

    # ----------------------------------------------------------------------------------------------
    # Reveals or hides sliders depending on which tab selected and which dropdowns.
    # ----------------------------------------------------------------------------------------------

    # print("show_hide_element: ", model, tab, fluid, end_use, visibility_state)
    # fluid getting changed makes this run "twice" but it's working as expected.
    
    try:
        from utils_labels import canonicalize_param_label
        
        b = {'display': 'block'}
        n = {'display': 'none'}
        
        # Canonicalize the parameter label to handle both metric and imperial units
        param_key = canonicalize_param_label(visibility_state) if visibility_state else ""
        
        def style(show):  # helper
            return b if show else n

        if model == "HDF5":
            if tab == "about-tab":
                if fluid == "H2O" or end_use == "Heating":
                    return b, b, b, b, b, b, b,\
                            b, b, b, n, n, b

                else:
                    return b, b, b, b, b, b, b,\
                            b, b, b, b, b, b
            
            elif tab == "energy-time-tab":
                return b, b, b, b, b, b, b, \
                        n, n, n, n, n, b
            
            elif tab == "energy-tab":
                # Use canonicalized parameter keys instead of hardcoded param_list indices
                show_mdot = (param_key == "mdot")
                show_L2 = (param_key == "L2")
                show_L1 = (param_key == "L1")
                show_grad = (param_key == "grad")
                show_D = (param_key == "D")
                show_Tinj = (param_key == "Tinj")
                show_k = (param_key == "k")
                
                return (
                    style(show_mdot),  # mdot-select-div
                    style(show_L2),    # L2-select-div
                    style(show_L1),    # L1-select-div
                    style(show_grad),  # grad-select-div
                    style(show_D),     # diameter-select-div
                    style(show_Tinj),  # Tinj-select-div
                    style(show_k),     # k-select-div
                    n, n, n, n, n, b   # drillcost, discount-rate, lifetime, precool, turb-pout, wellbore-params
                )
            
            elif tab == "economics-time-tab":
                if fluid == "H2O":
                    if end_use == "All":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b
                    if end_use == "Heating":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b
                    if end_use == "Electricity":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b

                else:
                    if end_use == "All":
                        return b, b, b, b, b, b, b, \
                                b, b, b, b, b, b
                    if end_use == "Heating":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b
                    if end_use == "Electricity":
                        return b, b, b, b, b, b, b, \
                                b, b, b, b, b, b

            elif tab == "summary-tab":
                if fluid == "H2O":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b
                else:
                    return b, b, b, b, b, b, b, \
                            b, b, b, b, b, b
                        
        elif model == "SBT V1.0":
            if tab == "about-tab":
                if fluid == "H2O" or end_use == "Heating":
                    return b, b, b, b, b, b, b,\
                            b, b, b, n, n, b

                else:
                    return b, b, b, b, b, b, b,\
                            b, b, b, b, b, b
            
            elif tab == "energy-time-tab":
                return b, b, b, b, b, b, b, \
                        n, n, n, n, n, b
            
            elif tab == "energy-tab":
                # Use canonicalized parameter keys instead of hardcoded param_list indices
                show_mdot = (param_key == "mdot")
                show_L2 = (param_key == "L2")
                show_L1 = (param_key == "L1")
                show_grad = (param_key == "grad")
                show_D = (param_key == "D")
                show_Tinj = (param_key == "Tinj")
                show_k = (param_key == "k")
                
                return (
                    style(show_mdot),  # mdot-select-div
                    style(show_L2),    # L2-select-div
                    style(show_L1),    # L1-select-div
                    style(show_grad),  # grad-select-div
                    style(show_D),     # diameter-select-div
                    style(show_Tinj),  # Tinj-select-div
                    style(show_k),     # k-select-div
                    n, n, n, n, n, b   # drillcost, discount-rate, lifetime, precool, turb-pout, wellbore-params
                )
            
            elif tab == "economics-time-tab":
                if fluid == "H2O":
                    if end_use == "All":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b
                    if end_use == "Heating":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b
                    if end_use == "Electricity":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b

                else:
                    if end_use == "All":
                        return b, b, b, b, b, b, b, \
                                b, b, b, b, b, b
                    if end_use == "Heating":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b
                    if end_use == "Electricity":
                        return b, b, b, b, b, b, b, \
                                b, b, b, b, b, b

            elif tab == "summary-tab":
                if fluid == "H2O":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b
                else:
                    return b, b, b, b, b, b, b, \
                            b, b, b, b, b, b
                        
        elif model == "SBT V2.0":
            if tab == "about-tab":
                if fluid == "H2O" or end_use == "Heating":
                    return b, b, b, b, b, b, b,\
                            b, b, b, n, n, b
                else:
                    return b, b, b, b, b, b, b,\
                            b, b, b, b, b, b
            
            elif tab == "energy-time-tab":
                return b, b, b, b, b, b, b, \
                        n, n, n, n, n, b
            
            elif tab == "energy-tab":
                # Use canonicalized parameter keys instead of hardcoded param_list indices
                show_mdot = (param_key == "mdot")
                show_L2 = (param_key == "L2")
                show_L1 = (param_key == "L1")
                show_grad = (param_key == "grad")
                show_D = (param_key == "D")
                show_Tinj = (param_key == "Tinj")
                show_k = (param_key == "k")
                
                return (
                    style(show_mdot),  # mdot-select-div
                    style(show_L2),    # L2-select-div
                    style(show_L1),    # L1-select-div
                    style(show_grad),  # grad-select-div
                    style(show_D),     # diameter-select-div
                    style(show_Tinj),  # Tinj-select-div
                    style(show_k),     # k-select-div
                    n, n, n, n, n, b   # drillcost, discount-rate, lifetime, precool, turb-pout, wellbore-params
                )
            
            elif tab == "economics-time-tab":
                if fluid == "H2O":
                    if end_use == "All":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b
                    if end_use == "Heating":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b
                    if end_use == "Electricity":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b
                else:
                    if end_use == "All":
                        return b, b, b, b, b, b, b, \
                                b, b, b, b, b, b
                    if end_use == "Heating":
                        return b, b, b, b, b, b, b, \
                                b, b, b, n, n, b
                    if end_use == "Electricity":
                        return b, b, b, b, b, b, b, \
                                b, b, b, b, b, b

            elif tab == "summary-tab":
                if fluid == "H2O":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b
                else:
                    return b, b, b, b, b, b, b, \
                            b, b, b, b, b, b
        else:
            raise PreventUpdate
        
    except Exception as e:
        print("Style callback failed:", e)
        # Return same number of outputs with safe defaults
        N = 13  # match your outputs count
        return tuple({'display': 'block'} for _ in range(N))
