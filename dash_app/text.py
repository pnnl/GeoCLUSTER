#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------------
# Contains text content for dash app.
# -----------------------------------------------------------------------------


# Left Panel Text
scenario1_text = "Evaluates thermal and economic performance of a commercial-scale operation, \
                                which is representative of hot-dry-rock (HDR) geothermal reservoirs within the continental United States."

scenario2_text = "Calculates HDR results for operations that optimize power output."

scenario3_text = "Calculates near-minima levelized-cost-of-heat (LCOH) and levelized-cost-of-electricity (LCOE) results. \
                            For supercritical CO2 (sCO2), multiple LCOE minima exist."

disclaimer_text = "The information on this website is for general informational purposes only. \
                    Your use of the site is solely at your own risk."


# Graph Guidance Text

dropdown_text1 = "Hover over the legend to interact with the pop-up graphics panel. There are 5 options to choose from. \
                        To download a static PNG, click on the leftmost option. To zoom in, drag \
                        the curser across a graph."

dropdown_text2 = "Contours show CLGS HDR results where parameters can be maximized. \
                        Contour lines and color indicate elevation or depression where maximation occurs \
                        in areas with lighter colors and minimization in areas with darker colors. To visualize \
                        various contours, select a parameter that will be compared to mass flow rate. Then, edit \
                        remaining parameters in the 'Finetune System Values' pane to see how the contours change."

dropdown_econ_text1 = "To address the economic aspects of CLGSs, levelized cost of heating (LCOH) and levelized cost of electricity (LCOE), \
                        which represent net average present costs, were calculated as follows:"

dropdown_econ_markdown_text1 = """where, **_$C_t$_** is capital expenditures (CAPEX) in year **_t_**, \
                            **_$O_t$_** is operating and maintenance expenditures (OPEX) in year **_t_**, \
                            **_n_** is expected lifetime of the CLGS, \
                            **_r_** is discount rate, \
                            **_$E_h$_** is heat produced in year **_t_**, \
                            and **_$E_e$_** is electricity generated in year **_t_**."""

dropdown_econ_text2 = "Finally, different assumptions were defined based on the working fluid. \
                            When water (H2O) is the heat transfer fluid, an organic Rankine cycle was assumed. \
                            When supercritical carbon dioxide (sCO2) is the circulating fluid, a direct turbine expansion cycle was assumed \
                            where the circulating CO2 is also the working fluid through the cycle."

note = "This work explores closed-loop geothermal systems (CLGSs). These deep geothermal systems are \
designed to extract thermal energy at the commercial scale for applications, \
such as district heating and electricity generation. \
This research was funded by the Geothermal Technologies Office (GTO) within the Office of Energy Efficiency \
and Renewable Energy (EERE) at the U.S. Department of Energy (DOE) to form a collaborative study of CLGSs involving four national \
laboratories and two universities."


note2 = "CLGSs extract thermal energy from hot-dry-rock (HDR) by circulating a working fluid through \
a system of sealed boreholes originating and terminating \
at the ground surface and passing through HDR. Two basic borehole configurations exist: \
1) u-shaped and 2) co-axial. In the u-shaped \
configuration, the borehole inlet and outlet collars are separate while, in the coaxial configuration, \
the inlet and outlet fluids are separated within a single collar via co-axial tubing. For HDR, thermal \
conduction is the sole heat transport mechanism by which energy is transferred from the rock to the \
borehole exterior."

note3 = "Explore thermal performance and economic analyses across multiple parameters for a hot-dry-rock reservoir. \
All parameters are pinned under 'Finetune System Values' and several parameters are demonstrated in the visuals. \
To learn more about HDR subsurface and economic results, click on the tabs above."
